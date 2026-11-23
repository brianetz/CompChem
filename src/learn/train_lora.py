"""LoRA fine-tuning of the plan-board Planner on its own corrected mistakes.

Stage 1, SFT: imitate accepted, oracle-confirmed module attempts.
Stage 2, DPO: prefer the correct attempt over the planner's mistake for the same
first-try prompt; the reference model is the base (LoRA disabled). A small SFT
term on the chosen reply keeps the policy from drifting.

Loss covers only the assistant turn (reasoning + tool call), never the prompt.
The model is spread over the node's GCDs (device_map="auto"); LoRA sits on the
attention and Mamba projections.

usage: python -m learn.train_lora --base <hf dir> --sft sft.jsonl --dpo dpo.jsonl --out <adapter dir>
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import time

import torch
import torch.nn.functional as F
from peft import LoraConfig, get_peft_model
from transformers import AutoModelForCausalLM, AutoTokenizer

TARGETS = ["q_proj", "k_proj", "v_proj", "o_proj", "in_proj", "out_proj"]


def load_jsonl(path):
    return [json.loads(l) for l in open(path) if l.strip()] if path and os.path.exists(path) else []


def render(tok, ex: dict, reply: dict, max_len: int):
    """Token ids for prompt + reply, and the index where the reply starts."""
    if "text" in reply:   # a format failure: the model wrote text instead of calling the tool
        turn = {"role": "assistant", "content": reply["text"], "reasoning_content": reply.get("thinking") or ""}
    else:
        turn = {"role": "assistant", "content": "", "reasoning_content": reply.get("thinking") or "",
                "tool_calls": [{"type": "function", "function": {"name": reply["name"], "arguments": reply["arguments"]}}]}
    msgs = ex["messages"] + [turn]
    full = tok.apply_chat_template(msgs, tools=ex["tools"], tokenize=False)
    prompt = tok.apply_chat_template(ex["messages"], tools=ex["tools"], add_generation_prompt=True, tokenize=False)
    if not full.startswith(prompt):  # template differences: fall back to the common prefix
        n = 0
        while n < min(len(full), len(prompt)) and full[n] == prompt[n]:
            n += 1
        prompt = full[:n]
    p_ids = tok(prompt, add_special_tokens=False)["input_ids"]
    f_ids = tok(full, add_special_tokens=False)["input_ids"]
    if len(f_ids) > max_len:
        cut = len(f_ids) - max_len
        f_ids, p_len = f_ids[cut:], max(0, len(p_ids) - cut)
    else:
        p_len = len(p_ids)
    return torch.tensor([f_ids]), p_len


def seq_logp(model, ids, p_len, device):
    ids = ids.to(device)
    start = max(0, p_len - 1)
    # only the reply's positions are scored: convert just that slice to float32
    logits = model(input_ids=ids).logits[:, start:-1]
    tgt = ids[:, start + 1:].to(logits.device)
    lp = torch.log_softmax(logits.float(), -1).gather(-1, tgt.unsqueeze(-1)).squeeze(-1)
    return lp.sum(), torch.tensor(lp.numel(), device=lp.device)


class Log:
    def __init__(self, path):
        self.fh = open(path, "a")

    def __call__(self, **rec):
        rec["t"] = round(time.time(), 1)
        self.fh.write(json.dumps(rec) + "\n")
        self.fh.flush()
        print(rec, flush=True)


def lr_lambda(total, warm=0.1):
    def f(step):
        w = max(1, int(total * warm))
        return (step + 1) / w if step < w else 0.5 * (1 + math.cos(math.pi * (step - w) / max(1, total - w)))
    return f


def run_sft(model, tok, rows, a, dev, log):
    random.Random(0).shuffle(rows)
    nval = max(1, len(rows) // 10)
    val, train = rows[:nval], rows[nval:]
    params = [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(params, lr=a.lr, weight_decay=0.0)
    steps = math.ceil(len(train) * a.sft_epochs / a.accum)
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lr_lambda(steps))
    k = 0
    for ep in range(a.sft_epochs):
        model.train()
        for i, ex in enumerate(train):
            ids, p_len = render(tok, ex, ex["reply"], a.max_len)
            lp, n = seq_logp(model, ids, p_len, dev)
            (-(lp / n) / a.accum).backward()
            if (i + 1) % a.accum == 0 or i == len(train) - 1:
                torch.nn.utils.clip_grad_norm_(params, 1.0)
                opt.step(); sched.step(); opt.zero_grad(set_to_none=True); k += 1
                log(stage="sft", epoch=ep, step=k, loss=round(float(-(lp / n)), 4), lr=sched.get_last_lr()[0])
        model.eval()
        with torch.no_grad():
            vl = [float(-(lambda r: r[0] / r[1])(seq_logp(model, *render(tok, ex, ex["reply"], a.max_len), dev))) for ex in val]
        log(stage="sft_val", epoch=ep, val_loss=round(sum(vl) / len(vl), 4), n=len(vl))


def run_dpo(model, tok, rows, a, dev, log):
    random.Random(1).shuffle(rows)
    nval = max(1, len(rows) // 10)
    val, train = rows[:nval], rows[nval:]
    params = [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(params, lr=a.dpo_lr, weight_decay=0.0)
    steps = math.ceil(len(train) * a.dpo_epochs / a.accum)
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lr_lambda(steps))

    def pair(ex, train_mode):
        c_ids, c_p = render(tok, ex, ex["chosen"], a.max_len)
        r_ids, r_p = render(tok, ex, ex["rejected"], a.max_len)
        with torch.no_grad(), model.disable_adapter():
            rc, _ = seq_logp(model, c_ids, c_p, dev)
            rr, _ = seq_logp(model, r_ids, r_p, dev)
        ctx = torch.enable_grad() if train_mode else torch.no_grad()
        with ctx:
            pc, nc = seq_logp(model, c_ids, c_p, dev)
            pr, _ = seq_logp(model, r_ids, r_p, dev)
            margin = a.beta * ((pc - rc) - (pr - rr))
            loss = -F.logsigmoid(margin) + a.sft_alpha * (-(pc / nc))
        return loss, float(margin)

    k = 0
    for ep in range(a.dpo_epochs):
        model.train()
        for i, ex in enumerate(train):
            loss, m = pair(ex, True)
            (loss / a.accum).backward()
            if (i + 1) % a.accum == 0 or i == len(train) - 1:
                torch.nn.utils.clip_grad_norm_(params, 1.0)
                opt.step(); sched.step(); opt.zero_grad(set_to_none=True); k += 1
                log(stage="dpo", epoch=ep, step=k, loss=round(float(loss), 4), margin=round(m, 3))
        model.eval()
        ms = [pair(ex, False)[1] for ex in val]
        log(stage="dpo_val", epoch=ep, reward_acc=round(sum(x > 0 for x in ms) / len(ms), 3),
            mean_margin=round(sum(ms) / len(ms), 3), n=len(ms))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True)
    ap.add_argument("--sft"), ap.add_argument("--dpo")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-len", type=int, default=7168)
    ap.add_argument("--rank", type=int, default=16)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--dpo-lr", type=float, default=2e-5)
    ap.add_argument("--sft-epochs", type=int, default=1)
    ap.add_argument("--dpo-epochs", type=int, default=1)
    ap.add_argument("--accum", type=int, default=8)
    ap.add_argument("--beta", type=float, default=0.1)
    ap.add_argument("--sft-alpha", type=float, default=0.2)
    ap.add_argument("--max-gib", type=int, default=52)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    log = Log(os.path.join(a.out, "train_log.jsonl"))
    tok = AutoTokenizer.from_pretrained(a.base, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        a.base, trust_remote_code=True, torch_dtype=torch.bfloat16, device_map="auto",
        max_memory={i: f"{a.max_gib}GiB" for i in range(torch.cuda.device_count())})
    model.gradient_checkpointing_enable(gradient_checkpointing_kwargs={"use_reentrant": False})
    model.enable_input_require_grads()
    model.config.use_cache = False
    model = get_peft_model(model, LoraConfig(r=a.rank, lora_alpha=2 * a.rank, lora_dropout=0.05,
                                             target_modules=TARGETS, task_type="CAUSAL_LM"))
    dev = model.get_input_embeddings().weight.device
    sft, dpo = load_jsonl(a.sft), load_jsonl(a.dpo)

    def fits(ex, key):   # skip, rather than cut, examples longer than --max-len
        return render(tok, ex, ex[key], 10 ** 6)[0].shape[1] <= a.max_len
    n_sft, n_dpo = len(sft), len(dpo)
    sft = [r for r in sft if fits(r, "reply")]
    dpo = [r for r in dpo if fits(r, "chosen") and fits(r, "rejected")]
    log(stage="start", sft=len(sft), dpo=len(dpo), skipped_too_long={"sft": n_sft - len(sft), "dpo": n_dpo - len(dpo)},
        args=vars(a))
    if sft:
        run_sft(model, tok, sft, a, dev, log)
        model.save_pretrained(os.path.join(a.out, "after_sft"))
        log(stage="saved_sft", out=os.path.join(a.out, "after_sft"))
    if dpo:
        run_dpo(model, tok, dpo, a, dev, log)
    model.save_pretrained(a.out)
    log(stage="saved", out=a.out)


if __name__ == "__main__":
    main()
