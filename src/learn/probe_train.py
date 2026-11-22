"""Feasibility probe: can nemotron-3-nano (bf16 HF weights) be LoRA-trained on one
Frontier node? Loads the model over the node's GCDs, adds LoRA to the attention
and Mamba projections, runs forward + backward + one optimizer step on a board-
sized sequence, saves the adapter and generates a few tokens.

usage: python -m learn.probe_train --base <hf dir> --out <adapter dir> [--tokens 3000]
"""
from __future__ import annotations

import argparse
import json
import time

import torch
from peft import LoraConfig, get_peft_model
from transformers import AutoModelForCausalLM, AutoTokenizer

TARGETS = ["q_proj", "k_proj", "v_proj", "o_proj", "in_proj", "out_proj"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tokens", type=int, default=3000)
    ap.add_argument("--max-gib", type=int, default=52)
    a = ap.parse_args()
    report = {"devices": torch.cuda.device_count(), "torch": torch.__version__, "hip": torch.version.hip}
    t = time.time()
    tok = AutoTokenizer.from_pretrained(a.base, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        a.base, trust_remote_code=True, torch_dtype=torch.bfloat16, device_map="auto",
        max_memory={i: f"{a.max_gib}GiB" for i in range(torch.cuda.device_count())})
    report["load_s"] = round(time.time() - t, 1)
    report["devices_used"] = sorted({str(p.device) for p in model.parameters()})
    model.gradient_checkpointing_enable(gradient_checkpointing_kwargs={"use_reentrant": False})
    model.enable_input_require_grads()
    model.config.use_cache = False
    cfg = LoraConfig(r=16, lora_alpha=32, lora_dropout=0.0, target_modules=TARGETS, task_type="CAUSAL_LM")
    model = get_peft_model(model, cfg)
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total = sum(p.numel() for p in model.parameters())
    report["trainable_params"], report["total_params"] = trainable, total
    msgs = [{"role": "system", "content": "You plan bond dissociation energy calculations. " * 40},
            {"role": "user", "content": ("Board state: target accepted. Publish the 'fragments' module now. " * 60)},
            {"role": "assistant", "content": json.dumps({"fragments": [{"smiles": "[CH2]c1ccccc1", "charge": 0,
                                                                          "multiplicity": 2},
                                                                         {"smiles": "[H]", "charge": 0,
                                                                          "multiplicity": 2}]})}]
    ids = tok.apply_chat_template(msgs, tokenize=True, return_tensors="pt")
    ids = ids[:, : a.tokens]
    report["seq_len"] = int(ids.shape[1])
    dev = model.get_input_embeddings().weight.device
    ids = ids.to(dev)
    labels = ids.clone()
    labels[:, : ids.shape[1] - 60] = -100       # loss only on the assistant tail
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=1e-4)
    model.train()
    for step in range(2):
        torch.cuda.synchronize()
        t = time.time()
        out = model(input_ids=ids, labels=labels)
        out.loss.backward()
        opt.step()
        opt.zero_grad(set_to_none=True)
        torch.cuda.synchronize()
        report[f"step{step}_s"] = round(time.time() - t, 1)
        report[f"step{step}_loss"] = round(float(out.loss), 4)
    report["max_mem_gib"] = [round(torch.cuda.max_memory_allocated(i) / 2**30, 1) for i in range(torch.cuda.device_count())]
    model.save_pretrained(a.out)
    model.eval()
    prompt = tok.apply_chat_template([{"role": "user", "content": "In one word, what is the SMILES of methane?"}],
                                     add_generation_prompt=True, return_tensors="pt").to(dev)
    with torch.no_grad():
        gen = model.generate(input_ids=prompt, max_new_tokens=24, do_sample=False)
    report["generation"] = tok.decode(gen[0, prompt.shape[1]:], skip_special_tokens=True)
    print("PROBE", json.dumps(report))


if __name__ == "__main__":
    main()
