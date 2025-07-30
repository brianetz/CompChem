# CompChem

This repository contains a computational chemistry workflow for calculating bond dissociation energies (BDEs), built in Python and using NWChem as the electronic structure engine. The workflow automates the process of preparing, running, and analyzing quantum chemistry calculations for both a parent molecule and its molecular fragments. 

This workflow is detailed in [publication coming soon].

The workflow starts from a user-provided SMILES string and proceeds through the following key stages:

1. Conformer Search to identify the lowest-energy molecular conformation using conformer generation and energy ranking tools.

2. Bond Fragmentation to break every single bond in the molecule to generate fragment pairs representing potential BDE pathways.

3. DFT Input Preparation to generate NWChem input files for both the parent molecule and each fragment.

4. Job Submission to submit geometry optimizations, single-point energy calculations, and vibrational frequency analyses to an HPC cluster or local environment.

5. Output Parsing to extract final energies, thermodynamic corrections (ZPE, enthalpy, entropy), and frequencies from NWChem output files.

6. Post-Processing to compute bond dissociation enthalpies (BDEs) and generate structured data products for further analysis or visualization.
