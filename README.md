# structural-evolution

This repository scripts for running the analysis described in the paper ["Unsupervised evolution of antibody and protein complexes with a structure-informed language model"](https://www.science.org/doi/10.1126/science.adk8946).

<p align="center">
  <img src="./structural-evolution-overview.png" alt="structural-evolution-overview"  width="700px"/>
</p>

## Setup/Installation
1. Clone this repository
```
git clone https://github.com/varun-shanker/structural-evolution.git

```
2. Install Conda Environment and Required Dependencies
```

```
3. Download and unzip the model weights from ["here"](https://zenodo.org/records/12631662), then insert them in torch checkpoints.
```
wget -P ~/.cache/torch/hub/checkpoints https://zenodo.org/records/12631662/files/esm_if1_20220410.zip
unzip ~/.cache/torch/hub/checkpoints/esm_if1_20220410.zip
```
4. Navigate to repository
```
cd structural-evolution
```

## Generating Predictions

To evaluate this model on a new protein or protein complex structure, run
```bash
python bin/recommend.py [pdb/cif file] --chain [X]
```
where `[pdb file]` is the file path to the pdb/cif structure file of the protein or protein complex and `[X]` is the target chain you wish to evolve. The default script will output the top `n`=10 predicted substitutions at unique residue positions (`maxrep=1`), where `n` and `maxrep` can be modified using the arguments (see below).

To recommend mutations to antibody variable domain sequences, we have simply run the above script separately on the heavy and light chains.

Additional arguments:

```
--seqpath: filepath where fasta with dms library should be saved (defaults to new subdirectory in outputs directory)
--outpath: output filepath for scores of variant sequences (defaults to new subdirectory in outputs directory)
--chain: chain id for the chain of interest
--n: number of top recommendations to be output (default: n=10)
--maxrep: maximum representation of a single site in the output recommendations (default: maxrep = 1 is a unique set of recommendations where each mutation of a given wildtype residue is recommended at most once)
--upperbound: only residue positions less than the user-defined upperbound are considered for recommendation in the final output (but all positions are still conditioned for scoring)
--order: for multichain conditioning, provides option to specify the order of chains
--offset: integer offset or adjustment for labeling of residue indices encoded in the structure file
--multichain-backbone: use the backbones of all chains in the input for conditioning (default is True)
--singlechain-backbone: use the backbone of only the target chain in the input for conditioning
--nogpu: Do not use GPU even if available
```

## Paper analysis scripts

To reproduce the analysis in the paper, first download and extract data with the commands:
```bash
wget https://zenodo.org/record/6968342/files/data.tar.gz
tar xvf data.tar.gz
```

## Citation

Please cite the following publication when referencing this work.

```
@article {Shanker-struct-evo,
	author = {Shanker, Varun and Bruun, Theodora and Hie, Brian and Kim, Peter},
	title = {Unsupervised evolution of antibody and protein complexes with a structure-informed language model},
	year = {2024},
	doi = {10.1126/science.adk8946},
	publisher = {American Association for the Advancement of Science},
	URL = {https://www.science.org/doi/10.1126/science.adk8946},
	journal = {Science}
}
```

## License
This project is licensed under the MIT License - see the LICENSE file for details.
