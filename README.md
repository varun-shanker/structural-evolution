# structural-evolution

This repository scripts for running the analysis described in the paper ["Unsupervised evolution of antibody and protein complexes with a structure-informed language model"](https://www.science.org/doi/10.1126/science.adk8946).
<img width="875" alt="structural-evolution-overview" src="https://github.com/varun-shanker/structural-evolution/blob/main/structural-evolution-overview.png">

### Running the model

To evaluate the model on a new sequence, clone this repository and run
```bash
python bin/recommend.py [pdb file] --chain X
```
where `[pdb file]` is the file path to the pdb structure file of the protein or protein complex and `[X]` is the wildtype target chain you wish to evolve. The script will output the top N predicted substitutions, where N can be modified using the arguments (more below).

To recommend mutations to antibody variable domain sequences, we have simply run the above script separately on the heavy and light chains.

We have also made a [Google Colab](https://colab.research.google.com/drive/18QLOmi5yNb1i9wztAzv981Wgk2E4IP4q?usp=sharing) notebook available. However, this notebook requires a full download and installation of the language models for each run and requires Colab Pro instances with a higher memory requirement than the free version of Colab. When making many predictions, we recommend the local installation above, as this will allow you to cache and reuse the models.

### Paper analysis scripts

To reproduce the analysis in the paper, first download and extract data with the commands:
```bash
wget https://zenodo.org/record/6968342/files/data.tar.gz
tar xvf data.tar.gz
```

To acquire mutations to a given antibody, run the command
```bash
bash bin/eval_models.sh [antibody_name]
```
where `[antibody_name]` is one of `medi8852`, `medi_uca`, `mab114`, `mab114_uca`, `s309`, `regn10987`, or `c143`.

DMS experiments can be run with the command
```bash
bash bin/dms.sh
```
