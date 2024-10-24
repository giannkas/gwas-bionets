# gwas-bionets

gwas-bionets is a repository to run different biological network methods, namely [Heinz](https://academic.oup.com/bioinformatics/article/24/13/i223/231653?login=true), [HotNet2](https://www.nature.com/articles/ng.3168), and [SigMod](https://academic.oup.com/bioinformatics/article/33/10/1536/2874362). The input is a typical fileset for GWAS analysis in PLINK 1.9 format (.bim, .bed and .fam). The general workflow for constructing a consensus network is as follows:

[<img src="img/consensus_pipeline.svg" width="700"/>](img/consensus_pipeline.svg)

And for generating a stable consensus network (when the _k_ parameter is greater than 1) the correspoding pipeline is:

[<img src="img/stable_consensus_pipeline.svg" width="700"/>](img/stable_consensus_pipeline.svg)

For the latter image, we used _k_=5 as illustration but you can change to a greater number which makes sense in your experiments. Also, 1H,2H, ... corresponds to the solutions outputted by Heinz; 1N, 2N, ... solutions of HotNet2; and 1S, 2S, ... solutions of SigMod. We omitted the filtering step in the pipeline, so please ensure your data is filtered prior. We used BioGRID as a source for the PPI, but you can choose another network of reference, considering a two-column header indicating the connection between two molecules (genes, proteins, ...) as the format, i.e., 'Official Symbol Interactor A' and 'Official Symbol Interactor B'.

## Software requirements

This code uses software already outdated so it is encouraged to follow the installation process; otherwise, it may not work. An update of the methods and requirements is planned but still needs to be carried out. Due to not having administrative rights or to avoiding conflicts, most software needs to be installed locally, and paths to these installations must be redirected.

Ideally, create a folder in your home directory to store all software. For example:

```bash
  mkdir ~/bin
```

**Install Java (required for installing nextflow)**

Create a "java" folder in the home directory and navigate to it.

```bash
mkdir ~/bin/java
cd ~/bin/java
```


1. This script works with the raw data for splitting it if parametrized with the _k_ parameter. The script has the needed parameters to be filled by the user, clearly, you can run each of the steps within the script separately.

`bionets_construction_from_data.sh`

2. This script works with the scores previously computed using a software like [MAGMA](https://cncr.nl/research/magma/) for the gene P-values. Again, a _k_ parameter greater than 1 generates k-fold solutions. As above, we conceived the script to be modified to provide the parameters.

`bionets_construction_from_scores.sh`

