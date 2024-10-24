# gwas-bionets

gwas-bionets is a repository to run different biological network methods, namely [Heinz](https://academic.oup.com/bioinformatics/article/24/13/i223/231653?login=true), [HotNet2](https://www.nature.com/articles/ng.3168), and [SigMod](https://academic.oup.com/bioinformatics/article/33/10/1536/2874362). The input is a typical fileset for GWAS analysis in PLINK 1.9 format (.bim, .bed and .fam). The general workflow for constructing a consensus network is as follows:

[<img src="img/consensus_pipeline.svg" width="700"/>](img/consensus_pipeline.svg)

And for generating a stable consensus network (when the _k_ parameter is greater than 1) the correspoding pipeline is:

[<img src="img/stable_consensus_pipeline.svg" width="700"/>](img/stable_consensus_pipeline.svg)

For the latter image, we used _k_=5 as illustration but you can change to greater number which makes sense in your experiments. Also, 1H,2H, ... corresponds to the solutions outputted by Heinz; 1N, 2N, ... solutions of HotNet2; and 1S, 2S, ... solutions of SigMod. Filtering step is not included in the pipeline, then please assure your data is a priori filtered. BioGRID is used as a source for the PPI but you can choose another network of reference, be advised that a two-column header indicating the connection between two molecules is expected, i.e., 'Official Symbol Interactor A' and 'Official Symbol Interactor B'.

1. This script works with the raw data for splitting it if parametrized with the _k_ parameter. The script is conceived to be modified to provide the parameters.

`bionets_construction_from_data.sh`

2. This script works with the scores previously computed using a software like [MAGMA](https://cncr.nl/research/magma/) for the gene P-values. Again, a _k_ parameter greater than 1 generates k-fold solutions. The script is conceived to be modified to provide the parameters.

`bionets_construction_from_scores.sh`

