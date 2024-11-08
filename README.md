# gwas-bionets

gwas-bionets is a repository to run different biological network methods, namely [Heinz](https://academic.oup.com/bioinformatics/article/24/13/i223/231653?login=true), [HotNet2](https://www.nature.com/articles/ng.3168), and [SigMod](https://academic.oup.com/bioinformatics/article/33/10/1536/2874362). The input is a typical fileset for GWAS analysis in PLINK 1.9 format (.bim, .bed and .fam). The general workflow for constructing a consensus network is as follows:

[<img src="img/consensus_pipeline.svg" width="700"/>](img/consensus_pipeline.svg)

And for generating a stable consensus network (when the _k_ parameter is greater than 1) the correspoding pipeline is:

[<img src="img/stable_consensus_pipeline.svg" width="700"/>](img/stable_consensus_pipeline.svg)

For the latter image, we used _k_=5 as illustration but you can change to a greater number which makes sense in your experiments. Also, 1H,2H, ... corresponds to the solutions outputted by Heinz; 1N, 2N, ... solutions of HotNet2; and 1S, 2S, ... solutions of SigMod. We omitted the filtering step in the pipeline, so please ensure your data is filtered prior. We used BioGRID as a source for the PPI, but you can choose another network of reference, considering a two-column header indicating the connection between two molecules (genes, proteins, ...) as the format, i.e., 'Official Symbol Interactor A' and 'Official Symbol Interactor B'.

## Software requirements

This code uses software already outdated so it is encouraged to follow the installation process; otherwise, it may not work. An update of the methods and requirements is planned but still needs to be carried out. Due to not having administrative rights or to avoiding conflicts, most software needs to be installed locally or within an environment, and paths to these installations must be redirected. It is assumed that you run these commands on a Unix-like machine.

Ideally, create a folder in your home directory to store all software. For example:

```bash
mkdir ~/bin
```

**Install Java (required for installing nextflow)**

Create a "java" folder in the software directory and navigate to it.

```bash
mkdir ~/bin/java
cd ~/bin/java
```

Download the `x64` version of Java as a `tar.gz` file from the Oracle website into your machine and decompress it.

```bash
wget https://download.oracle.com/java/17/archive/jdk-17.0.10_linux-x64_bin.tar.gz
tar xzfv jdk-17.0.10_linux-x64_bin.tar.gz 
rm jdk-17.0.10_linux-x64_bin.tar.gz
```

Export the path to the `bin` directory of this folder into the system variable `$PATH` to make Java executable. Also, export the `$JAVA_HOME` variable indicating the root directory. Ideally, add these to `~/.bashrc` to avoid repeating the process on each server connection or reboot, eg.

```bash
export PATH=/home/username/bin/java/jdk-17.0.10/bin:$PATH
export JAVA_HOME=/home/username/bin/java/jdk-17.0.10
```

You may need to source `.bashrc` file before checking installation, so type:

```bash
source ~/.bashrc
```

Test the installation:

```bash
java -version
```

You should see something like:

```bash
java version "17.0.10" 2024-01-16 LTS
Java(TM) SE Runtime Environment (build 17.0.10+11-LTS-240)
Java HotSpot(TM) 64-Bit Server VM (build 17.0.10+11-LTS-240, mixed mode, sharing)
```

**Install Nextflow**

Change directory to our software directory.

```bash
cd ~/bin
```

Download Nextflow version 22.10.4 and decompress:

```bash
wget https://github.com/nextflow-io/nextflow/archive/refs/tags/v22.10.4.tar.gz
tar -xzvf /gwas-bionets/nextflow/v-22.10.4.tar.gz
```

Change folder name and navigate to it.

```bash
mv nextflow-22.10.4 nextflow
cd nextflow
```

Compile and install it:

```bash
make compile
make pack
make install
```

Add to your path:

```bash
export PATH=$PATH:/home/username/bin/nextflow
```

Test the installation:

```bash
nextflow -version
```

You should see something like:

```bash
      N E X T F L O W
      version 22.10.4 build 5837
      created 26-10-2024 08:56 UTC (10:56 CEST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

**Install MAGMA**

You can follow the installation instructions for MAGMA at its website (version 1.10): [Multi-marker Analysis of GenoMic Annotation](https://cncr.nl/research/magma/). According to the documentation MAGMA: is a self-contained executable and does not need to be installed. 

Go to our software directory and create a 'magma' folder where the binaries will be located.

```bash
cd ~/bin
mkdir magma
```

Download the `zip` file.

```bash
curl -v "https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download" -H "Accept-Encoding: zip" > magma_v1.10.zip
```

Descompress `magma_v1.10.zip` file. We will decompress to the created folder `magma` (unsing `unzip -d` option)

```bash
unzip magma_v1.10.zip -d magma
```

Remove the `zip` file.

```bash
rm magma_v1.10.zip
```

Test MAGMA

```bash
./magma/magma
```

You should see something like:

```bash
No arguments specified. Please consult manual for usage instructions.

Exiting MAGMA. Goodbye.
```

You can decide to include MAGMA's location into the PATH variable so it is called system-wide under your session. Otherwise, you must indicate where to find MAGMA when calling `magma_calc.nf` script; it has a parameter named `magma` for this purpose.

**Install PLINK**

Similarly, you can install PLINK (version 1.9) from its website: [population linkage](https://www.cog-genomics.org/plink/1.9/). PLINK is also self-contained executable so either you add to your path or reference the executable when using it.

Again, please locate our software directory and create a 'plink' folder where the binaries will be saved.

```bash
cd ~/bin
mkdir plink
```

Download the `zip` file

```bash
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip
```

Decompress `plink_linux_x86_64_20241022.zip` file. We will decompress to the created folder `plink` (unsing `unzip -d` option)

```bash
unzip plink_linux_x86_64_20241022.zip -d plink
```

Remove the `zip` file.

```bash
rm plink_linux_x86_64_20241022.zip
```

Test PLINK

```bash
./plink/plink --version
```

You should see something like:

```bash
PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)
```

We have to include PLINK's location into our PATH variable because there is no parameter in the pipeline to reference its location. Conversely, we did include a PLINK parameter to indicate PLINK version, either 1 or 2. Then, to add plink to the environment variable, you proceed as follow:

```bash
export PATH=$PATH:/home/username/bin/plink
```

You can add it to your `.bashrc` file to make it permanent. And then, `source ~/.bashrc` to apply changes in your current session.


**Install R and some packages (required for the methods)**

If you dont't have R installed in your machine (add your superuse credentials if needed to install software), then proceed as follows (or you can check instructions from [The Comprehensive R archive Network](https://cloud.r-project.org/)):

```bash
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt-get -y install --no-install-recommends r-base r-base-dev
```

Setup the general CRAN repo.

```bash
echo 'local({
    r <- getOption("repos")
    r["CRAN"] = "https://cloud.r-project.org/"
    options(repos = r)
  })' >> /etc/R/Rprofile.site
```

Install Bioconductor (`BiocManager`), `twilight` and `BioNet` (the latter contains the necessary files to use Heinz method).

```bash
R -e "install.packages('BiocManager')"
R -e "BiocManager::install('BioNet')"
R -e "BiocManager::install('twilight')"
```

Install R packages, `tidyverse`, `cowplot`, `igraph` and `gprofiler2`:

```bash
R -e "install.packages(c('tidyverse', 'cowplot', 'igraph', 'gprofiler2'))" 
```

**Install Python2**

HotNet2 uses python2 to run some of its scripts; nowadays, it may be troublesome to install `python 2.7` so you may opt to use a conda environment (see below for instructions).

```bash
apt-get -y install python2-dev python2 python-pip
```

Add some python2 libraries needed for HotNet2:

```bash
pip2 install numpy==1.12.1 scipy==0.19.0 networkx==1.11 h5py==2.7.0
```

**Install Python2 using a Conda environment**

We assume you have Conda installed, if not please refer to [Conda installing instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), you should follow the steps under the section 'Regular installation'. After you have installed Conda, you can continue with the next instructions.

We will create a Conda environment called `gwas-bionets` with `Python 2.7` version installed and its corresponding `pip` version.

```bash
conda create -n gwas-bionets python=2.7 pip
```

You should see a last message saying:

```bash
# To activate this environment, use
#
#     $ conda activate gwas-bionets
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

So we will activate our environment.

```bash
conda activate gwas-bionets
```

Within the environment, we install the required packages for HotNet2:

```bash
conda install -c conda-forge numpy=1.12.1 scipy=0.19.0 networkx=1.11 h5py=2.7.0
```


**Install SigMod**

You can install SigMod (version 2) from this website: [Strongly Interconnected Gene MODule](https://github.com/YuanlongLiu/SigMod/tree/20c561876d87a0faca632a6b93882fcffd719b17). This is an R package and it suffices to assign the parameter `sigmod_path` when calling the `bionets.nf` script, eg. `--sigmod_path="~/bin/SigMod_v2"`.

You change directory to our software directory `bin` folder:

```bash
cd ~/bin
```

Download the `zip` file and decompress it (no need to create a folder since there is a folder inside including the code and manual):

```bash
wget https://github.com/YuanlongLiu/SigMod/raw/20c561876d87a0faca632a6b93882fcffd719b17/SigMod_v2.zip
unzip SigMod_v2.zip
```

Change folder name

```bash
mv SigMod_v2 sigmod
```


**Install HotNet2**

You can install HotNet2 from this website: [HotNet2](https://github.com/raphael-group/hotnet2). Save the code in a folder and name it `hotnet2` whose location can therefore reference in the parameter `hotnet2_path` when calling `bionets.nf` script, eg. `--hotnet2_path="~/bin/hotnet2"`

Move to out parent folder.

```bash
cd ~/bin
```

And to donwload the code, you can clone it using `git`. This will create folder called `hotnet2`

```bash
git clone https://github.com/raphael-group/hotnet2.git
```

_Or_, in case you cannot use `git`, you can use `wget` to download the `zip` file and decompress it.

```bash
wget --no-check-certificate -O hotnet2.zip https://github.com/raphael-group/hotnet2/archive/refs/heads/master.zip
unzip hotnet2.zip
```

Change folder name.

```bash
mv hotnet2-master hotnet2
```

And remove the `zip` file (we remove them so we have a clean repository).

```bash
rm hotnet2.zip
```

**Install Heinz**

You have already installed it when installing the `BioNet` package from Bioconductor :-)

After that, we are all set!

## Main Scripts

1. This script works with the raw data for splitting it if parametrized with the _k_ parameter. The script has the needed parameters to be filled by the user, clearly, you can run each of the steps within the script separately.

`bionets_construction_from_data.sh`

2. This script works with the scores previously computed using a software like [MAGMA](https://cncr.nl/research/magma/) for the gene P-values. Again, a _k_ parameter greater than 1 generates k-fold solutions. As above, we conceived the script to be modified to provide the parameters.

`bionets_construction_from_scores.sh`

