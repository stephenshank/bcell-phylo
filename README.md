# Build and visualize trees from antibody JSON data #

## Installation

### Requirements

- SnakeMake
- Anaconda  
- Python 
- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [yarn](https://yarnpkg.com/en/)

### Install instructions

```
conda env create -f environment.yml
source activate bcell
bash install.sh
yarn
```

## Pipeline

Obtain a copy of the compressed input data, `bcell-phylo_Ver2.tar.gz`, and place in the input data directory via

```
tar xvzf /path/to/bcell-phylo_Ver2.tar.gz -C data/input/
```

After installing requirements, run the pipeline from the bcell-phylo directory:

```
snakemake -j $NUMBER_OF_CONCURRENT_JOBS $TARGET
```

Or, to distribute jobs with Sun Grid Engine:


```
snakemake $TARGET --cluster "qsub -V" -j $NUMBER_OF_CONCURRENT_JOBS
```

Make sure that the `python` executable for the `conda` environment is on your `$PATH`.

## Visualization

After running the pipeline:

```
yarn start
```
