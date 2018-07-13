# Build and visualize trees from antibody JSON data #

## Installation

### Requirements

- Anaconda 
- R 
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

Obtain a copy of the compressed input data, `patients_clones.tar.gz`, and place in the input data directory via

```
tar xvzf /path/to/patients_clones.tar.gz -C data/input/
```

After installing requirements, run the pipeline from the bcell-phylo directory:

```
snakemake $TARGET --cluster "qsub -V" -j $NUMBER_OF_CONCURRENT_JOBS
```

Make sure that the `python` executable for the `conda` environment is on your `$PATH`.

## Visualization

After running the pipeline:

```
yarn start
```
