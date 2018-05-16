# bcell-phylo

Build and visualize phylogenetic trees from antibody JSON data.

## Requirements

- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [yarn](https://yarnpkg.com/en/)

Obtain a copy of the compressed input data, `bcells.tar.gz`, and place in the input data directory via

```
tar xvzf /path/to/bcells.tar.gz -C data/input/
```

Install JavaScript dependencies with

```
yarn
```

## Pipeline

After installing requirements, from the base of this directory run

```
bash pipeline/main.sh
```

## Visualization

After running the pipeline:

```
yarn start
```

