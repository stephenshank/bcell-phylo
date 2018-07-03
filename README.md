## Build and visualize trees from antibody JSON data. ##



First:
	Need to activate the virtual environment to run the pipeline.
	Copy and paste these commands:
		- if you have Mac:
		`conda create --name bcell --file spec-file.txt`
		`source activate bcell`
		
		- if you have Windows
		- if you have Linux

Second: 
	cd in the bcell-phylo directory and run the command 
	`bash bcell-pipe.sh`
	
## Requirements

- Anaconda 
- R 
- Python 
- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [yarn](https://yarnpkg.com/en/)

Obtain a copy of the compressed input data, `patients_clones.tar.gz`, and place in the input data directory via

```
tar xvzf /path/to/patients_clones.tar.gz -C data/input/
```

Install JavaScript dependencies with

```
yarn
```

## Pipeline

After installing requirements, from the base of this directory run

To run the pipeline, cd in the bcell-phylo directory:
```
bash bcell-pipe.sh
```

## Visualization

After running the pipeline:

```
yarn start
```

```
go to data/out/pretty_pictures  and enjoy
```

```
when you are finished, run bash cleaner.sh
```
