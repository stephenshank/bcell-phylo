# anti\_JSON_FASTA

Build trees from antibody JSON data.

## Requirements

- mafft
- FastTree
- yarn

Obtain JSON data and place in `data` directory.

## Pipeline

-for non-replicate data, put your data into the data_in folder
```
bash pipeline.sh
```
-for replicate data, put your data into the data_rep_in folder
```
bash pipeline_rep.sh
```
## Downstream files
-you can visualize both replicate and non replicate data in a convenient csv format
## Visualization

```
yarn start
```

