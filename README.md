# Nematode Protein Folding

- ### batch_align_pdb.py

This script parses two file lists of pdb files and performs pairwise structural alignments using FATCAT2. The output is stored in csv format.
The script can run alignments in parallel (-t option) to speed up the process.

Example command:

```{bash}
python batch_align_pdb.py --listA structure_A.pdb structure_B.pdb --listB structure_C.pdb structure_D.pdb -t 30 --out output.csv
```

- ### main_analysis.R
Contains the steps necesary to reproduce most of the results and figures. Uses the pairwise comparison table from `batch_align_pdb.py` as input files (https://doi.org/10.6084/m9.figshare.22340365)


- ### elegans_lookup.py

Takes the following csv files and outputs a table combinding those sources.
- the G. pallida/C. elegans alignment csv with corrected pvalues - https://doi.org/10.6084/m9.figshare.22339975
- A list containing unannotated G. pallida IDs - [knowable_genome/gpal_unannotated_gene_ids.csv](knowable_genome/gpal_unannotated_gene_ids.csv)
- C. elegans functional descriptions  - [celegans/elegans_functional_annotations.tsv](celegans/elegans_functional_annotations.tsv)
 
