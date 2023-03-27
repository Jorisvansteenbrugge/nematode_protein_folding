# Nematode Protein Folding

### Protein Structural alignment script
- batch_align_pdb.py

This script parses two file lists of pdb files and performs pairwise structural alignments using FATCAT2. The output is stored in csv format.
The script can run alignments in parallel (-t option) to speed up the process.

Example command:

```{bash}
python batch_align_pdb.py --listA structure_A.pdb structure_B.pdb --listB structure_C.pdb structure_D.pdb -t 30 --out output.csv
```
