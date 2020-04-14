# panFEAR
pangenome Functional Enrichment Analysis from Roary

Help message:
```
python panFEAR.py -h
usage: panFEAR.py [-h] [-gpa GENE_PRESENCE_ABSENCE] [-fhp] [-b N_BINS]
                  [-up UPPER_PERCENT] [-lp LOWER_PERCENT] [-c CORRECTION]
                  [-o OUTPUT_FOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -gpa GENE_PRESENCE_ABSENCE, --gene_presence_absence GENE_PRESENCE_ABSENCE
                        Roary output file "gene_presence_absence.csv" from
                        which information is retrieved
  -fhp, --filter_hypothetical_proteins
                        Flag for excluding "hypothetical protein" from the
                        analysed functions
  -b N_BINS, --n_bins N_BINS
                        Number of bins in the pangenome
  -up UPPER_PERCENT, --upper_percent UPPER_PERCENT
                        Upper percent among the pangenome of the group of
                        which panFEAR analyses functional enrichment
  -lp LOWER_PERCENT, --lower_percent LOWER_PERCENT
                        Lower percent among the pangenome of the group of
                        which panFEAR analyses functional enrichment
  -c CORRECTION, --correction CORRECTION
                        P-value correction method from multiple tests bias
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        ABSOLUTE! path where to store the produced outputs
```
