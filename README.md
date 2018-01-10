## pairwise-pipeline.py

#### Pairwise-pipeline uses mummer to align two closely related genomes to identify both single and structural variants. The outputs are passed to annovar for annotation of the variants.
#
### Required software:
#### Python
#### Python Libraries: os, sys, subprocess, pandas, SeqIO
#### Mummer
#### Annovar

#
### Run Instructions:
#### python ./pairwise-pipeline.py {Reference_Fasta} {Query_Fasta} {Reference_gff} {Prefix}

## reference-readmap.py

#### Reference-readmap.py aligns illumina reads to a provided complete reference genome. Alignment is passed through freebayes, bcftools and annovar for variant calling, filtering and annotation.
#
### Required software:
#### Python
#### Python Libraries: os, sys, subprocess, pandas, SeqIO
#### bwa
#### samtools
#### freebayes
#### bcftools
#### Annovar

#
## unique-variants.py

#### Unique-Variants.py takes the output files from pairwise-pipeline.py or reference-readmap.py and compares two against each other. The result is a list of SNVs that are unique to each sample.

### Required software:
#### Python
#### Python libraries

### Run Instructions:
#### python ./unqiue-variants.py {Reference Input File} {Query Input File} {Reference Prefix} {Query Prefix}

