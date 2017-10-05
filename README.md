## pairwise-pipeline

#### Pairwise-pipeline uses mummer to align two closely related genomes to identify both single and structural variants. The outputs are passed to annovar for annotation of the variants.

#### Required software:
#### Python
#### Python Libraries: os, sys, subprocess, pandas, SeqIO
#### Mummer
#### Annovar

#### Run Instructions:
#### python ./pairwise-pipeline.py {Reference_Fasta} {Query_Fasta} {Reference_gff} {Prefix}
