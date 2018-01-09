#!/bin/py

#  alignment-pipeline.py
#
#
#  Created by Kieran Chacko on 10/02/17.
#
#

# Import Modules
import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO

# Load arguments
firstarg=sys.argv[1]    # Reference File
seconarg=sys.argv[2]    # Query File
thirdarg=sys.argv[3]    # Reference Name
fourtarg=sys.argv[4]    # Query Name

# Load data
print ('\n' + 'Loading data...' + '\n')
df_1 = pd.read_csv(firstarg, sep='\t', header = None)
df_2 = pd.read_csv(seconarg, sep='\t', header = None)

columns = [
    "Ref-PROKKA_ID",
    "Mutation_Type",
    "Query-Contig",
    "Ref-Start",
    "Ref-End",
    "Ref-NT",
    "Query-NT",
    "Variant_Type",
    "GAP",
    "Description",
    "Mutation"
    ]

df_1.columns = columns
df_2.columns = columns
df_1['Isolate'] = str(thirdarg)
df_2['Isolate'] = str(fourtarg)


# Combine Data
print ('\n' + 'Finding unique variants...' + '\n')
df = pd.concat([df_1, df_2], ignore_index=True)
df = df.drop_duplicates(subset = ['Ref-PROKKA_ID', 'Ref-Start', 'Ref-End', 'Query-NT', 'Variant_Type', 'Mutation'], keep = False )

# Output Files
print ('\n' + 'Writing output files...' + '\n')
df.to_csv('./' + str(thirdarg) + '_vs_' + str(fourtarg) + '_Unique-Variants.txt', sep='\t', index=False, header=True)

# Finish
print ('\n' + 'Finding unique variants complete...' + '\n')
print ('\n' + 'Author: Kieran Chacko...' + '\n')
