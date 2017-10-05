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
firstarg=sys.argv[1]    # Reference Fasta
seconarg=sys.argv[2]    # Query Fasta
thirdarg=sys.argv[3]    # Reference Gff
fourtarg=sys.argv[4]    # Prefix

# Make Directories
print ('\n' + 'Making Folders...' + '\n')
subprocess.Popen('mkdir ./' + str(fourtarg), shell=True).wait()
subprocess.Popen('mkdir ./' + str(fourtarg) + '/NucDIFF', shell=True).wait()
subprocess.Popen('mkdir ./' + str(fourtarg) + '/Annovar', shell=True).wait()
subprocess.Popen('mkdir ./' + str(fourtarg) + '/tmp-files', shell=True).wait()
subprocess.Popen('cd ./' + str(fourtarg), shell=True).wait()

# Run NucDIFF
print ('\n' + 'Running NucDiff...' + '\n')
subprocess.Popen('python ~/dwns/NucDiff/nucdiff.py ' + str(firstarg) + ' ' + str(seconarg) + ' ' + str(fourtarg) + '/NucDIFF' + ' NucDIFF', shell=True).wait()

# NucDIFF Output to Annovar Input
print ('\n' + 'Converting NucDIFF to Annovar...' + '\n')
with open('./' + str(fourtarg) + '/NucDIFF/results/NucDIFF_ref_snps.gff') as inFile:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile:
        for line in inFile:
            
            line = line.replace(';', '\t')
            line = line.replace('ID=', '')
            line = line.replace('Name=', '')
            line = line.replace('subst_len=', '')
            line = line.replace('ins_len=', '')
            line = line.replace('del_len=', '')
            line = line.replace('query_dir=', '')
            line = line.replace('query_sequence=', '')
            line = line.replace('query_coord=', '')
            line = line.replace('query_bases=', '')
            line = line.replace('ref_bases=', '')
            line = line.replace('color=', '')
            
            if line.startswith('#'):
                pass
            
            else:
                 outFile.write(line)

df = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df.columns = [
    "Ref-Contig",
    "NucDiff-Version",
    "SO-Number",
    "Ref-Start",
    "Ref-End",
    "Gap_1",
    "Gap_2",
    "Gap_3",
    "SNP-ID",
    "Substitution-Type",
    "Substitution-Length",
    "Query-Dir",
    "Query-Contig",
    "Query-Coord=",
    "Query-Base",
    "Ref-Base",
    "Color"
    ]

df = df.loc[:, ['Ref-Contig','Ref-Start','Ref-End', 'Ref-Base', 'Query-Base', 'Substitution-Type', 'Substitution-Length']]            
df.to_csv('./' + str(fourtarg) + '/Annovar/Annovar_Input.txt', sep='\t', index=False, header = None)

# Gff to RefGene
print ('\n' + 'Making RefGene File...' + '\n')
with open(thirdarg) as inFile:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile:
        for line in inFile:
            
            line = line.replace(';', '\t', 1)
            line = line.replace('ID=', '')
            
            if line.startswith('#'):
                pass
            if line.startswith('##FASTA'):
                break
            elif not line.startswith('#'):
                outFile.write(line)

df = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df.columns = [
    "Ref-Contig",
    "Prodigal-Version",
    "Coding-Region",
    "Gene-Start",
    "Gene-End",
    "Gap_1",
    "Strand",
    "Gap_2",
    "PROKKA-ID",
    "Description"
    ]

df['Gene-Start'] = df['Gene-Start'] - 1
df['Gene-Start-2'] = df['Gene-Start']
df['Gene-End-2'] = df['Gene-End']
df['Gene-Start-3'] = df['Gene-Start']
df['Gene-End-3'] = df['Gene-End']
df['Number'] = 1
df['Number-2'] = 1
df['Gap_2'] = 0
df['PROKKA-ID-2'] = df['PROKKA-ID']
df['Start-Stat'] = 'cmpl'
df['End-Stat'] = 'cmpl'
df['Gap_3'] = df['Gap_2']

df.ix[(df['Coding-Region'] == 'tRNA'), ['Gene-Start-2', 'Gene-End-2']] = df['Gene-End']
df.ix[(df['Coding-Region'] == 'tRNA'), ['Start-Stat', 'End-Stat']] = 'none'
df.ix[(df['Coding-Region'] == 'tRNA'), ['Gap_3']] = -1
df.ix[(df['Coding-Region'] == 'tRNA'), ['Coding-Region']] = "mRNA"
df.ix[(df['Coding-Region'] == 'tmRNA'), ['Coding-Region']] = "mRNA"

df = df[['Number',
         'PROKKA-ID',
         'Ref-Contig',
         'Strand',
         'Gene-Start',
         'Gene-End',
         'Gene-Start-2',
         'Gene-End-2',
         'Number-2',
         'Gene-Start-3',
         'Gene-End-3',
         'Gap_2',
         'PROKKA-ID-2',
         'Start-Stat',
         'End-Stat',
         'Gap_3'
         ]]

df = df.sort_values(by = ['PROKKA-ID'], ascending = [0])
df.to_csv('./' + str(fourtarg) + '/Annovar/' + str(fourtarg) + '_refGene.txt', sep='\t', index=False, header = None)

# Split contigs
fasta_sequences = SeqIO.parse(open(firstarg),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    with open('./' + str(fourtarg) + '/Annovar/' + name + '.fa', 'w') as outFile:
        outFile.write('>' + name + '\n' + sequence)
        
# Finalize RefGene File
subprocess.Popen('retrieve_seq_from_fasta.pl -format refGene -outfile ./' + str(fourtarg) + '/Annovar/' + str(fourtarg) + '_refGeneMrna.fa -seqdir ./' + str(fourtarg) + '/Annovar/ ./' + str(fourtarg) + '/Annovar/' + str(fourtarg) + '_refGene.txt', shell=True).wait()

# Gff to refLink
print ('\n' + 'Making Ref-Link File...' + '\n')
with open(thirdarg) as inFile:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile:
        for line in inFile:
            
            line = line.replace(';', '\t', 1)
            line = line.replace('ID=', '')
            
            if line.startswith('#'):
                pass
            if line.startswith('##FASTA'):
                break
            elif not line.startswith('#'):
                outFile.write(line)

df = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df.columns = [
    "Ref-Contig",
    "Prodigal-Version",
    "Coding-Region",
    "Gene-Start",
    "Gene-End",
    "Gap_1",
    "Strand",
    "Gap_2",
    "PROKKA-ID",
    "Description"
    ]

df['Gene'] = df['PROKKA-ID'] + str('.gene')

df = df.loc[:, ['PROKKA-ID',
                'Gene'
               ]]            

df.to_csv('./' + str(fourtarg) + '/Annovar/' + str(fourtarg) +  '_refLink.txt', sep='\t', index=False, header= None)


# Run Annovar
print ('\n' + 'Running Annovar...' + '\n')
subprocess.Popen('annotate_variation.pl -buildver ' + str(fourtarg) + ' ./' + str(fourtarg) + '/Annovar/Annovar_Input.txt ./' + str(fourtarg) + '/Annovar/', shell=True).wait()

# Creating Output Files
print ('\n' + 'Creating Output Files...' + '\n')

# Load Exonic Variant Function files
with open('./' + str(fourtarg) + '/Annovar/Annovar_Input.txt.exonic_variant_function') as inFile:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile:
        for line in inFile:
            
            line = line.replace(':', '\t', 1)
            line = line.replace(':c.', '\t', 1)
            line = line.replace(',', '', 1)
            outFile.write(line)
            
df = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df.columns = [
    "Line",
    "Variant-Effect",
    "PROKKA-ID",
    "Exon",
    "Variant-in-Gene",
    "Ref-Contig",
    "Ref-Start",
    "Ref-End",
    "Start-NT",
    "End-NT",
    "Variant-Type",
    "Variant-Length"
    ]

df['Codon'] = df['Variant-in-Gene'].str.split(':p.').str[-1]

df = df[['PROKKA-ID', 
    "Variant-Effect",
    "Ref-Contig",
    "Ref-Start",
    "Ref-End",
    "Start-NT",
    "End-NT",
    "Variant-Type",
    "Variant-Length",
    "Codon",
        ]]

# Load GFF files
with open(thirdarg) as inFile2:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile2:
        for line in inFile2:
            
            line = line.replace(';', '\t', 1)
            line = line.replace('ID=', '')
            line = line.replace('product=', '\t')
            
            if line.startswith('#'):
                pass
            if line.startswith('##FASTA'):
                break
            elif not line.startswith('#'):
                outFile2.write(line)

df2 = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df2.columns = [
    "Ref-Contig",
    "Prodigal-Version",
    "Coding-Region",
    "Gene-Start",
    "Gene-End",
    "Gap_1",
    "Strand",
    "Gap_2",
    "PROKKA-ID",
    "Prokka-Description",
    "Description"
    ]

df2 = df2[['PROKKA-ID', 'Description']]

# Load Variant Function file
with open('./' + str(fourtarg) + '/Annovar/Annovar_Input.txt.variant_function') as inFile3:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile3:
        for line in inFile3:
            if line.startswith('exonic'):
                pass
            else:
                outFile3.write(line)

df3 = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df3.columns = ['Variant-Effect', 
    "PROKKA-ID",
    "Ref-Contig",
    "Ref-Start",
    "Ref-End",
    "Start-NT",
    "End-NT",
    "Variant-Type",
    "Variant-Length"
        ]

df3['Description'] = 'Non-Coding Region'
df3['Codon'] = 'n/a'

df3 = df3[['PROKKA-ID',
           'Variant-Effect',
           'Ref-Contig',
           'Ref-Start',
           'Ref-End',
           'Start-NT',
           'End-NT',
           'Variant-Type',
           'Variant-Length',
           'Description',
           'Codon'
    ]]

df = pd.merge(df, df2, on=['PROKKA-ID'], how='inner')
df = pd.concat([df, df3])

df = df[['PROKKA-ID',
           'Variant-Effect',
           'Ref-Contig',
           'Ref-Start',
           'Ref-End',
           'Start-NT',
           'End-NT',
           'Variant-Type',
           'Variant-Length',
           'Description',
           'Codon'
    ]]

df = df.sort_values(by = ['Ref-Contig', 'Ref-Start'], ascending = [1, 1])
df.to_csv('./' + str(fourtarg) + '/' + str(fourtarg) +  '_SNVs.txt', sep='\t', index=False, header= None)


# Preparing Annovar Input for SV
print ('\n' + 'Annotating SV Variants...' + '\n')
print ('\n' + 'Converting NucDIFF to Annovar...' + '\n')
with open('./' + str(fourtarg) + '/NucDIFF/results/NucDIFF_ref_struct.gff') as inFile:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile:
        for line in inFile:

            line = line.replace(';', '\t', 3)
            line = line.replace('ID=', '')
            line = line.replace('Name=', '')
            line = line.replace('subst_len=', '')
            line = line.replace('ins_len=', '')
            line = line.replace('del_len=', '')
            line = line.replace('overlap_len=', '')            
            line = line.replace('query_dir=', '')
            line = line.replace('query_sequence=', '')
            line = line.replace('query_coord=', '')
            line = line.replace('query_bases=', '')
            line = line.replace('ref_bases=', '')
            line = line.replace('color=', '')
                  
            if line.startswith('#'):
                pass
            
            else:
                outFile.write(line)

df = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df.columns = [
    "Ref-Contig",
    "NucDiff-Version",
    "SO-Number",
    "Ref-Start",
    "Ref-End",
    "Gap_1",
    "Gap_2",
    "Gap_3",
    "SV-ID",
    "Substitution-Type",
    "Substitution-Length",
    "Other"
    ]

df['Ref-Base'] = 0
df['Query-Base']= 0

df.ix[(df['Substitution-Type'] == 'insertion'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'duplication'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'relocation-insertion'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'relocation-insertion'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'relocation-overlap'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'translocation-overlap'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'translocation-insertion'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'collapsed_repeat'), ['Ref-Base', 'Query-Base']] = ['0', '-']
df.ix[(df['Substitution-Type'] == 'deletion'), ['Ref-Base', 'Query-Base']] = ['0', '-']
df.ix[(df['Substitution-Type'] == 'substitution'), ['Ref-Base', 'Query-Base']] = ['0', '0']
df.ix[(df['Substitution-Type'] == 'unaligned_beginning'), ['Ref-Base', 'Query-Base']] = ['-', '0']
df.ix[(df['Substitution-Type'] == 'unaligned_end'), ['Ref-Base', 'Query-Base']] = ['-', '0']

df = df.loc[:, ['Ref-Contig','Ref-Start','Ref-End', 'Ref-Base', 'Query-Base', 'Substitution-Type', 'Substitution-Length']]            
df.to_csv('./' + str(fourtarg) + '/Annovar/Annovar_SV_Input.txt', sep='\t', index=False, header = True)

# Run Annovar
print ('\n' + 'Running Annovar...' + '\n')
subprocess.Popen('annotate_variation.pl -buildver ' + str(fourtarg) + ' ./' + str(fourtarg) + '/Annovar/Annovar_SV_Input.txt ./' + str(fourtarg) + '/Annovar/', shell=True).wait()

# Creating Output Files
print ('\n' + 'Creating Output Files...' + '\n')
with open('./' + str(fourtarg) + '/Annovar/Annovar_SV_Input.txt.exonic_variant_function') as inFile:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile:
        for line in inFile:
            
            line = line.replace(':', '\t', 1)
            outFile.write(line)

df = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df.columns = [
    "Line",
    "Variant-Effect",
    "PROKKA-ID",
    "Variant-in-Gene",
    "Ref-Contig",
    "Ref-Start",
    "Ref-End",
    "Start-NT",
    "End-NT",
    "Variant-Type",
    "Variant-Length"
    ]

df['Codon'] = df['Variant-in-Gene'].str.split(':p.').str[-1]

df = df[['PROKKA-ID', 
    "Variant-Effect",
    "Ref-Contig",
    "Ref-Start",
    "Ref-End",
    "Start-NT",
    "End-NT",
    "Variant-Type",
    "Variant-Length",
    "Codon"
        ]]

# Load Variant Function file
with open('./' + str(fourtarg) + '/Annovar/Annovar_SV_Input.txt.variant_function') as inFile3:
    with open('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', 'w') as outFile3:
        for line in inFile3:
            if line.startswith('exonic'):
                pass
            else:
                outFile3.write(line)

df2 = pd.read_csv('./' + str(fourtarg) + '/tmp-files/tmp-file.txt', sep='\t', header = None)

df2.columns = ['Variant-Effect', 
    "PROKKA-ID",
    "Ref-Contig",
    "Ref-Start",
    "Ref-End",
    "Start-NT",
    "End-NT",
    "Variant-Type",
    "Variant-Length"
        ]

df2['Description'] = 'Non-Coding Region'
df2['Codon'] = 'n/a'

df2 = df2[['PROKKA-ID',
           'Variant-Effect',
           'Ref-Contig',
           'Ref-Start',
           'Ref-End',
           'Start-NT',
           'End-NT',
           'Variant-Type',
           'Variant-Length',
           'Description',
           'Codon'
    ]]

df = pd.concat([df, df2])

df = df[['PROKKA-ID',
           'Variant-Effect',
           'Ref-Contig',
           'Ref-Start',
           'Ref-End',
           'Start-NT',
           'End-NT',
           'Variant-Type',
           'Variant-Length',
           'Description',
           'Codon'
    ]]

df = df.sort_values(by = ['Ref-Contig', 'Ref-Start'], ascending = [1, 1])
df.to_csv('./' + str(fourtarg) + '/' + str(fourtarg) +  '_SVs.txt', sep='\t', index=False, header= True)

# Cleaning
print ('\n' + 'Removing temp files...' + '\n')
subprocess.Popen('rm -r ./' + str(fourtarg) + '/tmp-files', shell=True).wait()

# Cleaning
print ('\n' + 'Pairwise Pipeline Complete...' + '\n')
print ('\n' + 'Author: Kieran Chacko...' + '\n')
