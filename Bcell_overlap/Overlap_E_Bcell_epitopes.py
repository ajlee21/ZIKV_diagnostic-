
# coding: utf-8

# In[1]:

#--------------------------------------------------------------------------------
# READ ME
# Overlap of sites within surface exposed diagnostic peptides and E B-cell epitopes
# For all FLAV species
#--------------------------------------------------------------------------------

import csv
from scipy import *
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy
import os


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# In[2]:

#--------------------------------------------------------------------------------
# Import files:
# E protein alignment FASTA
# Immune epitope data
# 15-mer data for all FLAV species
#--------------------------------------------------------------------------------
file_MSA = 'J:\Informatics\BPickett\ZIKV_manuscript\submission_initial\Revision_v2\Bcell_overlap\Flavi_E-reduced_combined1.fasta'
file_epitope = 'J:\Informatics\BPickett\ZIKV_manuscript\submission_initial\Revision_v2\Bcell_overlap\Flav_Bcell_epitope_continuous.csv'
file_epitope_dis = 'J:\Informatics\BPickett\ZIKV_manuscript\submission_initial\Revision_v2\Bcell_overlap\Flav_Bcell_epitope_discontinuous.csv'
file_15mer = 'J:\Informatics\BPickett\ZIKV_manuscript\submission_initial\Revision_v2\Bcell_overlap\Count_surfDiagPep_sites_E.xls'
file_out = 'J:\Informatics\BPickett\ZIKV_manuscript\submission_initial\Revision_v2\Bcell_overlap\Flav_Bcell_epitope_coord.csv'

protein = 'E'

input = open(file_MSA, "rU")
alignment = SeqIO.parse(input, "fasta")

epitope_data = pd.read_csv(file_epitope, sep=',',header=None)
source = epitope_data.values[1:,1]
epitope_list=epitope_data.values[1:,2]

epitope_dis_data = pd.read_csv(file_epitope_dis, sep=',',header=None)
discontinuous_epitope = epitope_dis_data.values[1:,1]

xls = pd.ExcelFile(file_15mer)

if protein == 'E':
    numSitesTot = 509


# In[3]:

#--------------------------------------------------------------------------------
# Create vector of virus type
#--------------------------------------------------------------------------------

msa_id=[]
seq=[]
for record in alignment:
    seq_id = record.id.split('|')[2]
    msa_id.append(seq_id)
    seq.append(record.seq)


#--------------------------------------------------------------------------------
# Find coordinates of immune epitopes within aligned sequences
#--------------------------------------------------------------------------------
start_ls = []
end_ls = []
epitope_ls = []
missing_epitope = 0

for i in range(len(epitope_list)):
    array_id =  array(msa_id)
    ind = where(array_id == source[i])[0]
    # Check for substring "Dengue_virus" in string "Dengue_virus_1"
    if len(ind) == 0:
        ind = where(numpy.char.find(array_id,source[i]) == 0)[0]
    start_all = []
    end_all = []
    for j in range(len(ind)):
        aligned_seq = str(seq[j])
        epitope = epitope_list[i]
        if (aligned_seq.find(epitope) != -1):
            start = aligned_seq.index(epitope)+1
            end = start+len(epitope) -1
            start_all.append(start)
            end_all.append(end)

    uniq_start = unique(start_all)
    uniq_end = unique(end_all)
    if (len(uniq_start)==1) and (len(uniq_end)==1):
        print(epitope, start, end)
        start_ls.append(uniq_start)
        end_ls.append(uniq_end)
        epitope_ls.append(epitope)
    elif (len(uniq_start)>1) or (len(uniq_end)>1):
        print('Warning: multiple positions for this epitope were found ', epitope)
        for k in range(len(uniq_start)):
             print(epitope, uniq_start[k], uniq_end[k])
    else:
        print('Warning: epitope not found')
        missing_epitope+=1
        
print "No. linear epitopes: %d" % len(start_ls)
print "No. missing epitopes: %d" % missing_epitope    
assert(len(start_ls) == len(end_ls))    
   
 


# In[4]:

#--------------------------------------------------------------------------------
# Create list of AA positions within B-cell epitope
#--------------------------------------------------------------------------------
aa = range(1,numSitesTot +1, 1)
epitope_count = zeros(numSitesTot + 1)
num_continuous = len(start_ls)

for position in aa:
    for j in range(num_continuous):
        if position in range(start_ls[j], end_ls[j] + 1, 1):
            epitope_count[position] = 1
    for jj in discontinuous_epitope:
        epitope_count[int(jj)] = 1

epitope_count_sum = sum(epitope_count)

print "No. of B-cell eptitopes: %d " % epitope_count_sum

#--------------------------------------------------------------------------------
# Output .csv file with boolean 1/0 if site is located within a B-cell epitope
#--------------------------------------------------------------------------------
with open(file_out,'wb') as f:
    writer = csv.writer(f)
    for epitope_bool in epitope_count:
        writer.writerow([epitope_bool])
        

#--------------------------------------------------------------------------------
# Create list of AA positions within diagnostic 15-mer peptides
#--------------------------------------------------------------------------------
sheet_names=xls.sheet_names
num_sheets = len(sheet_names)
for iSheet in range(num_sheets):
    data = xls.parse(iSheet)
    nrow = data.shape[0]
    peptide_count = zeros(nrow+1)
    print(len(peptide_count))
    start_pep_ls = []
    end_pep_ls = []
    
    selected_15mer = data.icol(2)
    aa_15mer = data.icol(0)
    for i in range(len(selected_15mer)):
        if int(selected_15mer[i])>0:
            start_pep = int(aa_15mer[i])
            start_pep_ls.append(start_pep)
            end_pep_ls.append(start_pep + 14)
        
    
    
    for position in aa:
        for j in range(len(start_pep_ls)):
            if end_pep_ls[j]>= numSitesTot:
                if position in range(start_pep_ls[j], numSitesTot + 1, 1):
                    peptide_count[position] = 1
            else:
                if position in range(start_pep_ls[j], end_pep_ls[j] + 1, 1):
                    peptide_count[position] = 1

    peptide_count_sum = sum(peptide_count)
    assert(len(start_pep_ls) == len(end_pep_ls))
    
    #--------------------------------------------------------------------------------
    # Calculate the number of sites that overlap between the selected 15-mer peptides 
    # and the B-cell epitopes
    #--------------------------------------------------------------------------------
    overlap_count = 0

    for position in aa:
        if peptide_count[position] == 1 and epitope_count[position] == 1:
            overlap_count +=1
            
    print "No. of sites within surface diagnostic 15-mers in %s: %d " % (sheet_names[iSheet], peptide_count_sum)
    print "No. of sites within surface diagnostic 15-mers AND B-cells in %s: %d " % (sheet_names[iSheet], overlap_count)


# In[ ]:




# In[ ]:



