
# coding: utf-8

# In[4]:

#------------------------------------------------------------------------------------------------------------------------------
# READ ME
# Calculate sensitivity/specificity scores for each position
#------------------------------------------------------------------------------------------------------------------------------
from scipy import *
import csv
import matplotlib.pyplot as plt
import glob
import os

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment


# In[5]:

#------------------------------------------------------------------------------------------------------------------------------
# Import all reformated sequence .txt files from directory
#------------------------------------------------------------------------------------------------------------------------------
dir_ls = {'E': ('J:\Informatics\BPickett\ZIKV_manuscript\submission_initial\Revision_v2\MetaCATS_Results\E\*rResultChisqTest.txt'),
         'NS1': ('J:\Informatics\BPickett\ZIKV_manuscript\submission_initial\Revision_v2\MetaCATS_Results\NS1\*rResultChisqTest.txt')}

dir_pick = 'NS1'
dir_path = dir_ls[dir_pick]


# In[6]:

#------------------------------------------------------------------------------------------------------------------------------
# Identify diagnostic residue: majority residue 
#     if "-" is dominant than sensitivity and specificity scores = NA
#     "-" treated as a different residue in calculations
# Calculate sensitivity and specificity based on the diagnostic residue
#------------------------------------------------------------------------------------------------------------------------------
for file in glob.glob(dir_path):
    fileName = (file.split('\\')[-1])
    virusName = fileName.split('-')[0]
    print(fileName)

    ifile = open(file,"rb")
    reader = csv.reader(ifile, delimiter ='\t')

    header = []
    region = []
    position = []
    chi = []
    pval = []
    df = []
    group1 = []
    group2 = []
    diagnostic_residue = []
    sensitivity = []
    specificity = []
    F1 = []

    rownum = 1
    reader.next()
    for row in reader:
        if ('\n' not in row):
            TP = 0.0
            FN = 0.0
            TN = 0.0
            FP = 0.0

            if rownum == 1:
                header = row
                print header        
            else:
                colnum = 0
                for col in row:
                    if header[colnum] == '1 "Site Number':
                        if col != '':
                            num = col.replace('Site','')
                            position.append(num)
                        print "Position: ", num
                    elif header[colnum] == "Chi-Square Score":
                        if col != '':
                            pval.append(col)
                        print "P-Value: ", col
                    elif header[colnum] == "P-Value":
                        if col != '':
                            chi.append(col)
                        print "Chi: ", col
                    elif header[colnum] == "Degrees of Freedom":
                        if col != '':
                            df.append(col)
                        print "DF: ", col
                    elif header[colnum] == "Residue Diversity Between Groups":
                        if col == "":
                            col = "1000 -"
                        grp_ls = col.split("|")
                        
                        grp1 = grp_ls[0]
                        grp1_start = grp1.find('(')
                        grp1_end = grp1.find(')')
                        grp1_str = grp1[grp1_start+1:grp1_end]
                        group1.append(grp1_str.replace('_',' '))
                        
                        grp2 = grp_ls[1]
                        grp2_start = grp2.find('(')
                        grp2_end = grp2.find(')')
                        grp2_str = grp2[grp2_start+1:grp2_end]
                        group2.append(grp2_str.replace('_',' '))
                        
                        #------------------------------
                        # Calculate scores for Group 1
                        #------------------------------
                        residue_list = grp1_str.split(",")
                        print "Group 1: ", residue_list
                        
                        grp1_residue_count_dict = {}
                        for i in range(len(residue_list)):
                            residue_count = residue_list[i].split("_")
                            for j in range(len(residue_count)):
                                if residue_count[j].isdigit():
                                    count_tmp = int((residue_count[j]))
                                elif residue_count[j].isalpha():
                                    grp1_residue_count_dict[residue_count[j]] = count_tmp
                                elif residue_count[j] == '-':
                                    grp1_residue_count_dict[residue_count[j]] = count_tmp

                        dominant = max(grp1_residue_count_dict.values())
                        dominant_residue = ''
                        for key, value in grp1_residue_count_dict.items():
                            if value == dominant:
                                dominant_residue = dominant_residue+key
                        if len(dominant_residue) == 1:
                                diagnostic_residue.append(dominant_residue)
                        elif len(dominant_residue) > 1:
                                dominant_residue = 'NA'
                                diagnostic_residue.append('NA')
                        print "Diagnostic Residue: ", dominant_residue                

                        if dominant_residue == '-' or dominant_residue == 'NA':
                            sensitivity.append('NA')
                        else:
                            TP = float(grp1_residue_count_dict[dominant_residue])
                            for key, value in grp1_residue_count_dict.items():
                                if key != dominant_residue:
                                    FN += float(grp1_residue_count_dict[key])
                            sensitivity.append(TP/(TP+FN))
                        
                        #------------------------------    
                        # Calculate scores for Group 2
                        #------------------------------
                        residue_list = grp2_str.split(",")
                        print "Group 2: ", residue_list
                        
                        grp2_residue_count_dict = {}
                        for i in range(len(residue_list)):
                            residue_count = residue_list[i].split("_")
                            for j in range(len(residue_count)):
                                if residue_count[j].isdigit():
                                    count_tmp = int((residue_count[j]))
                                elif residue_count[j].isalpha():
                                    grp2_residue_count_dict[residue_count[j]] = count_tmp
                                elif residue_count[j] == '-':
                                    grp2_residue_count_dict[residue_count[j]] = count_tmp

                        if dominant_residue == '-' or dominant_residue == 'NA':
                            specificity.append('NA')   
                        else:
                            if dominant_residue in grp2_residue_count_dict.keys():
                                FP = float(grp2_residue_count_dict[dominant_residue])
                            elif dominant_residue not in grp2_residue_count_dict.keys():
                                FP = 0.0
                            for key, value in grp2_residue_count_dict.items():
                                if key != dominant_residue:
                                    TN += float(grp2_residue_count_dict[key])
                            specificity.append(TN/(FP+TN))
                           
                        #-----------------------------------------------------------
                        # Compute F-measure score
                        # Best = 1, 0 = worst
                        # harmonic mean of precision (TP/number of positives) and 
                        # recall (TP/number that should be positives)
                        #-----------------------------------------------------------
                        if sensitivity[-1] != "NA" or specificity[-1] != "NA":
                            f1 = (2*TP)/(2*TP + FP + FN)
                            F1.append(f1)
                        else:
                            F1.append("NA")
                        
                    colnum += 1
            rownum += 1

    ifile.close()
    
    #print(len(position), len(chi), len(pval), len(df), len(group1), len(group2), len(diagnostic_residue), len(sensitivity), len(specificity))
    assert(len(position) == len(chi) == len(pval) == len(df) == len(group1) == len(group2) == len(diagnostic_residue) == len(sensitivity) == len(specificity) ==len(F1))

    header = []
    header[0:5] = ['Site', 'Chi-square value', 'p-value', 'Degree of freedom','Group1 Residues', 'Group2 Residues'] 
    header.append("Diagnostic Residue")
    header.append("Sensitivity")
    header.append("Specificity")
    header.append("F-measure")

    data = ndarray(shape=(rownum-1, len(header)), dtype=object)
    data[0,:] = header
    data[1:,0] = position
    data[1:,1] = chi
    data[1:,2] = pval
    data[1:,3] = df
    data[1:,4] = group1
    data[1:,5] = group2
    data[1:,6] = diagnostic_residue
    data[1:,7] = sensitivity
    data[1:,8] = specificity
    data[1:,9] = F1
    #------------------------------
    # Write to tab delimited file
    #------------------------------
    
    out_file = "%s_%s_stats.csv" % (virusName,dir_pick)
    ofile = open(out_file,"wb")
    writer = csv.writer(ofile, delimiter=',')
    for jj in range(len(data)):
        writer.writerow(data[jj])

    ofile.close()


# In[ ]:



