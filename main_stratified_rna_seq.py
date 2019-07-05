#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3

# Author : Debojyoti Das
# Date   : 28th May, 2019

import sys
import time
import pandas as pd
import numpy as np
import statsmodels.stats.contingency_tables as sm

start = time.time()

filename_in = str(sys.argv[1])
file_in = open(filename_in, mode = 'rt')
# output file
file_parts = filename_in.split("/")
filename_out = "sorted_data_frame_" + str(file_parts[len(file_parts)-1])
file_out = open(filename_out, mode = 'wt')

line = file_in.readline()
columns = line[:len(line)-1].split(",")
#print(len(columns))

# data.frame containing data variant counts and meta-data
data = pd.DataFrame()

# Number of strata
sk = 10
# same_row = reference and allele counts in the same line/row
#if(str(sys.argv[1]) == "same_row"):
number_of_elements_of_stratified_table = 2 * 2 * sk

for line in file_in:
    l = list() 
    if (line.find('"') != -1):
        l_list = line[:len(line)-1].split('"')
        l0 = l_list[0]
        l0 = l0[:len(l0)-1].split(',')
        l.extend(l0)
        l1 = l_list[1]
        l1 = l1[1:len(l1)]
        l.append(l1)
        l2 = l_list[2]
        l2 = l2[1:].split(',')
        l.extend(l2)
    else :
        l_list = line[:len(line)-1].split(',')
        l.extend(l_list)
    l = l[1:]
    row = pd.DataFrame(l).T
    data = data.append(row, ignore_index=True)
#   print(l)

data.columns = columns[1:]
#print(data.columns)

currTemp = data.iloc[0,0]
currDay = data.iloc[0,1]
currPop = data.iloc[0,2]
currMapping = data.iloc[0,3]
currPosition = data.iloc[0,4]
currType = data.iloc[0,5]
currRef = data.iloc[0,6]
currAllele = data.iloc[0,7]
currCount = data.iloc[0,12]
nrow = data.shape[0]
currRow = 0
'''
print("Population         : " + currPop)
print("Chromosome/ gene   : " + currMapping)
print("loci/ region       : " + currPosition)
print("Type of variant    : " + currType) 
print("Reference Allele   : " + currRef) 
print("ALternate Allele   : " + currAllele) 
print("Allele Count       : " + currCount)
print("Rows in data.frame : " + str(nrow))
#print(data.shape)
#print(data)
'''

data.sort_values(by = ['Mapping', 'Reference.Position', 'Temp', 'Day'], inplace = True)
data.to_csv(file_out)

while(currRow < nrow):
    matchindex = currRow
    while((matchindex != nrow) and (data.iloc[matchindex,3] == currMapping) 
            and (data.iloc[matchindex,4] == currPosition)):
        matchindex = matchindex + 1
#
    if ((matchindex - currRow) == number_of_elements_of_stratified_table):
        ctable = np.zeros(shape=(2,2,sk))
        ict = currRow 
        for k in range(0,sk):
#           print("=======================================")
#           print("strata :" + str(k+1))
#           print("=======================================")
            for i in range(0,2):
                for j in range(0,2):
                    ctable[i,j,k] = data.iloc[ict,12]
#                   print(" c(" + str(i+1) + "," + str(j+1) + "," + str(k+1) + ") = " + str(ctable[i,j,k]), end = '')
                    ict = ict + 1
#               print("\n", end = '')
#           print(ctable[:,:,k])
        cmh = sm.StratifiedTable(ctable)
        t = cmh.test_null_odds()
        print(currTemp + " " + currDay + " " + currMapping + " " + currPosition + "\t" + currRef + "\t" + currAllele + "\t" + str(t.pvalue) + "\t"+ str(t.statistic))
#       print(cmh.test_null_odds())
#       print(cmh.test_equal_odds())
#       print(cmh.summary())
#   while loop check    
    currRow = matchindex
    if(currRow < nrow):
        currTemp = data.iloc[currRow,0]
        currDay = data.iloc[currRow,1]
        currPop = data.iloc[currRow,2]
        currMapping = data.iloc[currRow,3]
        currPosition = data.iloc[currRow,4]
        currType = data.iloc[currRow,5]
        currRef = data.iloc[currRow,6]
        currAllele = data.iloc[currRow,7]
        currCount = data.iloc[currRow,12]
    else:
        print(matchindex)
        break

print('It took', time.time()-start, 'seconds.')
