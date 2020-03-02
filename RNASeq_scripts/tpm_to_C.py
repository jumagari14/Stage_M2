#!/usr/bin/python3
"""
Juan Manuel Garcia 
Calculate concentration of each transcript from TPM values 
Based on equations explained in https://arxiv.org/pdf/1104.3889.pdf

"""



import argparse
import pandas as pd 
import numpy as np
import os 


def tpmToC(filepath,filename,exp_values): 
    data=pd.read_csv(filepath)
    suma=data.loc[:,2].sum()
    #for i in 
    return True

parser=argparse.ArgumentParser(prog='Calculate concentration (fmol/g F.W) based on TPM values', description='This script performs the calculation of concentration considering the 5 TPM values of the spikes.') 

parser.add_argument('-d', action='store', dest='pathdir', type=str, help='Path to directory where gene and isoform tpm values are stored in one subdirectory per each sample')
args=parser.parse_args()

dirpath=args.pathdir

## Load data about spikes obtained experimentally 

path_file1=pd.ExcelFile("dataKentaro/samples concentration.xlsx")

kiwi_data=pd.read_excel(path_file1,"Kiwifruit")
conc_data=pd.read_excel(path_file1,"Concombre")


spike_mol={'M8': [[2.574E-5,1.901E-4,1.259E-3,1.043E-2,9.931E-2,8.108E-1,5.132,37.44],[9.975e-06, 7.304e-05, 0.0004923999999999999, 0.004072, 0.03801, 0.3144, 1.99, 14.540000000000001],[1.022e-05, 7.508e-05, 0.0005027, 0.0041600000000000005, 0.03937, 0.3221, 2.04, 14.89]], 
'M5': [[5.148e-05, 0.0003802, 0.002518, 0.02086, 73.46, 0.08108, 0.687, 3.7809999999999997],[1.995e-05, 0.0001461, 0.0009847999999999999, 0.008144, 28.12, 0.031439999999999996, 0.2664, 1.468],[7.661e-06, 7.633e-05, 0.0008072, 0.008157999999999999, 29.12, 0.032209999999999996, 0.2731, 1.503]],
'M6': [[5.148e-05, 0.0003802, 0.002518, 0.02086, 0.1986, 59.98, 0.5061, 5.132],[1.995e-05, 0.0001461, 0.0009847999999999999, 0.008144, 0.07601000000000001, 23.26, 0.19649999999999998, 1.99],[7.661e-06, 7.633e-05, 0.0008072, 0.008157999999999999, 0.08030000000000001, 23.830000000000002, 0.2012, 2.04]],
'M7': [[2.574e-05, 0.0001901, 0.001259, 0.01043, 0.09931000000000001, 0.8108, 50.82, 3.7809999999999997],[9.975e-06, 7.304e-05, 0.0004923999999999999, 0.004072, 0.03801, 0.3144, 19.709999999999997, 1.468],[3.831e-06, 3.816e-05, 0.0004036, 0.004078999999999999, 0.040150000000000005, 0.4012, 39.71, 4.013]]}



os.makedirs("Conc")
for root,dirs,files in os.walk(dirpath): 
    for direc in dirs: 
        os.makedirs("Conc/"+direc)
    if root[len(dirpath):].count(os.sep)<2:       
        for f in files: 
            if ("genes" in f or "isoforms" in f): 
                filepath=os.path.join(root,f)
                values=0
                tpmToC(filepath,f,values)



