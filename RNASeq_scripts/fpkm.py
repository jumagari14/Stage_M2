#!/usr/bin/python3.6
"""
Juan Manuel Garcia 
Convert FPKM to TPM 
Based on equations explained in https://arxiv.org/pdf/1104.3889.pdf

"""
## Functions 
def getTPM(filepath, filename,folder): 
    data=pd.read_table(filepath)

    sum_fpkm=data.sum().iloc[-4]
    TPM=[]
    for i in range(data.shape[0]): 
        tpm=(data.iloc[i]["FPKM"]/sum_fpkm)*10**6
        TPM.append(tpm)

    TPM_data=pd.DataFrame(data={"tracking_id": data["tracking_id"], "TPM": TPM})
    header=filename.split(".")[0]
    path=filepath.rsplit("/",2)[1]
    TPM_data.to_csv(folder+path+"/"+header+"_TPM.csv",mode="w",index=False)
    #TPM_data.to_csv

## MAIN


import argparse
import pandas as pd 
import numpy as np
import os 

parser=argparse.ArgumentParser(prog='TPM calculator', description='This script extracts the TPM values from a tab-delimited file based on FPKM values') 

parser.add_argument('-d', action='store', dest='pathdir', type=str, help='Path to directory where gene and isoform measurements are stored in one folder per each sample')
args=parser.parse_args()

dirpath=args.pathdir

new_direc=dirpath.split("/")[0]+"tpm"
if os.path.isdir(new_direc)==False: 
    os.makedirs(new_direc)
for root,dirs,files in os.walk(dirpath): 
    for direc in dirs: 
        if os.path.isdir(new_direc+direc)==False: 
            os.makedirs(new_direc+direc)
    if root[len(dirpath):].count(os.sep)<2:       
        for f in files: 
            if (f=="genes.fpkm_tracking" or f=="isoforms.fpkm_tracking"): 
                filepath=os.path.join(root,f)
                getTPM(filepath,f,new_direc)