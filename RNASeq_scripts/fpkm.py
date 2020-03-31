#!/usr/bin/python3.6
"""
Juan Manuel Garcia 
Convert FPKM to TPM 
Based on equations explained in https://arxiv.org/pdf/1104.3889.pdf

"""
## Functions 
def getTPM(filepath, filename,folder,count_file): 
    data=pd.read_table(filepath)
    count_data=pd.read_table(count_file,header=None,sep="\t",names=['tracking_id','count'])
    for i in range(data.shape[0]): 
        start_pos=re.search('[A-Za-z|0-9]*:([0-9]*)-([0-9]*)',str(data.iloc[i]["locus"])).group(1)
        end_pos=re.search('[A-Za-z|0-9]*:([0-9]*)-([0-9]*)',str(data.iloc[i]["locus"])).group(2)
        longueur=int(end_pos)-int(start_pos)+1
        count_data.loc[count_data["tracking_id"]==data.iloc[i]["tracking_id"],"Length"]=longueur
        count_data.loc[count_data["tracking_id"]==data.iloc[i]["tracking_id"],"Norm"]=count_data.loc[count_data["tracking_id"]==data.iloc[i]["tracking_id"],"count"]/longueur
    sum_norm=count_data.sum()["Norm"]
    # TPM=[]
    # tpm=list(count_data.apply(lambda line: (line["Norm"]/sum_norm)*10**6,axis=1))
    # TPM.extend(tpm)
    TPM=[]
    for i in range(count_data.shape[0]): 
        tpm=(count_data.iloc[i]["Norm"]/sum_norm)*10**6
        TPM.append(tpm)

    TPM_data=pd.DataFrame(data={"tracking_id": count_data["tracking_id"], "TPM": TPM})
    TPM_data.drop(TPM_data.tail(5).index,inplace=True)
    header=filename.split(".")[0]
    path=filepath.rsplit("/",2)[1].replace("FPKM","/tpm")
    TPM_data.to_csv(folder+"/"+path+"/"+header+"_TPM.csv",mode="w",index=False)
    #TPM_data.to_csv

## MAIN


import argparse
import pandas as pd 
import numpy as np
import os 
import re

parser=argparse.ArgumentParser(prog='TPM calculator', description='This script extracts the TPM values from a tab-delimited file based on FPKM values') 

parser.add_argument('-d', action='store', dest='pathdir', type=str, help='Path to directory where gene and isoform measurements are stored in one folder per each sample')
args=parser.parse_args()

dirpath=args.pathdir

new_direc=dirpath.split("/")[0]+"/tpm"
print(new_direc)
if os.path.isdir(new_direc)==False: 
    os.makedirs(new_direc)
for root,dirs,files in os.walk(dirpath): 
    for direc in dirs: 
        # print(direc)
        # direc=direc.replace("fpkm","/tpm/")
        if os.path.isdir(new_direc+direc)==False: 
            os.makedirs(new_direc+"/"+direc)
    if (root[len(dirpath):].count(os.sep)<3 and root[len(dirpath):]!=""):
        print(root[len(dirpath):])
        count_filename=root[len(dirpath):].replace("FPKM_","HTSeq_") 
        count_file=dirpath.split("/")[0]+"/readCount/HTSeq_"+count_filename+".txt"  
        print(count_file)
        for f in files: 
            if (f=="isoforms.fpkm_tracking"): 
                filepath=os.path.join(root,f)
                getTPM(filepath,f,new_direc,count_file)

#getTPM("resultsConcombre/fpkm/FPKM_5012_CCGCGGTT-AGCGCTAG-BHKJVVDSXX_L003/genes.fpkm_tracking","genes.fpkm_tracking","resultsConcombre/tpm","resultsConcombre/readCount/readCount-5012_CCGCGGTT-AGCGCTAG-BHKJVVDSXX_L003.txt")