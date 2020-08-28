"""
Garcia, Juan Manuel
Convert FPKM to TPM 
Based on equations explained in https://arxiv.org/pdf/1104.3889.pdf

"""


## Functions 
def getTPM(filepath, filename,folder,count_file): 
    data=pd.read_table(filepath)
    count_data=pd.read_table(count_file,header=None,sep="\t",names=['tracking_id','count'])
    data["longueur"]=[int(re.search('[A-Za-z|0-9]*:([0-9]*)-([0-9]*)',str(x)).group(2)) - int(re.search('[A-Za-z|0-9]*:([0-9]*)-([0-9]*)',str(x)).group(1))+1 for x in data["locus"]]
    count_data=pd.merge(data,count_data,on="tracking_id")
    count_data["Norm"]=count_data["count"]/count_data["longueur"]
    sum_norm=count_data.sum()["Norm"]
    count_data["TPM"]=(count_data["Norm"]/sum_norm)*10**6
    TPM_data=pd.DataFrame(data={"tracking_id": count_data["tracking_id"], "TPM": count_data["TPM"]})
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

# parser=argparse.ArgumentParser(prog='TPM calculator', description='This script extracts the TPM values from a tab-delimited file based on FPKM values') 

# parser.add_argument('-d', action='store', dest='pathdir', type=str, help='Path to directory where gene and isoform measurements are stored in one folder per each sample')
# args=parser.parse_args()

dirpath="Red5Kiwi/fpkm/"

new_direc=dirpath.split("/")[0]+"/tpm"
print(new_direc)
if os.path.isdir(new_direc)==False: 
    os.makedirs(new_direc)
for root,dirs,files in os.walk(dirpath): 
    for direc in dirs: 
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
