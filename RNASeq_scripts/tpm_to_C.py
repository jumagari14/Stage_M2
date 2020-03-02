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
from skmisc.loess import loess


def tpmToC(filepath,filename,exp_values): 
    data=pd.read_csv(filepath)
    y=exp_values
    x=list(data.loc[data['tracking_id'].str.contains('spike'),"TPM"]) ## Second column 
    #suma=data.loc[:,2].sum()
    l = loess(x,y)
    l.fit()
    x_pred=list(data.loc[data['tracking_id'].str.contains("."),"TPM"])
    pred = l.predict(x_pred, stderror=True) ## Extrapolation problem
    conf = pred.confidence()

    lowess = pred.values
    ll = conf.lower
    ul = conf.upper
    #for i in 
    return {'Column name': data['tracking_id'], 'Concentration (fmol)': pred}

parser=argparse.ArgumentParser(prog='Calculate concentration (fmol/g F.W) based on TPM values', description='This script performs the calculation of concentration considering the 5 TPM values of the spikes.') 

parser.add_argument('-d', action='store', dest='pathdir', type=str, help='Path to directory where gene and isoform tpm values are stored in one subdirectory per each sample')
args=parser.parse_args()

dirpath=args.pathdir

## Load data about spikes obtained experimentally 

path_file1=pd.ExcelFile("./dataKentaro/samples concentration.xlsx")

kiwi_data=pd.read_excel(path_file1,"Kiwifruit")
conc_data=pd.read_excel(path_file1,"Concombre")


spike_mol={'M8': [[2.574E-5,1.901E-4,1.259E-3,1.043E-2,9.931E-2,8.108E-1,5.132,37.44],[9.975e-06, 7.304e-05, 0.0004923999999999999, 0.004072, 0.03801, 0.3144, 1.99, 14.540000000000001],[1.022e-05, 7.508e-05, 0.0005027, 0.0041600000000000005, 0.03937, 0.3221, 2.04, 14.89]], 
'M5': [[5.148e-05, 0.0003802, 0.002518, 0.02086, 73.46, 0.08108, 0.687, 3.7809999999999997],[1.995e-05, 0.0001461, 0.0009847999999999999, 0.008144, 28.12, 0.031439999999999996, 0.2664, 1.468],[7.661e-06, 7.633e-05, 0.0008072, 0.008157999999999999, 29.12, 0.032209999999999996, 0.2731, 1.503]],
'M6': [[5.148e-05, 0.0003802, 0.002518, 0.02086, 0.1986, 59.98, 0.5061, 5.132],[1.995e-05, 0.0001461, 0.0009847999999999999, 0.008144, 0.07601000000000001, 23.26, 0.19649999999999998, 1.99],[7.661e-06, 7.633e-05, 0.0008072, 0.008157999999999999, 0.08030000000000001, 23.830000000000002, 0.2012, 2.04]],
'M7': [[2.574e-05, 0.0001901, 0.001259, 0.01043, 0.09931000000000001, 0.8108, 50.82, 3.7809999999999997],[9.975e-06, 7.304e-05, 0.0004923999999999999, 0.004072, 0.03801, 0.3144, 19.709999999999997, 1.468],[3.831e-06, 3.816e-05, 0.0004036, 0.004078999999999999, 0.040150000000000005, 0.4012, 39.71, 4.013]]}

"""
if os.path.isdir("Conc")==False: 
    os.makedirs("Conc")
for root,dirs,files in os.walk(dirpath): 
    concen=[]
    cont=0
    for direc in dirs: 
        cont+=1
        iden=direc.split("_")[1]
        mix_n=list(conc_data.loc[conc_data["Sample"]==float(iden),"Mix"])[0]
        if os.path.isdir("Conc/"+direc): 
            os.makedirs("Conc/"+direc)
    if root[len(dirpath):].count(os.sep)<2:       
        for f in files: 
            if "genes" in f: 
                filepath=os.path.join(root,f)
                mol_values=spike_mol[mix_n]
                concen.append(tpmToC(filepath,f,mol_values))
        if cont%3==0: 
            print(concen)
            print("Jaja")"""
            ## rowSums of conc values from 3 consecutive dictionnaries 
            ## Set concen list as empty 

conc=tpmToC("/run/user/1000/gvfs/sftp:host=curta.mcia.fr,user=juagarcia/gpfs/home/juagarcia/git_repo/RNASeq_scripts/tpm/FPKM_5012_CCGCGGTT-AGCGCTAG-BHKJVVDSXX_L003/genes_TPM.csv","genes_TPM.csv",spike_mol['M5'][0])
print(conc)