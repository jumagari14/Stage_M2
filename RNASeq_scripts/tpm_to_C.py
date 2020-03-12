#!/usr/bin/python3.6
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
from sklearn.linear_model import LinearRegression


def tpmToC(filepath,filename,exp_values): 
    data=pd.read_csv(filepath)
    y=exp_values
    x=list(data.loc[data['tracking_id'].str.contains('spike'),"TPM"]) ## Second column 
    # suma=data.loc[:,2].sum()
    x=np.array(x).reshape(-1,1)
    y=np.array(y).reshape(-1,1)
    reg = LinearRegression().fit(x,y)
    if reg.score(x,y)!=0: 
        x_pred=data.loc[~data['tracking_id'].str.contains("spike"),"TPM"]
        x_pred=np.array(x_pred).reshape(-1,1)
        pred=reg.predict(x_pred)
        return pd.DataFrame({'Column name': data['tracking_id'][~data['tracking_id'].str.contains("spike")], 'Concentration (fmol)': pred.flatten()})
    else : 
        return False
    

parser=argparse.ArgumentParser(prog='Calculate concentration (fmol/g F.W) based on TPM values', description='This script performs the calculation of concentration considering the 5 TPM values of the spikes.') 

parser.add_argument('-d', action='store', dest='pathdir', type=str, help='Path to directory where gene and isoform tpm values are stored in one subdirectory per each sample')
parser.add_argument('-f', action='store', dest='fruit', type=str, help='Fruit to be analyzed. This value is equal to the sheet of the excel file from which data is obtained')
args=parser.parse_args()

dirpath=args.pathdir

## Load data about spikes obtained experimentally 

path_file1=pd.ExcelFile("./dataKentaro/samples concentration.xlsx")


fruit_data=pd.read_excel(path_file1,args.fruit)

if args.fruit=="Concombre": 
    fw_ind="Weight extraction"
else: 
    fw_ind="Weigth for extraction "

spike_mol={'M8': [[2.574E-5,1.901E-4,1.259E-3,1.043E-2,9.931E-2,8.108E-1,5.132,37.44],[9.975e-06, 7.304e-05, 0.0004923999999999999, 0.004072, 0.03801, 0.3144, 1.99, 14.540000000000001],[1.022e-05, 7.508e-05, 0.0005027, 0.0041600000000000005, 0.03937, 0.3221, 2.04, 14.89]], 
'M5': [[5.148e-05, 0.0003802, 0.002518, 0.02086, 73.46, 0.08108, 0.687, 3.7809999999999997],[1.995e-05, 0.0001461, 0.0009847999999999999, 0.008144, 28.12, 0.031439999999999996, 0.2664, 1.468],[7.661e-06, 7.633e-05, 0.0008072, 0.008157999999999999, 29.12, 0.032209999999999996, 0.2731, 1.503]],
'M6': [[5.148e-05, 0.0003802, 0.002518, 0.02086, 0.1986, 59.98, 0.5061, 5.132],[1.995e-05, 0.0001461, 0.0009847999999999999, 0.008144, 0.07601000000000001, 23.26, 0.19649999999999998, 1.99],[7.661e-06, 7.633e-05, 0.0008072, 0.008157999999999999, 0.08030000000000001, 23.830000000000002, 0.2012, 2.04]],
'M7': [[2.574e-05, 0.0001901, 0.001259, 0.01043, 0.09931000000000001, 0.8108, 50.82, 3.7809999999999997],[9.975e-06, 7.304e-05, 0.0004923999999999999, 0.004072, 0.03801, 0.3144, 19.709999999999997, 1.468],[3.831e-06, 3.816e-05, 0.0004036, 0.004078999999999999, 0.040150000000000005, 0.4012, 39.71, 4.013]]}


new_direc=dirpath.split("/")[0]+"/Conc"
if os.path.isdir(new_direc)==False: 
    os.makedirs(new_direc)

cont=-1
save_dir=[]
for root,dirs,files in os.walk(dirpath): 
    dirs.sort()
    if len(dirs)!=0: 
        for direc in dirs: 
            direc=direc.split("_",1)[1]
            save_dir.append(direc)
            if os.path.isdir(new_direc+"/"+direc)!=True: 
                os.makedirs(new_direc+"/"+direc)
    if len(files)!=0:  
        cont=cont+1
        iden=root.split("_")[1]
        mix_n=list(fruit_data.loc[fruit_data["Sample"]==float(iden),"Mix"])[0]  
        fw= list(fruit_data.loc[fruit_data["Sample"]==float(iden),fw_ind])[0] 
        for f in files: 
            if "genes" in f: 
                filepath=os.path.join(root,f)
                mol_values=spike_mol[mix_n]
                concen=tpmToC(filepath,f,mol_values[cont%3])
                if isinstance(concen,pd.DataFrame) : 
                    #concen["Concentration (fmol)"]=concen["Concentration (fmol)"].flatten()
                    concen["Concentration (fmol)"]=concen["Concentration (fmol)"]/(float(fw)/1000)
                    concen.to_csv(new_direc+"/"+save_dir[cont]+"/genes.csv",mode="w",index=False)
            ## rowSums of conc values from 3 consecutive dictionnaries 
            ## Set concen list as empty 

# concen=tpmToC("resultsConcombre/tpm/tpm_5012_CCGCGGTT-AGCGCTAG-BHKJVVDSXX_L003/genes_TPM.csv","genes_TPM.csv",spike_mol["M5"][0])
# print(concen)