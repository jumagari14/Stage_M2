"""Mercator parser
Garcia, Juan Manuel
Extract first-level definitions and their associated gene IDs obtained from Mercator
17/04/2020
"""
if __name__=="__main__": 
    import re 
    import argparse
    import fileinput
    import csv
    import pandas as pd
    from pathlib import Path

    parser=argparse.ArgumentParser(prog="MapMan parser", description="Small script that extract genes related to every first-level annotation")

    parser.add_argument('-r',action='store',dest='fichier_r',type=Path,help='File to be read')
    parser.add_argument('-w',action='store',dest='newfile',type=Path,help='New csv file')

    args=parser.parse_args()
    parsed={}

    with open(args.fichier_r,'r') as f: 
        next(f)
        for line in f: 
            line=line.strip("\n")
            line=line.split("\t")
            if re.search("'[0-9]+\.",line[0])==None: 
                bin_id=re.search("'[0-9]+",line[0]).group(0)
                annot=line[1].replace("'","")
                parsed[annot]=[]
            else: 
                string=line[2].replace("'","").capitalize()
                if (re.search("'[0-9]+",line[0]).group(0)==bin_id) and (string not in parsed[annot]) and (string !=""): 
                    parsed[annot].append(string)
    f.close()
    # data=pd.DataFrame((parsed),columns=list(parsed.keys()))
    data=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in parsed.items() ]))
    # print(data)
    data.to_csv(args.newfile,sep=",")
