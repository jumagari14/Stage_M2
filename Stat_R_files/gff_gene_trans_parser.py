"""Gene-mRNA-protein parser. 
This script creates gene-mRNA-protein relations from a gff and protein info file. Protein file must have the information on the 5th column, where the associated mRNA IDs are found.
Garcia, Juan Manuel
17/04/2020
"""

if __name__=="__main__":
    import argparse
    import pandas as pd 
    import re 
    import csv

    from pathlib import Path

    parser=argparse.ArgumentParser(prog='GFF gene-transcrits parser', description='Get a list of genes and their associated transcrits from a gff file and proteins from a separate file.') 

    parser.add_argument('-f', action='store', dest='gff', type=Path, help='GFF file to be analyzed')
    parser.add_argument('-protein', action='store', dest='protein_data', type=Path, help='Protein csv data to be analyzed')
    parser.add_argument('-out',action='store',dest='out_file',type=str,help='Output csv filename where the list is stored.')
    
    args=parser.parse_args()

    f_gff=open(args.gff,"r")
    lines_gff=f_gff.readlines()

    dict_genes=[]

    for line in lines_gff: 
        if "##" in line.split("\t")[0]: 
            continue
        feat=line.split("\t")[2]
        info=line.split("\t")[8]
        if feat=="gene": 
            ID_gene=re.search("(?<=ID=)[A-Za-z|0-9]*",info).group(0)
        if feat=="mRNA": 
            ID_mrna=re.search("(?<=ID=)[A-Za-z|0-9|\.]*",info).group(0)
            Parent=re.search("(?<=Parent=)[A-Za-z|0-9]*",info).group(0)
            if Parent==ID_gene:
                dict_genes.append({"Gene":ID_gene,"Transcrit":ID_mrna})

    f_gff.close()
    header= dict_genes[0].keys()
    with open('`gene_mrna_Red5`.csv', 'w',encoding='utf8',newline='') as output_file:
        dict_writer = csv.DictWriter(output_file,fieldnames=header)
        dict_writer.writeheader()
        dict_writer.writerows(dict_genes)
    list_prot=[]
    with open(args.protein_data,"r") as lines_pro: 
        next(lines_pro)
        for line in lines_pro: 
            line=line.strip("\n")
            id_prot=line.split(",")[1].replace('"','')
            info_prot=line.split(",")[4]
            gene_prot=re.search("(?<=CEY00_)Acc[0-9]+",info_prot).group(0)
            list_prot.append({"Gene": gene_prot,"Protein": id_prot, "Accession": line.split(",")[3].replace('"','')})
    lines_pro.close()
    for k in range(len(list_prot)): 
        for i in range(len(dict_genes)): 
            if list_prot[k]["Gene"]==dict_genes[i]["Gene"]: 
                dict_genes[i].update(list_prot[k])
    dict_genes=[x for x in dict_genes if "Protein" in list(x.keys())]

    header= dict_genes[0].keys()
    with open(args.out_file, 'w',encoding='utf8',newline='') as output_file:
        dict_writer = csv.DictWriter(output_file,fieldnames=header)
        dict_writer.writeheader()
        dict_writer.writerows(dict_genes)
