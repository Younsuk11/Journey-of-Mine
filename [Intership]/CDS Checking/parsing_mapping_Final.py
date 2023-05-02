# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:02:19 2021

@author: curig
"""



from tkinter import filedialog
from tkinter import *
import sys, os
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import gffpandas.gffpandas as gffpd
from Bio import Entrez, SeqIO
from time import sleep
from tqdm import tqdm
from collections import defaultdict
import datetime

def directory():
    root = Tk()
    root.dirName = filedialog.askdirectory() # directory 저장
    root.destroy()
    root.mainloop()
    
    
    return root.dirName

#gz file 압축해제
def unzip(f_n, directory):
    os.chdir(directory)
    with gzip.open(f"{f_n}.gz","rb") as f_in:
        with open(f"{f_n}", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
            
            
def Transcript_ID(x):
    num = x.count("-")
    if num == 2:
        idx = x.rfind("-")
        return x[4:idx]
        
    elif num == 1:
        return x[4:]

def find_start_end_idx(a, b):
    start = a.find(b)
    end = int(start) + len(b) -1
    return int(start) + 1, int(end) + 1


#for parsing
def find_gene(x):
    idx1 = x.find("gene")
    idx2 = x.find("db_xref")
    
    return x[idx1+5:idx2-3]

def find_ID(x):
    idx1 = 4
    idx2 = x.find("_", 7)
    seq_id = x[idx1:idx2]
    idx3 = x.find("mrna_") +5
    idx4 = x.find("_",idx3+3)
    t_id = x[idx3:idx4]
    idx5 = x.find("cds_") + 4
    idx6 = x.find("_", idx5 + 3)
    c_id =  x[idx5:idx6]
    
    return seq_id, t_id, c_id


def making_ref(x, y):
    
    if x == y:
        
        if int(x) != 0:
            return "Mapping"
    
        else:
            return "Unmapping"
        
    else:
        if int(y) == 0:
            return "NCBI-0"
        
        else: 
            return "NCBI"
        
def NCBI_search(t_id):
    
   
   Entrez.email = "yuemyoun@gmail.com"
   hd = Entrez.efetch(db = "nucleotide", id =t_id   , rettype = "gb", retmode = "text")
   recs= list(SeqIO.parse(hd, "gb"))
        
   temp = [f for f in recs[0].features if f.type == "CDS"][0]
   return int(temp.location.start +1), int(temp.location.end)
  



def Entrez_search(index_list, map_df):
    
    for i, index in enumerate(tqdm(index_list)):
        (start, end) = NCBI_search(map_df.loc[index, "Transcript_ID"])
        map_df.loc[index,"start"] = start
        map_df.loc[index,"end"] = end
        map_df.loc[index, "CDS_seq"] = map_df.loc[index, "RNA_seq"][start-1 : end]# ***
        
        

def find_UTR(RNA_seq, CDS_start, CDS_end):
    return RNA_seq[:int(CDS_start)-1], RNA_seq[int(CDS_end):]


def parsing_fasta(file_name):


   
    fasta_di = {}
    fasta_di["rna_g"] = ""
    fasta_di["cds_g"] = ""
    fasta_di["rna"] = ""
    
    for key in fasta_di.keys():
    
        seq = file_name[key]
    
        
        ID = []
        description = []
        sequence = []
    
        for i,seq_record in enumerate(tqdm(SeqIO.parse(seq,"fasta"))):
           
            
            ID.append(str(seq_record.id))
            
            if key == "rna":
                temp = seq_record.description.split(",")[0]
                t = temp.split(" ")[1:-1]
                description.append(" ".join(t))
            else:
                description.append(seq_record.description)
                
            sequence.append(seq_record.seq)
        
        
        fasta_di[key] = pd.DataFrame({"ID": ID, "Description":description, "Sequence": sequence})
        
        if key == "rna":
            fasta_di[key].columns = ["Transcript_ID", "Description", "Sequence"]
            fasta_di[key]["Tag"] = fasta_di[key]["Transcript_ID"].map(lambda x: x[:3])
            fasta_di[key] = fasta_di[key].loc[fasta_di[key]["Tag"].isin(["NM_","XM_"])]
            fasta_di[key].reset_index(drop = True, inplace = True)
        
        if key == "rna_g":
            fasta_di[key] = fasta_di[key].loc[fasta_di[key]["ID"].str.contains("_mrna_")]
            fasta_di[key].reset_index(drop = True, inplace = True)
            fasta_di[key]["Gene"] = fasta_di[key]["Description"].map(lambda x: find_gene(x))
            fasta_di[key]["seq_id"] =  fasta_di[key]["ID"].map(lambda x : find_ID(x)[0])
            fasta_di[key]["Transcript_ID"] = fasta_di[key]["ID"].map(lambda x : find_ID(x)[1])
            fasta_di[key]["Tag"] = fasta_di[key]["Transcript_ID"].map(lambda x: x[:3])
            fasta_di[key] = fasta_di[key].loc[fasta_di[key]["Tag"].isin(["NM_","XM_"])]
            fasta_di[key].reset_index(drop = True, inplace = True)
            
            
            
        elif key == "cds_g":
            
            fasta_di[key]["Gene"] = fasta_di[key]["Description"].map(lambda x: find_gene(x))
            fasta_di[key]["seq_id"] =  fasta_di[key]["ID"].map(lambda x : find_ID(x)[0])
            fasta_di[key]["Protein_ID"] = fasta_di[key]["ID"].map(lambda x : find_ID(x)[2])
            fasta_di[key]["Tag"] = fasta_di[key]["Protein_ID"].map(lambda x: x[:3])
            fasta_di[key] = fasta_di[key].loc[fasta_di[key]["Tag"].isin(["NP_","XP_"])]
            fasta_di[key].reset_index(drop = True, inplace = True)
            
            
        
        fasta_di[key]["Length"] = fasta_di[key]["Sequence"].map(lambda x : len(x))
        
    return fasta_di



    
def annotation_file(file_name):

    
    annotation = gffpd.read_gff3(file_name)

    #어차피 cds annotation에서 parent(rna) 매칭 되어있음
    cds_ref = annotation.filter_feature_of_type(["CDS"]).attributes_to_columns()## no NC_filtering - id_116250
    cds_ref.dropna(subset = ["Name"], inplace = True)
    i = cds_ref.index[cds_ref["Name"].str.contains("YP_")].tolist()
    cds_ref.drop(i, inplace = True)
    ####전체 위치 가져오기.
    
    #seq_id, gene, Parent, Name, 모두 같은데 서열 위치가 다른 경우가 있음.(각각 서열이 조금씩 다르다는 이야기)
    #일단, rna-cds matching을 위해 모두 같은 것들은 중복 제거함. -- seq_id에 따라 rna id 달라질것임.
    cds_ref = cds_ref[["seq_id","gene","Parent", "Name", "strand" ,"tag"]]
    cds_ref.drop_duplicates(["seq_id","gene","Parent", "Name" ], keep = "first", inplace = True) ####
    
    num_cds = len(cds_ref["Name"].drop_duplicates())
    g_cds = cds_ref["gene"].drop_duplicates()
    
    cds_ref.sort_values(by = "Name", ascending = True, inplace = True)
    cds_ref.reset_index(drop = True, inplace = True)
    
    #cds vs rna checking : 1.gene pool check 2.seq number check
    rna_ref = annotation.filter_feature_of_type(["mRNA"]).attributes_to_columns()
    rna_ref.dropna(subset = ["Name"], inplace = True)
    rna_ref.sort_values(by = "seq_id", ascending = True, inplace = True)
    rna_ref.reset_index(drop = True, inplace = True)
    num_rna = len(rna_ref["Name"].drop_duplicates())
    g_rna = rna_ref["gene"].drop_duplicates()
    
    if num_cds == num_rna:
        print(f"Yes, there is full match in cds annotation : {num_cds} sequences")
    else:
        print("No. check again")
        
    if len(g_cds) == len(g_rna):
        print("Yes, genes in both annotation are also matched! : ",len(g_cds))
    else: 
        print("No full match, check again!")
        
        
    cds_ref.insert(4, "Transcript_ID", cds_ref["Parent"].map(lambda x : Transcript_ID(x)))
    cds_ref.columns = ["seq_id", "Gene", "Parent", "Protein_ID", "Transcript_ID","pos_neg","tag"]
    
    return cds_ref





def RNA_CDS_Mapping(rna, cds,rna_seq, annotation):
    
  
    raw = {}
    raw["cds"] = cds[["Gene", "seq_id", "Protein_ID", "Sequence"]]
    raw["rna"] = rna[["Gene", "seq_id", "Transcript_ID", "Sequence"]]
    
    
    col = {}
    col["cds"] = ["Gene","seq_id",  "Protein_ID"]
    col["rna"] = ["Gene","seq_id",  "Transcript_ID"]
    
    
    mapping = annotation[["Gene", "seq_id", "Transcript_ID", "Protein_ID", "pos_neg","tag"]]
    
    
    for key in col.keys():
       
        
        mapping.sort_values(by = col[key], ascending = True,  ignore_index = True, inplace = True)
        
        raw[key].sort_values(by = col[key], ascending = True, ignore_index = True, inplace = True)
    
        
        df = pd.concat([mapping[col[key]], raw[key][col[key]]])
        df_grp = df.groupby(col[key])
        df_di = df_grp.groups
        idx = [x[1] for x in df_di.values() if len(x) == 2] #annotaion에 맞게 행 재배열 인덱스
        
        raw_rows = raw[key].loc[idx].reset_index(drop = True)
        
   
        
        if key == "rna":
            
            mapping["Full_seq"]= raw_rows["Sequence"]
        else:
            mapping["CDS_seq"] = raw_rows["Sequence"]
            
            
        

    #Length 계산 및 Full_seq에서 CDS_seq 위치 mapping(index 시작은 1)
    mapping["Full_length"]= mapping["Full_seq"].map(lambda x: len(x))
    mapping["CDS_length"]= mapping["CDS_seq"].map(lambda x: len(x))
    mapping["CDS_start"] = mapping.apply(lambda x : find_start_end_idx(x["Full_seq"], x["CDS_seq"])[0], axis = 1)
    mapping["CDS_end"] = mapping.apply(lambda x : find_start_end_idx(x["Full_seq"], x["CDS_seq"])[1], axis = 1)
    
 
    
    #vlookup
    mapping["Description"] = mapping["Transcript_ID"].map(rna_seq.set_index("Transcript_ID")["Description"])
    mapping["RNA_seq"] = mapping["Transcript_ID"].map(rna_seq.set_index("Transcript_ID")["Sequence"])
    
    #mapping = mapping.join(rna_seq.set_index("Transcript_ID")["Sequence"], on = "Transcript_ID)
    mapping["RNA_length"] = mapping["RNA_seq"].map(lambda x : len(x))
    mapping["start"] = mapping.apply(lambda x : find_start_end_idx(x["RNA_seq"], x["CDS_seq"])[0], axis = 1)
    mapping["end"] = mapping.apply(lambda x : find_start_end_idx(x["RNA_seq"], x["CDS_seq"])[1], axis = 1)
    mapping["REF"] = mapping.apply(lambda x: making_ref(x["CDS_start"], x["start"]), axis = 1)
    
    

    
    # mapping 안된것들 다시 mapping = NCBI data base 
    Unmapping = list(mapping[mapping["REF"] == "Unmapping"]["Transcript_ID"].index)
    N_0 =list(mapping[mapping["REF"] == "NCBI-0"]["Transcript_ID"].index)


    Entrez_search(Unmapping, mapping)
    Entrez_search(N_0, mapping)
    
    mapping = mapping[["Gene", "Description", "seq_id", "tag","Transcript_ID", "Protein_ID", "pos_neg", "RNA_seq", "CDS_seq", "RNA_length", "REF", "start","end"]]
    
        
    
    mapping["5_UTR"] = mapping.apply(lambda x: find_UTR(x["RNA_seq"], x["start"], x["end"])[0], axis = 1)
    mapping["3_UTR"] = mapping.apply(lambda x: find_UTR(x["RNA_seq"], x["start"], x["end"])[1], axis = 1)
    

    
    mapping.sort_values(by = ["Gene", "seq_id", "Transcript_ID", "Protein_ID"], ascending = True, ignore_index = True, inplace = True)
    
    return mapping








file_name = {}
file_name["annotation"] = "GCF_000001405.39_GRCh38.p13_genomic.gff"
file_name["rna_g"] = "GCF_000001405.39_GRCh38.p13_rna_from_genomic.fna" 
file_name["cds_g"] = "GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna" 
file_name["rna"] = "GCF_000001405.39_GRCh38.p13_rna.fna" 
d = datetime.datetime.now()

start = True
while(start):
    
    cmd = str(input("프로그램을 시작하시겠습니까? : (start or end)\n")).lower()
    if cmd == "start":
        start = True
        print("\"Let me know where the \"Raw data files\" are...!\"\n")
        pathway = directory()
        os.chdir(pathway)
            
        a = True
        while(a):
            cmd2 = str(input("원하시는 기능을 고르세요 : 1.parsing_fasta | 2.making_annotation | 3.mapping | 4.all_at_once |end\n")).lower()
            
            if cmd2 =="1":
                print("                     \"I AM parsing the fasta files\"\n")
                fasta_di = parsing_fasta(file_name)
                a = True
                
            elif cmd2 == "2":
                print("\"I AM making annotation file\"\n")
                annotation = annotation_file(file_name["annotation"])
                a = True
                
            elif cmd2 =="3":
                print("\"I AM Mapping the rna fasta file and the cds file with annotation file\"\n")
                mapping = RNA_CDS_Mapping(fasta_di["rna_g"], fasta_di["cds_g"], fasta_di["rna"], annotation)
                print("\"I AM saving the completed mapped file\"\n")
                
                mapping.to_csv("RNA_CDS({}_{}_{}).csv".format(d.year, d.month, d.day), sep = ",", index = False)
                a = True
                
            elif cmd2 == "4":
                print("\"I AM Doing all Process\"\n")
                fasta_di = parsing_fasta(file_name)
                annotation = annotation_file(file_name["annotation"])
                mapping = RNA_CDS_Mapping(fasta_di["rna_g"], fasta_di["cds_g"], fasta_di["rna"], annotation)
                
                mapping.to_csv("RNA_CDS({}_{}_{}).csv".format(d.year, d.month, d.day), sep = ",", index = False)
                a = False
                start = False
                break
            
            elif cmd2 == "end":
                a = False
                start = False
                break
            
            else:
                print("정확한 명령 숫자를 입력해주세요....\n")
                a = True
            
        
    elif cmd == "end":
        start = False
        break
    
    else: 
        print("명령어를 다시 입력해주세요 : (start or end)\n")
        start = True
        








