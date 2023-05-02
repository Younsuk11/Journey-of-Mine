# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 20:13:47 2021

@author: user
"""

import pandas as pd
import os
from tkinter import filedialog
from tkinter import *
from itertools import combinations
import time
import selenium
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import shutil 
import openpyxl
import tqdm

def save_enrichr_txt(webdriver_dir,key,original, pathways, pathway_dir, download_dir): #UP_DOWN original에서 webcrolling으로 pathway txt파일들 가져오기
    

    url = "https://maayanlab.cloud/Enrichr/"
    
    os.chdir(webdriver_dir)
    
    # options = webdriver.ChromeOptions()
    # # options.headless =True
    
    driver = webdriver.Chrome(executable_path = webdriver_dir+"/chromedriver.exe")
    
    driver.get(url = url)
    
    elem = driver.find_element_by_name("list")
        
    #gene list 복사해서 붙여넣기
    original.to_clipboard(index = False, header = False)
    elem.send_keys(Keys.CONTROL + "v")
        
    driver.find_element_by_class_name("proceed-button").click()
        
    for pathway in pathways:
       
        elem = WebDriverWait(driver,30).until(EC.presence_of_element_located((By.XPATH,"//*[@id='Pathways-link']")))
        elem.click()
        time.sleep(1)
        driver.find_element_by_xpath(f"//*[@id='{pathway}-link']").click()
        time.sleep(1)
        driver.find_element_by_id(f"{pathway}-Table-link").click()
        time.sleep(3)
        driver.find_element_by_id(f"{pathway}-TableExport-link").click()
        time.sleep(4)
        
    driver.quit()  
    ## Reactome pathway만 저장이 안됌....
    for pathway in pathways:
        shutil.move(f"{download_dir}"+f"/{pathway}_table.txt",pathway_dir+f"/{pathway}_table.txt")
        # except FileNotFoundError:
        #     print(f"[Erorr] {pathway} table file was not downloaded...")
        #     print(f"[Error] you have to prepare this {pathway}table by your own hand! and then running the program again!")
            
    
def removeAllFile(dir):
    if os.path.exists(dir):
        for file in os.scandir(dir):
            os.remove(file.path)
    
def directory():
    root = Tk()
    root.dirName = filedialog.askdirectory() # directory 저장
    root.destroy()
    root.mainloop()
    
    
    return root.dirName

def making_pathway_lists(pathway):
    

    #pathway 받은 것을 dictionary로 바꿔서 각각의 key 값에 dataframe 형성
    p_name = [p[:p.find("_")] for p in pathway]
    Dic_pathway = {}
    for i in range(len(pathway)):
        Dic_pathway[p_name[i]] = pd.read_csv(f"{pathway[i]}.txt", delimiter = "\t")
        Dic_pathway[p_name[i]].insert(0,"Enrichr",f"{p_name[i]}")
        
    
    # data 전처리 이후 파일 만들
    ## NCI, Panther은 합치기 전 "Term"문자열 전처리 먼저 하기
    #뒤에 NCI, Pather regex ] "Homo sapiens" 찾고, 뒤에 문자열 모두 지우기.(단, 지우기 전에 열 복사 해서 original 남기기)
    #split으로 앞뒤 띄어쓰기 공백 없애기
    #"Overlap" m/d 형식으로 바꾸어주기
    
    #(Term열 복사 해서 original 남기기)
    # original 조건 :
    temp = ['NCI-Nature', 'Panther',"Reactome"]
    for key in Dic_pathway.keys():
        if key in temp:
            Dic_pathway[key].insert(1, "Original",Dic_pathway[key]["Term"])
        else:
            Dic_pathway[key].insert(1, "Original","")
            
            
    #raw data 파일 1차로 만들기
    df_raw = pd.concat(Dic_pathway.values(), axis =0)
    
    #합친 것 다시 인덱싱하기
    df_raw.reset_index(drop = True, inplace = True)
    
    #"Homo sapiens"를 가지고 있는 row에 대해....
    # 잡다한 것의 리스트
    junks = ["Homo sapiens"]
    
    for junk in junks:
        for key in Dic_pathway.keys():
            if key in temp:
                for i in range(len(Dic_pathway[key])):
                    if junk in Dic_pathway[key].loc[i,"Term"]:
                        idx = Dic_pathway[key].loc[i,"Term"].find(junk)
                        #뒤에 NCI, Pather regex ] "Homo sapiens" 찾고, 뒤에 문자열 모두 지우기.
                        #split으로 앞뒤 띄어쓰기 공백 없애기
                        Dic_pathway[key].loc[i,"Term"] = Dic_pathway[key].loc[i,"Term"][:idx].strip()
    
                    
        
    #지금까지 처리한 모든 dataframe합치기
    df_processed = pd.concat(Dic_pathway.values(), axis =0)
    #p-value 기준 오름차순 정렬하기
    df_processed.sort_values(by = ["P-value"],axis =0, ascending =True, inplace = True)
    #inplace: 변수 갱신 없이 자동으로 명령어으로 형성된 객체를 갱신
    #합친 것 다시 인덱싱하기
    df_processed.reset_index(drop = True, inplace = True)
    
    ##"Overlap" m/d 형식으로 바꾸어주기
    for idx in range(len(df_processed)):
        if str(type(df_processed.loc[idx,"Overlap"])) == "<class 'datetime.datetime'>":
            df_processed.loc[idx,"Overlap"]=df_processed.loc[idx,"Overlap"].strftime("%m/%d")
            # 버그 07-30 일 경우.... 7/30 이 아니라 7/3으로 됨...
            # 그냥 07/30 형식으로...
        else:
            continue
        
    ##"Term 같은 라인 중에 p-value 큰 row filtering(오름차순정리기 떄문에 첫번째 값이 가장 작은 값이므로 첫번째만 남기고 모두 drop)
    df_processed.drop_duplicates(["Term"], keep = "first",inplace = True)
    df_processed.reset_index(drop = True, inplace = True)
    
    
    #유의성 0.05 FILTERING을 할 경우
    df_final = df_processed.loc[df_processed["P-value"] < 0.05]
   
    #pathway filtering을 하지 않을 경우
    # df_final = df_processed

    return df_processed, df_final


def raw_UPDOWN_COUNT_file(keys, multiple):
    columns = ["Enrichr","Term","Overlap","P-value","Genes"]
    
    if len(keys) >1:
        df_dic = {}
        term_dic ={}
        
        for key in keys:
            df_dic[key] = pd.read_excel("Enrichr_raw.xlsx", sheet_name = f"{key}")[columns]
            term_dic[key] = df_dic[key]["Term"]
                
            
    
        #intersection 구하기
        for i in range(len(keys)):
            if i == 0 :
                term_intersec = set(term_dic[keys[i]])
                
            else:
                term_intersec = list(set(term_dic[keys[i]]).intersection(term_intersec))
                
        #교집합 term list의 index 구하기
        new_df_index = {}
        
        for key in keys:
            new_df_index[key] = []
            for idx in range(len(df_dic[key])):
                if df_dic[key].loc[idx,"Term"] in term_intersec:
                    new_df_index[key].append(idx)
                else: 
                    continue
                # for i in range(len(term_intersec)):
                #     if list(term_dic[key])[idx] == term_intersec[i]:
                #         new_df_index[key].append(idx)
        
        #index로 교집합 dataframe 구하기
        for key in keys:
            df_dic[key] = df_dic[key].loc[new_df_index[key]]
            
            df_dic[key].reset_index(drop = True, inplace = True)
            df_dic[key].sort_values(by = "Term", axis =0, ascending = True, inplace = True)
                
        #파일 저장
        with pd.ExcelWriter(f"{multiple}_UPDOWN_COUNT_raw.xlsx") as writer:
            for key in keys:
                df_dic[key].to_excel(writer, sheet_name = f"{key}", index = False)
   
    else:
        key = keys[0]
        df_dic = pd.read_excel("Enrichr_raw.xlsx", sheet_name = f"{key}")[columns]
        df_dic.to_excel(f"{multiple}_UPDOWN_COUNT_raw.xlsx", sheet_name = key, index =False)
    
#data2에서 logfc로 구분한 것들_ updown카운트를 위해 up, down구별해 주기
def UP_DOWN_LIST_file(keys):
    original = {}
    up = {}
    down = {}
    updown = {}
    columns = ["UP","updata","DOWN","downdata"]
    for key in keys:
        original[key] = pd.read_excel("UP_DOWN.xlsx", sheet_name = f"{key}")
        up[key] = original[key].loc[original[key].iloc[:,1]>=0]
        up[key].reset_index(drop = True, inplace = True)
        down[key] = original[key].loc[original[key].iloc[:,1]<0]
        down[key].reset_index(drop = True, inplace = True)
        
        updown[key] = pd.concat([up[key], down[key]], axis = 1)
        updown[key].columns = columns
        
    with pd.ExcelWriter("UP_DOWN.xlsx") as writer:
        for key in keys:
            original[key].to_excel(writer, sheet_name = f"ORIGINAL({key})", index = False)
            updown[key].to_excel(writer, sheet_name = f"{key}", index = False)
        
    

def counting_UPDOWN(up, down, data):
    gene_up_count = []
    gene_down_count = []
    num_up = []
    num_down = []
    
    for i in range(len(data)):
        gene_set = set(data.loc[i,"Genes"].split(";"))
        
        up_count= gene_set.intersection((up))
        down_count = gene_set.intersection((down))
        
        gene_up_count.append(up_count)
        gene_down_count.append(down_count)
        num_up.append(len(up_count))
        num_down.append(len(down_count))
        
        
        # print(f"{i}index_up : ",len(up_count))
        # print(f"{i}index_down :",len(down_count))
        
    return gene_up_count, gene_down_count, num_up, num_down
        


class UPDOWN:
    def __init__(self, keys, multiple):
        self.keys = keys
        self.multiple = multiple
        self.counted_dic = {}
        self.up = {} # up gene
        self.down = {}# down gene
        
        
        for key in keys:
            self.up[key] = set(pd.read_excel("UP_DOWN.xlsx", sheet_name = f"{key}")["UP"])
            self.down[key] = set(pd.read_excel("UP_DOWN.xlsx", sheet_name = f"{key}")["DOWN"])
            self.counted_dic[key]  = pd.read_excel(f"{multiple}_UPDOWN_COUNT_raw.xlsx", sheet_name = f"{key}")
           

    def UPDONW_DATAFRAME(self,gene_up_count, gene_down_count, num_up, num_down):
    
        
        for key in self.keys:
            #SET를 다시 한줄 쓰기로
            for i in range(len(gene_up_count[key])):
                
                gene_up_count[key][i] =list(gene_up_count[key][i])
                if len(gene_up_count[key][i]) >0: 
                    for idx in range(len(gene_up_count[key][i])):
                        if idx ==0:
                            temp = str(gene_up_count[key][i][idx])
                        else:
                            temp +=f";{str(gene_up_count[key][i][idx])}"
                        
                    
                    gene_up_count[key][i] =temp
                else: gene_up_count[key][i] = ""
                    
            for i in range(len(gene_down_count[key])):
                
                gene_down_count[key][i] =list(gene_down_count[key][i])
                if len(gene_down_count[key][i])>0:
                    for idx in range(len(gene_down_count[key][i])):
                        if idx ==0:
                            temp2 = str(gene_down_count[key][i][idx])
                        else: temp2+=f";{str(gene_down_count[key][i][idx])}"
                    
                    
                    gene_down_count[key][i] = temp2
                else: gene_down_count[key][i] = ""
                
            #DIC 재정의
            self.counted_dic[key]["UP"] = num_up[key]
            self.counted_dic[key]["DOWN"] = num_down[key]
            self.counted_dic[key]["UP_GENE"] = gene_up_count[key]
            self.counted_dic[key]["DOWN_GENE"] = gene_down_count[key]
            self.counted_dic[key].sort_values(by = ["Term"], axis = 0, ascending = True, inplace= True)        
        
        
        for i, key in enumerate(self.keys):
            if i ==0:
                counted_df = self.counted_dic[key]
            else: 
                counted_df = pd.concat([counted_df, self.counted_dic[key].iloc[:,2:]], axis = 1)
                
        return counted_df



print("[Setting] 작업할 폴더 위치를 지정해주세요[단, UP_DOWN.xlsx 파일도 같이 있어야합니다.]")
print("="*100+"\n")
main_dir = directory()


main_keys = openpyxl.load_workbook("UP_DOWN.xlsx").sheetnames
# main_keys = str(input("[Setting] Input your keys(ex. CON_102 CON_10G 102_10G): ")).split()

while(True):
    want_p_filtering = str(input("[Setting] Do you want p_value(0.05) filtering? : (yes or no): ")).lower()
    if want_p_filtering in ["yes", "no"]:
        break
    else:
        print("[Warning] Please answer with \"yes\", or \"no\"")
    
multiples = ["SINGLE", "DOUBLE", "TRIPLE", "QUADRA"]
main_keys_comb = {}
for i in range(len(main_keys)):
    main_keys_comb[multiples[i]] = [list(keys) for keys in list(combinations(main_keys,i+1))]



print("-"*100+"\n")
print(f"**Attention : 당신의 현재 key는 {main_keys}입니다. 기억해주세요!\n")
print("="*100+"\n")
#현재 진행된 상태 체크
print("[현재 준비된 파일- \"UP_DOWN.xlsx\" 상태 체크]")
print("조건 1 : 파일의 이름이 \"UP_DOWN.xlsx\"이다.")
print("조건 2 : 파일 내용은 DEG 실험 데이터를 logFC로 필터링한 Gene Symbol 리스트이다. ")
print("조건 3 : Gene Symbol 리스트의 첫번째 행은 gene symbol이 아닌 칼럼명이다.(칼럼명 내용은 무관)")

print("="*100+"\n")



answer = True
while(answer):
    
    INPUT_STATE = str(input("파일 준비 조건이 맞을 경우에만, Answer with \"yes\" : ")).lower()
    
    if INPUT_STATE == "yes":
        
        original_dic = {}
        os.chdir(main_dir)
        for key in main_keys:
            original_dic[key] = pd.read_excel("UP_DOWN.xlsx", sheet_name = f"{key}").iloc[:,0]
            #반드시 orginal sheet의 내용은 header가 있어야함. (내용 상관 없음.)
    
        UP_DOWN_LIST_file(main_keys)
        
        
        #case
        #case 1 : ENRICHR만 할 경우
        #CASE 2 : ENRICHR은 이미 했고, COMB_COUNT만 할 경우
        #CASE 3 : 한번에 둘다 할 경우
        #CASE 4 : 지금 안하고 프로그램을 종료하고 싶을 경우.
        print("")
        print("="*100+"\n")
        print("[What do you want me to do?]")
        print("-"*100+"\n")
        print("Mode : 1. ENRICHR(gene 리스트로 pathway raw 데이터 만들기")
        print("Mode : 2. COMB_COUNT(pathway raw 데이터로 double, triple... updown counting 하기\n")
        print("-"*100+"\n")
        print("Case 1 : Just do the 1. ENRICHR")
        print("Case 2 : I did 1.ENRICHR, I want you to do 2.COMB_COUNT")
        print("Case 3 : Do 1.ENRICHR & 2.COMB_COUNT at ONCE!")
        print("Case 4 : Let me think again, END this program\n")
        case = str(input("[YOUR CHOICE?] (1 : Case 1 | 2 : Case 2 | 3 : Case 3 | 4 : Case 4) = "))
        answer2 = True
        while(answer2):
            if case == "1":
                ENRICHR = True
                COMB_COUNT = False
                answer2 = False
                
            elif case == "2":
                ENRICHR = False
                COMB_COUNT = True
                answer2 = False
            
            elif case == "3":
                ENRICHR = True
                COMB_COUNT = True
                answer2 = False
                
            
            elif case == "4":
                ENRICHR = False
                COMB_COUNT = False
                answer2 = False
            
            else:
                ENRICHR = False
                COMB_COUNT = False
                answer2 = True
            
        
        answer = False
        break
        
        
    elif INPUT_STATE =="no":
        print("")
        print("[Warning] \"UP_DOWN.xlsx\"를 조건에 맞게 다시 준비해주세요!\n")
        print("[현재 준비된 파일- \"UP_DOWN.xlsx\" 상태 체크]")
        print("조건 1 : 파일의 이름이 \"UP_DOWN.xlsx\"이다.")
        print("조건 2 : 파일 내용은 DEG 실험 데이터를 logFC로 필터링한 Gene Symbol 리스트이다. ")
        print("조건 3 : Gene Symbol 리스트의 첫번째 행은 gene symbol이 아닌 칼럼명이다.(칼럼명 내용은 무관)")
        
        ENRICHR = False
        COMB_COUNT = False
        answer = False
        break
        
    else : 
        print("")
        print("[Warning] : Answer again!")
        ENRICHR = False
        COMB_COUNT = False
        answer = True
        
        
        
        
        
#ENRICHR : UP_DOWN 파일 GENE LIST를 가지고 pathway raw 데이터 만들기
if ENRICHR:

    print("")
    print("[PATHWAY EXAMPLES] \nBioPlanet_2019\nMSigDB_Hallmark_2020\nNCI-Nature_2016\nPanther_2016\nKEGG_2021_Human\nReactome_2016\n")
    all = ["BioPlanet_2019","MSigDB_Hallmark_2020", "NCI-Nature_2016","Panther_2016","KEGG_2021_Human","Reactome_2016" ]
    c = True
    while(c):
        if_all = str(input("모든 pathway에 대해 작업을 원하시면 \"all\"을, 원하는 pathway를 정하고 싶다면 \"no\"를 입력하세요 : ")).lower()
        
        if if_all == "all":
            pathways = all
            c = False
            
        elif if_all == "no":
            pathways = str(input("[Setting] Input your pathway file names(split with space) : \n")).split()
            c = False
            
        else:
            print("\n 답을 다시 입력해주세요! : \n")
            c = True
        
    pathways_table = [p+"_table" for p in pathways]
    
    print("\n[Setting] Chrome webdriver의 위치를 알려주세요!\n")
    webdriver_dir = directory() # 크롬 드라이버 위치 지정
    print("[Setting] 컴퓨터의 \"다운로드\" 위치를 알려주세요!\n")
    download_dir = directory()
    
    p_dir_dic ={}
    df_processed_dic = {}
    df_final_dic = {}
    #pathway web-crawling & pathway raw data 만들기
    for key in main_keys:
        # print(f"[Setting] <Before the processing, you should set \"pathways\" location of \"{key}\">\n")
        # print("[주의!] 선택하실 폴더는 아무 내용이 없는 새폴더를 만드셔야합니다.!\n")
        os.makedirs(main_dir+f"/{key}")
        p_dir_dic[key] = main_dir + f"/{key}"
        
        loop = True
        while(loop):
            try:
                print("[Loading] I AM DOWNLOADING THE PATHWYAS...\n")
                save_enrichr_txt(webdriver_dir,key,original_dic[key], pathways, p_dir_dic[key], download_dir)
                loop = False
                
            except SessionNotCreatedException:
                print("[Warning] please download right version of webdriver")
                
            except:
                print(f"[Error] During downloading \"{key}'s\" pathway table txt files, errors occurred...\n")
                print("[Error] But I'll do it again!\n")
                removeAllFile(p_dir_dic[key])
                loop = True
                
            
           
            # print("[Loading] I AM DOWNLOADING THE PATHWYAS...\n")
            # save_enrichr_txt(webdriver_dir,key,original_dic[key], pathways, p_dir_dic[key], download_dir)
            # loop = False
            
        
        
        print("[Loading] I AM MAKING THE ENRICHR PATHWAY RAW DATA...\n") 
        os.chdir(p_dir_dic[key])
        df_processed_dic[key], df_final_dic[key] = making_pathway_lists(pathways_table)
        print("="*100+"\n")
       
    #Enrichr_raw.xlsx 만들기
    print("[Loading] I AM SAVING \"Enrichr_raw.xlsx\"...\n")
    os.chdir(main_dir)
    
    
    
    
    with pd.ExcelWriter("Enrichr_raw.xlsx") as writer:
        if want_p_filtering == "yes":
            for key in main_keys:
                df_final_dic[key].to_excel(writer, sheet_name = f"{key}", index = False)
                
        elif want_p_filtering == "no":
            for key in main_keys:
                df_processed_dic[key].to_excel(writer, sheet_name = f"{key}", index = False)
        



#COMB_COUNT : DOUBLE, TRIPLE...COMBINATION FILE 만들기
# while(COMB_COUNT):
if COMB_COUNT:
    
    os.chdir(main_dir)
    
    # cmd_num = str(input("[Setting] Input combination number(1: SINGLE | 2 : DOUBLE | 3 : TRIPLE | 4 : Quadra | end : end program) > ")).lower()
    
    
    # if cmd_num == "1":
    #     # choice = str(input("Select one key from keys : "))
    #     # keys = [choice]
    #     keys_c = main_keys_comb["1"]
    #     multiple = "SINGLE"
    #     cmd = True
    #     COMB_COUNT = True
        
    # elif cmd_num == "2" :
    #     keys_c = main_keys_comb["2"]
    #     # # choice = str(input("Select two keys from keys : ")).split()
    #     # keys = choice
    #     multiple = "DOUBLE"
    #     cmd = True
    #     COMB_COUNT = True
        
    # elif cmd_num == "3":
        
    #     # keys = main_keys
    #     keys_c = main_keys_comb["3"]
    #     multiple = "TRIPLE"
    #     cmd = True
    #     COMB_COUNT = True
        
    # elif cmd_num =="4":
    #     keys_c = main_keys_comb["4"]
    #     multiple = "Quadra"
    #     cmd = True
    #     COMB_COUNT = True
        
    
    # elif cmd_num =="end":
    #     COMB_COUNT = False
    #     cmd = False
    #     break
    
    # else:
    #     print("[Warning] your input was wrong, do it again!")
    #     print("="*100+"\n")
    #     cmd = False
    #     COMB_COUNT = True
        
    
    for multiple, keys_c in main_keys_comb.items():
    # if cmd:
        
        print(f"\nyour current keys : {main_keys}\n")
        print(f"[LOADING] I am making {multiple} combinations of keys...\n")
        # multiple = str(input("Multiple(NONE OR DOUBLE OR TRIPLE) :")).upper()
        
        for keys in keys_c:
            
            comb_name = ",".join(keys)
            comb_name = f"({comb_name})"
            
            raw_UPDOWN_COUNT_file(keys, multiple)
            
            
            UD = UPDOWN(keys, multiple)
            gene_up_count ={} #up에 해당되는 gene들
            gene_down_count = {} # down에 해당되는 gene들 
            num_up = {} #upcount 개수
            num_down = {} # downcount 개수
                    
            for key in keys:
                gene_up_count[key], gene_down_count[key], num_up[key], num_down[key]= counting_UPDOWN(UD.up[key], UD.down[key], UD.counted_dic[key])
                
            
            counted_df = UD.UPDONW_DATAFRAME(gene_up_count, gene_down_count, num_up, num_down)
            
            with pd.ExcelWriter(f"{UD.multiple}_UPDOWN_COUNT_{comb_name}.xlsx") as writer:
                counted_df.to_excel(writer, sheet_name = f"{UD.multiple}_UPDOWN", index = False)
                for key in keys:
                    UD.counted_dic[key].to_excel(writer, sheet_name = f"{key}", index = False)
                #sheet에 덮어쓰기
            cmd = False
            
            
            
        

