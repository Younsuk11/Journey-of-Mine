# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 17:06:44 2021

@author: curig
"""

import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
from tqdm import tqdm
from collections import defaultdict
import os
from tkinter import filedialog
from tkinter import *
import datetime
from scipy import linalg
from scipy import stats
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from skbio.alignment import StripedSmithWaterman
from sklearn.metrics import confusion_matrix,accuracy_score, roc_auc_score, roc_curve
import statsmodels.api as sm


def directory():
    root = Tk()
    root.dirName = filedialog.askdirectory() # directory 저장
    root.destroy()
    root.mainloop()
    
    
    return root.dirName

class sequence_scoring:
    def __init__(self, seqA, seqB, shift):
        self.seqA = seqA
        self.seqB = seqB
        self.seq_with_s = ""
        self.seq_with_x = ""
        self.seqB_with_s = "" #3'~5'
        self.shift = shift
        
        

    def mismatch(self): #mismatch가 표시된(소문자), seqA, seqB,출력하고, 각각의 index정보 출력하기
        
        rev_seqB = list(str(Seq(self.seqB).reverse_complement()))
        l_seqA = list(self.seqA)
        seqA_with_s = list(self.seqA)
        seqA_with_x = list(self.seqA)
        seqB_with_s = list(str(Seq(self.seqB).reverse_complement().complement()))
        m_pos_A = []
        
        
        if self.shift >0:
            rev_seqB = rev_seqB[self.shift:]
            for i in range(len(rev_seqB)):
                if l_seqA[i] != rev_seqB[i]:
                    m_pos_A.append(i)
                    
        elif self.shift <0:
            rev_seqB = rev_seqB[:-self.shift]
            for i in range(len(rev_seqB)):
                if l_seqA[i-self.shift] != rev_seqB[i]:
                    m_pos_A.append(i-self.shift)
            
        else:
            for i in range(len(rev_seqB)):
                if l_seqA[i] != rev_seqB[i]:
                    m_pos_A.append(i)
                    
        
        for idx in m_pos_A:
            seqA_with_s[idx] = seqA_with_s[idx].lower()
            seqA_with_x[idx] = "x"
            seqB_with_s[idx] = seqB_with_s[idx].lower()
            
        
        self.seq_with_x = "".join(seqA_with_x)
        self.seq_with_s = "".join(seqA_with_s)
        self.seqB_with_s = "".join(seqB_with_s)
        
        

    def GC_ratio(self, start, end): #mismatch 고려, 고려X, - one strand 기준으로 출력.
        GC_content = {}
        GC_content["no_mismatch"] = GC(Seq(self.seqA[start:end]))
        GC_content["yes_mismatch"] = GC(Seq(self.seq_with_x[start:end]))
        
        return GC_content

    def AU_count(self, start, end):#mismatch 고려, 고려X, - one strand 기준으로 출력
        AU_content = defaultdict(int)
        t1 = self.seqA[start:end].replace("U", "A")
        t2 = self.seq_with_x[start:end].replace("U", "A")
        if t1.count("A") >= 3:
            AU_content["no_mismatch"] = 1
        else:AU_content["no_mismatch"] = 0
        
        if t2.count("A") >= 3:
            AU_content["yes_mismatch"] = 1
        else:AU_content["yes_mismatch"] = 0
                
        return AU_content
    
    def first_pos(self):
        
        if self.seqA[0] in ["A","U"]:

            return 1
        else: return 0
        
    def last_pos(self):
        if self.seqA[-1] in ["G", "C"]:
            return 1
        else: return 0
    
    def A_10(self):
        l = len(self.seqA)
        if l % 2 == 0:
            m_a = int((l/2) -1)
            m_b = int(l/2)
            
            if self.seqA[m_a] == "A" and self.seqA[m_b] == "A":
                return 1
            else: return 0
        else:
            m = int((l-1)/2)
            if self.seqA[m] == "A":
                return 1
            else: return 0
            
            
    
    def long_GC(self):
        
        temp = self.seqA.replace("U", "x")
        temp = temp.replace("A", "x")
        m = max(list(map(len, temp.split("x"))))
        if m >9 : 
            return 0
        else: return 1
    

    def Tm(self, start, end):#mismatch 고려, 고려X, hairpin structure(internal sequence)- one strand 기준으로 출력
        Tm_values = {}
        
        # shift고려 안함.
        Tm_values["no_mismatch"]  = mt.Tm_Wallace(self.seqA[start:end])
        Tm_values["yes_mismatch"] = mt.Tm_NN(self.seqA[start:end], c_seq = self.seqB_with_s[start:end], shift = self.shift, strict = False, nn_table = mt.RNA_NN3, de_table = mt.RNA_DE1)
        
        return Tm_values
    
    
def H_loop(seq):
    dic = defaultdict(list)
    internal_tm = []
    for i in range(len(seq)):
        test1 = sequence_scoring(seq[i:], seq[i:], 0) # shift<0
        test2 = sequence_scoring(seq[:len(seq)-i], seq[:len(seq)-i], 0) # shift>0
        test1.mismatch()
        test2.mismatch()
        
        for idx, t in enumerate([test1, test2]):

            l =len(t.seq_with_x)

            
            loop = 0
            

            if l %2 == 0: #짝수
                m = int(l/2)
            else: #홀수
                m = int((l -1)/2)

                
            # mismatch x로 이루어진 loop base counting    
            for n, c in enumerate(t.seq_with_x[m:]):
                if c == "x":
                    loop += 1
                else:
                    start = m+n
                    break
            for n, c in enumerate(reversed(t.seq_with_x[:m])):
                if c == "x":
                    loop +=1
                else: break

                    
            #연속적인 stem base세기
            stem = max(list(map(len, t.seq_with_x[m:].split("x"))))
                    
                    
            #조건 filtering
            if 4<= loop : 
                if stem >=4:
                    if idx == 0:
                        dic[-i].append(t.seq_with_x)
                        dic[-i].append(t.seq_with_s)
                        dic[-i].append(t.Tm(start, l)["yes_mismatch"])
                        internal_tm.append(t.Tm(start, l)["yes_mismatch"])
                        
                        
                    elif idx == 1:
                        if i!= 0:
                            dic[i].append(t.seq_with_x)
                            dic[i].append(t.seq_with_s)
                            dic[i].append(t.Tm(start, l)["yes_mismatch"])
                            internal_tm.append(t.Tm(start, l)["yes_mismatch"])
    if len(internal_tm) >0:                        
        m = max(internal_tm)
        if m >20:
            return 0
        else: return 1
    else: return 1
    

def aligned_query_sequence(Revcomp_seq, target_seq):
    # striped SW symbol 생성 및 SW score 계산
    Revcomp_seq = StripedSmithWaterman(Revcomp_seq) #유전자
    Result = Revcomp_seq(str(target_seq)) #sequence

    A_gene_seq = Result['aligned_query_sequence']
    B_target_seq = Result['aligned_target_sequence']


    Smt_score = int(Result['optimal_alignment_score'])
    # Smt_score = int(Result['suboptimal_alignment_score'])


    return Smt_score

def scoring_siRNA_thermodynamic(Df, criteria, cmd,n):
    
    df = Df.copy()
    
    
    cmd = cmd.lower() # cmd별로 yes vs no mismatch sheet 나누기s
    
    
    for i in tqdm(range(len(df))):
        score = sequence_scoring(df.loc[i,"A_siRNA"], df.loc[i,"B_siRNA"], shift = 0)
        score.mismatch()
        l = len(score.seqA)
        seed = int(round(l/3))
        
        df.loc[i, "first_A/U"] = score.first_pos()
        df.loc[i, "last_G/C"] = score.last_pos()
        df.loc[i, "U_10"] = score.A_10()
        df.loc[i, "GC_stretch"] = score.long_GC()
        df.loc[i, "Tm_hloop"] = H_loop(score.seqA)
        
        
        
        
        
        if cmd == "yes":
            #whole GC
            gc = score.GC_ratio(0, l)["yes_mismatch"]
            if gc >= criteria["GC_content_under"] and gc <= criteria["GC_content_upper"]:
                df.loc[i, "GC_content"] = 1
            else: df.loc[i, "GC_content"] = 0
            
            
            df.loc[i, "seed_3_A/U"] = score.AU_count(1, seed)["yes_mismatch"]
            
            # GC, Tm값들 조건에 따른 scoring하기    
            #seed TM <20
            if score.Tm(1, seed)["yes_mismatch"] <= criteria["seed_Tm"]:
                df.loc[i, "seed_Tm"] = 1
            else:  df.loc[i, "seed_Tm"] = 0
            
            
            #seed GC
            if int(round(score.GC_ratio(1,seed)["yes_mismatch"])) == criteria["GC_seed"]:
                df.loc[i, "GC_seed"] = 1
            else: df.loc[i, "GC_seed"] =0
            
            
            #non seed GC
            if int(round(score.GC_ratio(seed,l)["yes_mismatch"])) == criteria["GC_non_seed"]:
                df.loc[i, "GC_non_seed"] = 1
            else: df.loc[i, "GC_non_seed"] =0
            
            
        elif cmd == "no":
        
            #whole GC
            gc = score.GC_ratio(0, l)["no_mismatch"]
            if gc >= criteria["GC_content_under"] and gc <= criteria["GC_content_upper"]:
                df.loc[i, "GC_content"] = 1
            else: df.loc[i, "GC_content"] = 0
        
        
        
            df.loc[i, "seed_A/U"] = score.AU_count(1, seed)["no_mismatch"]
        
        
            # GC, Tm값들 조건에 따른 scoring하기    
                #seed TM <20
        
            if score.Tm(1, seed)["no_mismatch"] <criteria["seed_Tm"]:
                 df.loc[i, "seed_Tm"] = 1
            else:  df.loc[i, "seed_Tm"] = 0    
            
    
            #seed GC
            
            if int(round(score.GC_ratio(1, seed)["no_mismatch"]))==criteria["GC_seed"]:
                df.loc[i, "GC_seed"] = 1
            else: df.loc[i, "GC_seed"] =0   
        
    
            #non seed GC
        
            if int(round(score.GC_ratio(seed,l)["no_mismatch"])) ==criteria["GC_non_seed"]:
                df.loc[i, "GC_non_seed"] = 1
            else: df.loc[i, "GC_non_seed"] =0    
            
  
        
    df["score"] = df.apply(lambda x: np.sum(x[n:]), axis = 1)  ## mismatch를 고려하지 않은 총점
    
    cols = list(df.columns)[n:-1]
    
    for col in cols:
        a = np.sum(df[col]) / len(df[col])
        if a >= 0.9 or a <= 0.1:
            df.drop(col, axis =1, inplace =True)
            
    cols = list(df.columns)[n:-1]
    
    return df, cols


def reynold(Df, criteria, cmd, cmd2,n):
    
    
    df = Df.copy()
    
    
    cmd = cmd.lower() # cmd별로 yes vs no mismatch sheet 나누기s
    cmd2 = cmd2.lower()
    
    for i in tqdm(range(len(df))):
        score = sequence_scoring(df.loc[i,"A_siRNA"], df.loc[i,"B_siRNA"], shift = 0)
        score.mismatch()
        
        
        sense = str(Seq(score.seqA).reverse_complement())
        l = len(score.seqA)
        seed = int(round(l/3))
        
        ratio = float(l/19)
        
        df.loc[i, "Tm_hloop"] = H_loop(score.seqA)
        
        if cmd == "yes":
            gc = score.GC_ratio(0, l)["yes_mismatch"]
            if gc >=criteria["GC_content_under"] and gc <= criteria["GC_content_upper"]:
                df.loc[i, "GC_content"] = 1
            else: df.loc[i, "GC_content"] = 0
            
            #3' of sense  = seed of antisense
            df.loc[i, "seed_3_A/U"] = score.AU_count(1, seed)["yes_mismatch"]
            
            
        if cmd == "no":
            gc = score.GC_ratio(0, l)["no_mismatch"]
            if gc >=criteria["GC_content_under"] and gc <= criteria["GC_content_upper"]:
                df.loc[i, "GC_content"] = 1
            else: df.loc[i, "GC_content"] = 0
            
            
            df.loc[i, "seed_3_A/U"] = score.AU_count(1, seed)["no_mismatch"]
            
        
        if sense[l-1] == "A":
            df.loc[i, "A_19"] = 1
        else:
            df.loc[i, "A_19"] = 0
        
        if cmd2 == "ratio":

            if sense[round(2*ratio)] == "A":
                df.loc[i, "A_3"] = 1
            else: df.loc[i, "A_3"] = 0
            
            if sense[round(9*ratio)] == "U":
                df.loc[i, "U_10"] = 1
            else: df.loc[i, "U_10"] = 0
            
            if sense[round(12*ratio)] == "G":
                df.loc[i, "not_G_13"] =0
            else: df.loc[i, "not_G_13"] = 1
            
        elif cmd2 == "fixed":
            if sense[2] == "A":
                df.loc[i, "A_3"] = 1
            else: df.loc[i, "A_3"] = 0
            
            if l %2 == 0: # seqA 길이가 짝수
                
                if sense[int((l/2)-1)] == "U" and sense[round(l/2)] == "U":
                    df.loc[i, "U_10"] = 1
                else: df.loc[i, "U_10"] = 0
                
            
            
                if sense[int((l/2)-1)+3] == "G":
                    df.loc[i, "not_G_13"] =0
                else: df.loc[i, "not_G_13"] = 1
                
            else: #길이가 홀수
                if sense[int((l-1)/2)] == "U":
                    df.loc[i, "U_10"] = 1
                else: df.loc[i, "U_10"] = 0
                
                if sense[int((l-1)/2)+3] == "G":
                    df.loc[i, "not_G_13"] =0
                else: df.loc[i, "not_G_13"] = 1
        
        
        if sense[l-1] in ["G", "C"]:
            df.loc[i, "not_GC_19"] = 0
        else: df.loc[i, "not_GC_19"] = 1
        
    
    df["score"] = df.apply(lambda x : np.sum(x[n:]), axis = 1)
        
        
    cols = list(df.columns)[n:-1]     
    
    for col in cols:
        a = np.sum(df[col]) / len(df[col])
        if a >= 0.9 or a <= 0.1:
            df.drop(col, axis =1, inplace =True)
                    
    cols = list(df.columns)[n:-1]  
            
    return df,cols



def efficacy_classification(x, criteria):
    if x > criteria:
        return 0
    
    else:
        return 1
    
class Calculator:
    
    def __init__(self, df, criteria, efficacy_col):
        self.df = df
        
        self.value = np.array(df[efficacy_col])
        
        self.efficacy = np.array(self.df[efficacy_col].map(lambda x : efficacy_classification(x, criteria)))
        
        self.one_value = 1-self.value
        self.one_value_mean = self.one_value.mean()
        
        
        
    # TP, TN, FP, FN classification
    def basic(self, col):
        
        series = np.array(self.df[col])
        cfmat = confusion_matrix(self.efficacy, series)
        
        
        
        precision = cfmat[0,0] / (cfmat[0,0] + cfmat[0,1])
        
        recall = cfmat[0,0] / (cfmat[0,0] + cfmat[1,0])
        
        specificity = cfmat[1,1] / (cfmat[1,1] + cfmat[0,1])
            
        accuracy = (cfmat[0,0] + cfmat[1,1])/ (cfmat[0,0]+cfmat[0,1]+cfmat[1,0]+cfmat[1,1])
            
        
        
        
        return precision, recall, specificity, accuracy
    
    def f_score(self, precision, recall):
        if precision ==0 or recall == 0:
            return 0
        else:
            return (2*precision*recall) / (precision + recall)
        


def standardization(array):
    std = np.sum((array-array.mean())**2)/len(array)
    root = len(array)**(1/2)
    return (array-array.mean())/(std /root)

def scaler_normal(array):
    return (array-min(array)) / (max(array) - min(array))
    

def opp_tag(x):
    temp = x.split("_")
    temp = [temp[0], temp[2], temp[1], temp[3]]
    
    return "_".join(temp)


def linear_regression(x, y,k, special):
    coef = []
    R = []
    degree =[]
    Corr_p = []
    special_index = special.index
    
    
    print("")
    print("*"*100)
    a = 0
    while(True):
        
        while(True):
            print("")
                            
            print("[Setting] If you want to continue the fitting process, then input Linear fitting degree")
            print("[Setting] If you want to stop fitting, then input \"stop\" command")
                            
            n = str(input("[Setting] Input Linear Regression polynomial fitting degree (ex_ 1, 2, 3, 4... or stop) : ")).lower()
                           
                            
            if n == "stop":
                temp = False
                temp2 = True
                    
                break
            else: 
                try:
                    d = int(n)
                    temp = True
                    temp2 = False
                    break
                except:
                    print("[Warning] Please input Your integer degree again\n ")
                        
        if temp:
                    
            
            #all data model
            model = LinearRegression()
            model.fit(np.vander(x,d+1), y)
            p_y = model.predict(np.vander(x,d+1))
            coef_ = model.coef_[:-1]
            coef_ = np.append(coef_, model.intercept_)
            corr_p = stats.pearsonr(x,y)
            ssr = np.sum((y-p_y)**2)
            sst = np.sum((y-y.mean())**2)
            r=1-(ssr/sst)
                            
            x_many = np.arange(x.min(), x.max(), 0.1)
            y_many = np.vander(x_many, d+1)@coef_
                          
            
            
            print("[Result] Your fitting coefficient : {}".format(coef_))
            print("[Result] Determinant of Coefficient(Manual) : {}\n".format(r))    
            print("[Result] Pearson Coeffiecient : {}\n".format(corr_p))
                                
                            
            plt.suptitle("Scoring Criteria : {}".format(k), size = 10)
            plt.plot( x, y, "bo", label = "all seq")
            
            plt.plot(x[special_index], y[special_index],"yo", label = "both knockdown")
            for idx in special_index:
                plt.text(x[special_index][idx], y[special_index][idx]+0.02, special["Group"][idx], ha = "center", color = "black", fontsize = "small")
            plt.plot(x_many, y_many, "r-", label = "R^2 : {}\ncoefficient : {}".format(round(r,7), round(corr_p[0],7)))
            
            plt.xlabel("Predicted Score")
            plt.ylabel("Efficiency(1-VALUE)")
            plt.legend(loc  = "upper left")
            plt.savefig("Linear_{}_{}{}{}.png".format(k, Date.year, Date.month, Date.day), dpi = 100)
            plt.show()
                            
            a = 1
                    
                            
        if temp2:
            if a == 1:
                R.append(r)
                Corr_p.append(corr_p)
                coef.append(coef_)
                degree.append(d)
                
                break
            else:
                print("[Warning] You just skip {} scoring method\n".format(k))
                        
                    
                            
    coef.insert(0,"-")
    coef.insert(1,"-")
    degree.insert(0,"-")
    degree.insert(1,"-")
    R.insert(0,"-")
    R.insert(1,"-")
    Corr_p.insert(0,"-")
    Corr_p.insert(1,"-")
    
    
    return coef, R, degree, Corr_p, model

def cut_off(y, threshold):
    Y = y.copy()
    Y[Y>threshold] = 1
    Y[Y<threshold] = 0
    
    return Y.astype(int)

def acc(cfmat):
    return (cfmat[0,0] + cfmat[1,1])/ (cfmat[0,0]+cfmat[0,1]+cfmat[1,0]+cfmat[1,1])

def logistic_regression(x,y,k):
    x = sm.add_constant(x, has_constant = "add")
    model = sm.Logit(y,x)
    results =model.fit(method = "newton")
    params = results.params
    p_y = results.predict(x)
    
    
    
    
    threshold = np.arange(0,1, 0.1) #0~1까지 0.1씩 증가
    table = pd.DataFrame(columns = ["ACC"])
    for i in threshold:
        pred_Y = cut_off(p_y, i)
        cfmat = confusion_matrix(y, pred_Y)
        
        if acc(cfmat):
            table.loc[i] = acc(cfmat)
        else:
            table.loc[i] = 0
        
    table.index.name = "threshold"
    table.columns.name = "performance"
    
    print(table)
    
    
    p_Y = cut_off(p_y, table.loc[table["ACC"] == max(table["ACC"])].index.max())
    cfmat = confusion_matrix(y, p_Y)
    accuracy = acc(cfmat)
    
    
    print("*"*100)
    print("[Logistic Regression]")
    print(results.summary())
    print("")
    print("[Result] Parmaeter of regression")
    print(params)
    print("")
    print("[Result] Accuracy : {}".format(accuracy))
    
    #roc graph
    
    fpr, tpr, thresholds = roc_curve(y, p_y, pos_label = 1)
    auc = np.trapz(tpr, fpr)
    plt.suptitle("Scoring Criteria : {}".format(k), size = 10)
    plt.plot(fpr, tpr, label = "Accuracy : {}\nAUC : {}".format(accuracy, auc))
    plt.xlabel("False Positive Rate")
    plt.ylabel("Treu Positive Rate")
    plt.legend()
    plt.show()
    plt.savefig("Logistic_{}_{}{}{}.png".format(k, Date.year, Date.month, Date.day), dpi = 100)
    
    
    print("")
    print("[Result] AUC value : {}".format(auc))
    
    return p_y, p_Y
    
    


Date =datetime.datetime.now()
pathway = directory()
os.chdir(pathway)

cols = {}
scored_df = {}

print("\"Welcome\"\n")

          
start = True  
while(start):

    
    cmd0 = str(input("\"Start\" or \" Stop\" : ")).lower()
    
    #button1
    if cmd0 == "start":
        
        # 0. <Setting>
        print("[Setting] \"efficacy\" criteria")
        print("\t(example) \"efficacy=1\" : there is no regulation change")
        print("\t(example) \"efficacy=0.5\" : 50% down regulation")
        print("\t(example) \"efficacy=1.5\" : 50% up regulation")
      
        filtering_criteria = input("[Setting] \"efficacy\" criteria for fitering(no need to filter, then answer with \"no\") : ")
            
        
        efficacy_criteria = float(input("[Setting]\"good efficacy\" criteria : "))
        print("")
                                        
                
        criteria = defaultdict(int)
        #36-54-19-54-20
        print("[Setting] Article default values are : \"36-54-19-54-20\"")
        criteria["GC_content_under"] = int(input("[Setting] whole seq GC_content's under limit : "))
        criteria["GC_content_upper"] = int(input("[Setting] whole seq GC_content's upper limit : "))
        criteria["GC_seed"] = int(input("[Setting] seed seq GC_content limit : "))
        criteria["GC_non_seed"] = int(input("[Setting] non_seed seq GC_content limit : "))
        criteria["seed_Tm"] = int(input("[Setting] seed seq Tm limit : "))
        
        
        while(True):
            sure = input("[Setting] Are you sure about your setting values? (yes or no)\n").lower()
            if sure == "yes":
                
                break
            elif sure == "no":
                break
            else:
                continue
            
        
            
        if sure == "no":
            break
            #botton3 
            #점수 채점 결과
        
        
        #1. Input
        while(True):
            try:
                file = str(input("[Setting] Input your data file name!(format : xlsx): "))
                Df = pd.read_excel("{}.xlsx".format(file))
                
                #<filtering>
                try: 
                    filtering_criteria = float(filtering_criteria)
                    Df = Df.loc[Df["Value"] <= filtering_criteria].reset_index(drop =True)
                except:
                    print("[Setting] You answered with \"no\" for \"efficacy criteria\"\n")
                
                
                Df["Smt_score"] = Df.apply(lambda x : aligned_query_sequence(str(Seq(x["A_sequence"]).reverse_complement()), x["B_sequence"]) , axis = 1)
                basis = list(Df.columns)
                basis.append("score")
                n = len(Df.columns) #input data 기존 칼럼 수
                special = Df.loc[Df["dot"] == "O"]
                
                break
            #warning
            except:
                print("Wrong Directory\n")
                start = False
                break
                
        if not start:
            break
        
        
    
        #0. <Setting>
        while(True):
            print("-"*100)
            print("[Setting] Select your mode \n")
            print("\"yes\" : yes! consider sequence mismatches")
            print("\"no\" : no! consider sequence mismatches")
          
            cmd1 = str(input("[Setting] Select your mode_ (yes_mismatch | no_mismatch ) :   ")).lower()
            
            
            if cmd1 == "yes" or cmd1 == "no":
                scored_df["THdynamic"], cols["THdynamic"] = scoring_siRNA_thermodynamic(Df, criteria, cmd1,n)
                scored_df["reynold_ratio"],  cols["reynold_ratio"] = reynold(Df, criteria, cmd1, "ratio", n)
                scored_df["reynold_fixed"], cols["reynold_fixed"] = reynold(Df, criteria, cmd1,"fixed",n)
                
                scored_df["Comb_scoring(fixed)"], cols["Comb_scoring(fixed)"]= reynold(scored_df["THdynamic"].iloc[:,:-1], criteria, cmd1, "fixed", n)
                scored_df["Comb_scoring(ratio)"],cols["Comb_scoring(ratio)"] = reynold(scored_df["THdynamic"].iloc[:,:-1], criteria, cmd1, "ratio", n)
                
                break
            
            else: 
                print("[Warning] Your command was wrong, please selet again!\n")
        
        
        #2. Statistic Machine learning Weights for binary to continuous data form  
        weight = {} 
        for k in scored_df.keys():
            
            print("")
            print("[Current] current scoring criteria is \"{}\"\n".format(k))
            
            
            #초기 설정
            cal = Calculator(scored_df[k], efficacy_criteria, "Value") # criteria 1 : 1보다 작은게 efficiency
            columns = cols[k].copy()
            weight[k] = pd.DataFrame(columns = columns)
            
            
            #fscore , Accuracy, 
            f_weights = []
            A_weights = []
            # LR_weights = []
            for col in columns:
                precision, recall, specificity, accuracy = cal.basic(col)
                f_weights.append(cal.f_score(precision, recall))
                A_weights.append(accuracy)
               
            
            weight[k].loc["f_score"] = f_weights
            weight[k].loc["Accuracy"] = A_weights
            
            
            
            #least square weight(linear coefficient)
            # y = cal.one_value # continuous
            y = cal.efficiency
            df1 = sm.add_constant(scored_df[k][columns], has_constant = "add")
            
            #solution linear regression 모델 학습
            multi_model = sm.OLS(y, df1)
            fitted_multi_model = multi_model.fit()
            params = fitted_multi_model.params
            
            
            
            weight[k].insert(0,"const", "-")
            weight[k].loc["lstsq"] = params #beta weight
            
            #broad casting
            scored_df[k]["lstsq_score"] = pd.Series(np.sum(df1*params,axis = 1))
            
            
            
            #group mapping(Dualing)
            # scored_df[k]["tag"] = scored_df["THdynamic"].apply(lambda x : str(x["Group"]) +"_"+ str(x["pmol"]), axis = 1)
            # scored_df[k]["opp_tag"] = scored_df["THdynamic"]["tag"].map(lambda x : opp_tag(x))
            #점수 합산
            # for i in tqdm(range(len(scored_df[k]))):
            #     idx = scored_df[k].loc[scored_df[k]["tag"] == scored_df[k].loc[i,"opp_tag"]].index[0]
            #     scored_df[k].loc[i,"Total_sum"] = float(scored_df[k].loc[i,"score"]) + float(scored_df[k].loc[idx,"score"])
            #     scored_df[k].loc[i,"Total_lstsq"] = float(scored_df[k].loc[i,"lstsq_score"]) + float(scored_df[k].loc[idx,"lstsq_score"])
            
            
            
            
            
            
            
            
            
            
            #3. ML analysis (regression(linear, logistic), classification) -- by class Data_analysis
            
            #3-2. Regression
            #3-2-1. Linear Regression
            coef, R, degree, Corr_p, model= linear_regression(scaler_normal(scored_df[k]["lstsq_score"]), cal.one_value,k,special)
            try:
                weight[k]["coef[{}~{}]".format("High", "Low")] = coef
                weight[k]["degree "] =degree
                weight[k]["R^2_manual"] = R
                weight[k]["pearson coeff "] = Corr_p
            
            except:
                print("[Setting] You skip \"{}\" scoring criteria\n".format(k))
            
        
            #3-2-1. Logistic Regression - 단순 로지스틱 회귀
            p_y, p_Y = logistic_regression(standardization(scored_df[k]["lstsq_score"]), cal.efficacy, k)
            #Final inference
            scored_df[k]["classification"] = p_Y
            scored_df[k]["probability"] = p_y
        
                    
        
        #저장
        with pd.ExcelWriter("R_{}_mis({}_{}_{})_e{}.xlsx".format(cmd1, Date.year, Date.month, Date.day,efficacy_criteria)) as writer:
            for k in scored_df.keys():
                # scored_df[k].drop(cols[k], axis = 1, inplace = True)
                # scored_df[k].drop(["tag", "opp_tag"], axis = 1, inplace = True)
                scored_df[k].to_excel(writer, sheet_name = k, index = False)
                weight[k].to_excel(writer, sheet_name = k+"_weights")
                
        
            
        
        
        
        
        
        
        
    
        start = True
        
    elif cmd0 == "stop":
        print("\"End this program\"")
        start = False
        break

    else:
        print("\"Your command is wrong, please do it agian!\"")
        start = True

#%%
    

#21.12.29
# 실험값 없이 불가능.. clustering 
# 결론 : 결과가 있어야 예측.
#아니면 진짜 데이터베이스 뒤져야한다.
#근데 siRNA 발현 데이터베이스 잇나?
# dual targeting  문제는 또 어떻게 해결?  - random shuffling안되는데.... 상보성을 유지하면서 random? 
#
#각각 PSSM하고 (random score - 점수 계산(sum) - pvalue threshold- frequency profile - 실험 - 가중치 갱신 - 높은거 선별 - mismatch 상보)

#상보성 문제는 맨 마지막?


#21.1.05

#듀얼 모두 떨어진다(criteria 기준)
#모두 떨어지는 것 figure에 표시하기

#21.1.11
#모든 칼럼 한번에 regression하는게 나을듯
#input 데이터 전처리 상태 확인 중복들... 실험값 unique성 
#regression 전제 확인할 필요도 있음.









#%% smithwaterman
from skbio import TabularMSA, DNA
from skbio.alignment import local_pairwise_align_ssw
alignment, score, start_end_positions = local_pairwise_align_ssw(DNA("ACTAAGGCTCTCTACCCCTCTCAGAGA"),
                                                                  DNA("ACTAAGGCTCCTAACCCCCTTTTCTCAGA"))
#%%
alignment
TabularMSA[DNA]
------------------------------
Stats:
    sequence count: 2
    position count: 30
------------------------------
ACTAAGGCTCTC-TACCC----CTCTCAGA
ACTAAGGCTC-CTAACCCCCTTTTCTCAGA
>>> score
27
>>> start_end_positions
[(0, 24), (0, 28)]

#%%
alignment, score, start_end_postions = local_pairwise_align_ssw(DNA(str(Seq("GGCAaGCTGGGCgTcTTGC").reverse_complement()).upper()),
                                                                DNA("GCAAcAgGCCCAGCaTGCC".upper()))

#%%
from skbio.alignment import StripedSmithWaterman
query = StripedSmithWaterman("ACTAAGGCTCTCTACCCCTCTCAGAGA")
alignment = query("AAAAAACTCTCTAAACTCACTAAGGCTCTCTACCCCTCTTCAGAGAAGTCGA")
