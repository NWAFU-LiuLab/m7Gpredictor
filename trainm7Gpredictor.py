# -*- coding: utf-8 -*-
"""
Created on Wed Jun 05 21:59:53 2019

@author: Administrator
"""

from functools import reduce
import operator
import numpy as np
from sklearn import preprocessing
from sklearn.model_selection import cross_validate
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.ensemble import AdaBoostClassifier
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.svm import SVC
import subprocess
from repDNA.psenac import PseDNC
from repDNA.psenac import PseKNC
def kmer(seq):
    mer2={}
    mer3={}
    mer4={}
    for n1 in 'ATCG':
        for n2 in 'ATCG':
            mer2[n1+n2]=0
            for n3 in 'ATCG':
                mer3[n1+n2+n3]=0
                for n4 in 'ATCG':
                    mer4[n1+n2+n3+n4]=0
    seq=seq.replace('N','')
    seq_len=len(seq)
    for p in range(0,seq_len-3):
        mer2[seq[p:p+2]]+=1
        mer3[seq[p:p+3]]+=1
        mer4[seq[p:p+4]]+=1
    mer2[seq[p+1:p+3]]+=1
    mer2[seq[p+2:p+4]]+=1
    mer3[seq[p+1:p+4]]+=1
    v2=[]
    v3=[]
    v4=[]
    for n1 in 'ACGT':
        for n2 in 'ACGT':
            v2.append(mer2[n1+n2])
            for n3 in 'ACGT':
                v3.append(mer3[n1+n2+n3])
                for n4 in 'ACGT':
                    v4.append(mer4[n1+n2+n3+n4])
    v=v2+v3+v4
    return v
def ksnpf(seq):
    kn=5
    freq=[]
    v=[]
    for i in range(0,kn):
        freq.append({})
        for n1 in 'ATCGN':
            freq[i][n1]={}
            for n2 in 'ATCGN':
                freq[i][n1][n2]=0
    seq=seq.strip('N')
    seq_len=len(seq)
    for k in range(0,kn):
        for i in range(seq_len-k-1):
            n1=seq[i]
            n2=seq[i+k+1]
            freq[k][n1][n2]+=1
    for i in range(0,kn):
        for n1 in 'ATCG':
            for n2 in 'ATCG':
                v.append(freq[i][n1][n2])
    return v
def npf(seq):
    binary_dictionary={'A':[1,1,1],'T':[0,1,0],'G':[1,0,0],'C':[0,0,1],'N':[0,0,0]}
    cnt=[]
    for i in seq:
        cnt.append(binary_dictionary[i])
    return reduce(operator.add,cnt)
def ssc(seq):
    pname="D:\VRNA\RNAfold.exe"
    source=seq.replace('\n','')
    source=source+'N'
    seq=seq.replace('\n','')
    p=subprocess.Popen(pname,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    result=p.communicate(input=source)
    res=result[0].decode()[0:]
    length=len(seq)
    ssc={}
    ssc_vec={}
    for n1 in 'ATCG':
        for n2 in '.()':
            for n3 in '.()':
                for n4 in '.()':
                    ssc[n1+n2+n3+n4]=0
    res=res.split('N')
    res_str=res[1]
    res_str=res_str.encode()
    res_len=len(res_str)
    res_str=res_str[2:res_len-12]#from 2
    for p in range(0,length-2):
        ssc[seq[p]+res_str[p:p+3]]+=1
    for n1 in 'ATCG':
        ssc_vec[n1+'...']=ssc[n1+'...']
        ssc_vec[n1+'..(']=ssc[n1+'..(']+ssc[n1+'..)']
        ssc_vec[n1+'.(.']=ssc[n1+'.(.']+ssc[n1+'.).']
        ssc_vec[n1+'(..']=ssc[n1+'(..']+ssc[n1+')..']
        for n2 in '()':
            for n3 in '()':
                ssc_vec[n1+'.((']=ssc[n1+'.'+n2+n3]
        for n2 in '()':
            for n3 in '()':
                ssc_vec[n1+'(.(']=ssc[n1+n2+'.'+n3]
        for n2 in '()':
            for n3 in '()':
                ssc_vec[n1+'((.']=ssc[n1+n2+n3+'.']
        for n2 in '()':
            for n3 in '()':
                for n4 in '()':
                    ssc_vec[n1+'(((']=ssc[n1+n2+n3+n4]
    v=[]
    for n1 in 'ATCG':
        for n2 in '.(':
            for n3 in '.(':
                for n4 in '.(':
                    v.append(ssc_vec[n1+n2+n3+n4])
#    for n1 in 'ATCG':
#        for n2 in '.()':
#            for n3 in '.()':
#                for n4 in '.()':
#                    v.append(ssc[n1+n2+n3+n4])
    return v

def main():
    pseknc = PseKNC(k=2, lamada=1, w=0.05)
    pos_vec2 = pseknc.make_pseknc_vec(open('postrain.txt'))
    neg_vec2 = pseknc.make_pseknc_vec(open('negtrain.txt'))
    fea_vec2=[]
    fea_vec2.extend(pos_vec2+neg_vec2)
    psednc = PseDNC(lamada=1, w=0.05)
    pos_vec1 = psednc.make_psednc_vec(open('postrain.txt'))
    neg_vec1 = psednc.make_psednc_vec(open('negtrain.txt'))
    fea_vec1=[]
    fea_vec1.extend(pos_vec1+neg_vec1)
    feature_matrix=[]
    label_vector=[]
    train_samples=open('./data_new.txt','r')
    i=0
    for line in train_samples:
        feature_vector=[]
        if i<596:
           label_vector.append(1)
        else:
            label_vector.append(0)
        sequence=line
        feature_vector.extend(fea_vec1[i]+fea_vec2[i]+ssc(sequence))
        sequence=line.replace('\n','')
        feature_vector.extend(kmer(sequence)+ksnpf(sequence)+npf(sequence))

        feature_matrix.append(feature_vector)
        i=i+1
    train_samples.close()
    feature_array = np.array(feature_matrix,dtype=np.float32)
    min_max_scaler = preprocessing.MinMaxScaler(copy=True, feature_range=(-1, 1))
    feature_scaled= min_max_scaler.fit_transform(feature_array)
    X=feature_scaled
    y=label_vector
    X_train, X_test, y_train, y_test=train_test_split(X,y,test_size=0.4,random_state=0)
    clf = SVC(C=1.11,gamma=0.003,probability=True)
    clf.fit(X_train,y_train)
    print clf.score(X_test,y_test)
    predict_y_test = clf.predict(X_test)
        
    TP=0
    TN=0
    FP=0
    FN=0 
    for i in range(0,len(y_test)):
        if int(y_test[i])==1 and int(predict_y_test[i])==1:
            TP=TP+1
        elif int(y_test[i])==1 and int(predict_y_test[i])==0:
            FN=FN+1
        elif int(y_test[i])==0 and int(predict_y_test[i])==0:
            TN=TN+1
        elif int(y_test[i])==0 and int(predict_y_test[i])==1:
            FP=FP+1
    Sn=float(TP)/(TP+FN)
    Sp=float(TN)/(TN+FP)
    ACC=float((TP+TN))/(TP+TN+FP+FN)
    prob_predict_y_test = clf.predict_proba(X_test)
    predictions_test = prob_predict_y_test[:, 1]
#######generate combined negative scores        
        #combined_prob=predictions_test        
        
    y_validation=np.array(y_test,dtype=int)
    fpr, tpr, thresholds =metrics.roc_curve(y_validation, predictions_test,pos_label=1)
    roc_auc = auc(fpr, tpr)
        #print('AdaBoostClassifier AUC:%s'%roc_auc)
    F1=metrics.f1_score(y_validation, map(int,predict_y_test))
    MCC=metrics.matthews_corrcoef(y_validation,map(int,predict_y_test))
    print('SVM Accuracy:%s'%ACC)
    print('SVM AUC:%s'%roc_auc)
    print('SVM Sensitive:%s'%Sn)
    print('SVM Specificity:%s'%Sp)
    print('SVM F1:%s'%F1)
    print('SVM MCC:%s'%MCC)  

if __name__=='__main__':
    main()