# -*- coding: utf-8 -*-
"""
Created on Sat Nov 09 21:40:35 2019

@author: 80699
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
import subprocess
from repDNA.psenac import PseDNC
from sklearn.svm import SVC



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
def binary_code(seq):
    binary_dictionary={'A':[1,1,1],'T':[0,0,1],'G':[1,0,0],'C':[0,1,0],'N':[0,0,0]}
    nucleic_dictionary={'A':0,'T':0,'C':0,'G':0,'N':0}
    cnt=[]
    p=0
    for i in seq:
        temp=[]
        p=p+1
        nucleic_dictionary[i]+=1       
        temp=list(binary_dictionary[i])
        temp.append(nucleic_dictionary[i]/p)
        cnt.append(temp)
    return reduce(operator.add,cnt)
def main():
    psednc = PseDNC(lamada=8, w=0.8)
    pos_vec1train = psednc.make_psednc_vec(open('postrain.txt'))
    neg_vec1train = psednc.make_psednc_vec(open('negtrain.txt'))
    pos_vec1test = psednc.make_psednc_vec(open('postest.txt'))
    neg_vec1test = psednc.make_psednc_vec(open('negtest.txt'))
    fea_vec1train=[]
    fea_vec1test=[]
    fea_vec1train.extend(pos_vec1train+neg_vec1train)
    fea_vec1test.extend(pos_vec1test+neg_vec1test)
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
        feature_vector.extend(fea_vec1train[i]+ssc(sequence))
        sequence=line.replace('\n','')
        feature_vector.extend(binary_code(sequence))
        feature_matrix.append(feature_vector)
        i=i+1
    train_samples.close()
    feature_array = np.array(feature_matrix,dtype=np.float32)
    min_max_scaler = preprocessing.MinMaxScaler(copy=True, feature_range=(-1, 1))
    feature_scaled= min_max_scaler.fit_transform(feature_array)
    X=feature_scaled
    y=label_vector
    clf = SVC(C=0.98,gamma=0.001,probability=True)
    clf.fit(X,y)
    
    feature_matrix=[]
    test_label_vector=[]
    test_samples=open('./data_test.txt','r')
    i=0
    for line in test_samples:
        feature_vector=[]
        if i<149:
           test_label_vector.append(1)
        else:
            test_label_vector.append(0)
        sequence=line
        feature_vector.extend(fea_vec1test[i]+ssc(sequence))
        sequence=line.replace('\n','')
        feature_vector.extend(binary_code(sequence))
        feature_matrix.append(feature_vector)
        i=i+1
    test_samples.close()
    test_feature_array = np.array(feature_matrix,dtype=np.float32)
    X_test=min_max_scaler.transform(test_feature_array)
    y_test=test_label_vector    
    
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
    a=[]
    for i in range(0,298):
        a.append(predictions_test[i])
    file=open('./score1.txt','a')
    for i in range(len(a)):
        s=str(a[i]).replace('[','').replace(']','')
        s=s.replace("'",'').replace(',','')+'\n'
        file.write(s)
    file.close
        
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