# m7Gpredictor
An improved machine-learning based prediction tool for identifying RNA m7G modifications 

m7Gpredictor was implemented in python 2.7


# trainm7Gpredictor.py

#This script generates m7Gpredictor model file.

Usage:trainm7Gpredictor.py data_train.txt postrain.txt negtrain.txt

The m7Gpredictor trained with the Support Vector Machine classifier.


# testm7Gpredictor.py

#This script evaluates the model performance on the indepedent testing dataset.

Usage:testm7Gpredictor.py data_test.txt postrain.txt negtest.txt data_train.txt postrain.txt negtrain.txt

The m7Gpredictor trained with the Support Vector Machine classifier.

Evaluate the performance of m7Gpredictor on the indepedent testing dataset.


# trainiRNA-m7G.py 

#This script generates iRNAm7G model file.

Usage:trainiRNAm7G.py data_train.txt postrain.txt negtrain.txt

The iRNAm7G trained with the Support Vector Machine classifier.


# testiRNAm7G.py

#This script evaluates the model performance on the indepedent testing dataset.

Usage:testiRNAm7G.py data_test.txt postest.txt negtest.txt data_train.txt postrain.txt negtrain.txt

The iRNAm7G trained with the Support Vector Machine classifier.

Evaluate the performance of iRNAm7G on the indepedent testing dataset.
