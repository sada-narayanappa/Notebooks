from invx import *
import logging as log
import sys
import os
import pandas as pd
import numpy as np;
import datetime
import getopt
from collections import defaultdict
import math
import re;
import collections
import gc
import json
from Jupytils import Map
import Jupytils
from Jupytils.DataFrameUtils import *

'''
This will compare their generated XML FIle to our generated XML FIle and ensure they match;
Very very important Test case
'''

def AreTwoRowsSame(dfi1, dfi2, loc =0):
#    u1,v1,f1,c1,_,t1,n1,m1,k1,_,_,_,theta1=dfi1.loc[loc]
    u1,v1,f1,c1,n1,m1,k1,th1,t1,theta1=dfi1.loc[loc]
    
    r = dfi2[(dfi2.uName == u1) & ( dfi2.yName == v1)]
    #u,v,f,c,_,t,n,m,k,thr,theta=dfi2.loc[loc]
    u2,v2,f2,c2,n2,m2,k2,th2,t2,theta2=r.values[0]
    theta2 = np.array(theta2)
    theta2[0:n2] *= -1;
    
    ret = True
    if ( not np.allclose(f1, f2) ):
        print(f"*ERROR* Location: {loc} => Fitness score not match {f1} != {f2}")
        ret = False
        
    if (n1 != n2 or m1 != m2 or k1 != k2):
        print(f"*ERROR* Location: {loc} => nmk not match {n1, m1, k1} != {n2, m2, k2}")
        ret = False

    if ( not np.allclose(theta1, theta2) ):
        print(f"*ERROR* Location: {loc} => Parameters didnt match {theta1} != {theta2}")
        ret = False

    if ( not ret):
        print (len(r), u1, v1, r)
    else:
        pass; #print(loc, ' ', end='')    
    return ret;

def LoadInvFile(file):
    if ( file.endswith("xml")):
        dfi1=LoadDataSet( file, xmlTag="Invariant")
        dfi1 = dfi1["uName,yName,fitness,correlation,n,m,k,threshold,theta".split(',')];
    else:
        dfi1=pd.read_csv(file)
        
    dfi1['theta1']= dfi1.apply (lambda row: [ float(i.strip()) for i in  row['theta'].split(',') if i.strip()], axis=1)
    dfi1.n = dfi1.n.astype(int)
    dfi1.m = dfi1.m.astype(int)
    dfi1.k = dfi1.k.astype(int)
    dfi1.fitness = dfi1.fitness.astype(float)
    dfi1.threshold = dfi1.threshold.astype(float)    
    return dfi1
        
def TestCompareInvCSV(invFile, invCSV, log=True):
    dfi1 = LoadInvFile( invFile)
    dfi2 = LoadInvFile( invCSV)


    if( len(dfi1) != len(dfi2) ):
        print(f"*MISMACTH Files have different Length: {invFile}:{len(dfi1)}, {invCSV}:{len(dfi2)}")
        return dfi1, dfi2
    else:
        print(f"* PASS: Lengths Match: {invFile}:{len(dfi1)}, {invCSV}:{len(dfi2)}")

    ret = True;
    for i in range(len(dfi1)):
        ret = AreTwoRowsSame(dfi1, dfi2, i)
        if ( not ret):
            break
            
    if (ret):
        print(f'* PASS Compared {i+1} Rows')    
    return dfi1, dfi2

def Usage():
    print('''Usage: sys.argv[0]} file1 file2 
          Ex: sys.argv[0]}  test.inv.xml test.inv.csv
          ''')

def main():
    global GLOBAL_ARGS
    args = GLOBAL_ARGS['__ARGS__']
    if (len(args) < 2 ): 
        Usage();
        print(f"**ERROR: REQUIRED Parameters not given {len(args)} *\n\n");
        return;
    file1 = args[0]
    file2 = args[1]
    dfi1, dfi2 = TestCompareInvCSV (file1, file2)
    #cFrom = int(GLOBAL_ARGS['-f']) if ('-f' in GLOBAL_ARGS) else 0
    #cTo = int(GLOBAL_ARGS['-t']) if ('-t' in GLOBAL_ARGS) else 100000
    #CreateInvariants( csvp, outp, cFrom, cTo)
    
if __name__ == '__main__':
    if (not inJupyter()):
        t1 = datetime.datetime.now()
        getargs("-hf:t:")
        main()
        t2 = datetime.datetime.now()
        print(f"All Done in {str(t2-t1)} ***")
    else:
        #dfi1, dfi2 = testCompareInvCSV ("/EXTDATA/app1/test/rtest.inv.xml",'/EXTDATA/app1/test/rtest.model')
        dfi1, dfi2 = TestCompareInvCSV ("data/test.inv.xml", "data/test.csv.inv.csv")