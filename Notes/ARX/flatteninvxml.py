#!/usr/local/bin/python 
import Jupytils
import sys
from Jupytils.DataFrameUtils import *

def thetaF1(ret):
    t = np.array([ float(i.strip()) for i in  ret['theta'].split(',') if i.strip()])
    t[0:ret.n] *= -1;
    return t

def thetaF2(ret):
    t = [ float(i.strip()) for i in  ret['theta'].split(',') if i.strip()]
    return t

def thetaF1Str(ret):
    t = np.array([float(i.strip()) for i in  ret['theta'].split(',') if i.strip()])
    t[0:ret.n] *= -1;
    return str(list(t))[1:-1]

def LoadInvFile(file, needTheta1=True):
    if ( file.endswith("xml")):
        dfi1=LoadDataSet( file, xmlTag="Invariant")
        dfi1 = dfi1["uName,yName,fitness,correlation,n,m,k,threshold,theta".split(',')];
    else:
        dfi1=pd.read_csv(file)
        
    dfi1.uName = dfi1.uName.str.strip()
    dfi1.yName = dfi1.yName.str.strip()
    dfi1.n = dfi1.n.astype(int)
    dfi1.m = dfi1.m.astype(int)
    dfi1.k = dfi1.k.astype(int)
    dfi1.fitness = dfi1.fitness.astype(float)
    dfi1.threshold = dfi1.threshold.astype(float)     
   
    if ( needTheta1):
        if ( file.endswith("xml")):
            dfi1['theta1']= dfi1.apply (lambda row: thetaF1(row),axis=1)  # Convert to float array
        else:
            dfi1['theta1']= dfi1.apply (lambda row: thetaF2(row),axis=1)  # Convert to float array

    dfi1.sort_values(['uName', 'yName'], inplace=True)
    dfi1.reset_index(inplace=True, drop= True)
    return dfi1

def inJupyter():
    try:get_ipython; return True
    except: return False;


if __name__ == '__main__':
    if (not inJupyter()):
        if (len(sys.argv) < 2 ): 
            print("flatteninv.xml <inv.xml file>")
            sys.exit(1)
        print("processing ", sys.argv[1])
        df=LoadInvFile(sys.argv[1], needTheta1=False);
        df['theta'] = df.apply (lambda row: thetaF1Str(row),axis=1)
        df.to_csv( sys.argv[1]+".csv", index=False)