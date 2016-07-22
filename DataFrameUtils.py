import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
import pandas as pd
import os
from IPython.display import display
from IPython.display import HTML
import dateutil;
import json;
import urllib.request;

import matplotlib
#matplotlib.style.use('ggplot')

np.set_printoptions(precision=2, linewidth=100)
pd.set_option('display.width', 1000)
pd.set_option('precision',5);
pd.set_option('display.precision', 5)
pd.set_option('mode.sim_interactive', True);
pd.set_option('display.max_rows', 9)
pd.set_option('display.float_format', lambda x: '%.3f' % x)

'''
Loads data set for person comparative analysis
Problems and ISSUES - Check the following:
==========================================
** Does your data have headers! If not you need more complex call
If so, Pass "headers=None" and pass names

'''

def DetermineSeperator(line):
    sep = ","    
    split2  = line.split("\t");
    if (len(split2) > 1 ):
        sep = "\t";   
    return sep;
    
def getAuraDF(link):
    f = urllib.request.urlopen(link)
    js = f.read().decode('UTF-8')
    if (js.find("$rs=") > 0):
        js = js[js.find("$rs=")+5:]
        #print js[0:100]
        data = json.loads(js)
        df=pd.DataFrame(data['rows'],columns=data['colnames'])
        return df
    else:
        return js
        
def getDF(fileName, debug=False, headers=0, names=None, usecols=None):
    if (    not (fileName.startswith("http://"))  and
            not (fileName.startswith("https://")) and
            not os.path.exists(fileName)):
        #raise Exception( fileName + " does not exist")
        print ("ERROR: *** " +fileName + " does not exist");
        return None;
 
    if fileName.endswith(".xlsx"): 
       df1 = pd.read_excel(fileName)
    elif ("/aura/" in fileName):
        df1 = getAuraDF(fileName);
        return df1
    else:
        sep = ","
        if not fileName.endswith(".csv"):
            with open(fileName, 'r') as f:
                line    = f.readline();    
                #split1  = line.split(",");
                sep = DetermineSeperator(line);
                
        df1 = pd.read_csv(fileName, sep=sep, header=headers, low_memory=False,
                          skipinitialspace =True, names=names, comment='#', usecols=usecols)
    return df1;

    
def LoadDataSet(fileOrString, columns=None, 
                debug=False, headers=0, names=None, checkForDateTime=False, usecols=None, seperator=None):
    if (fileOrString.find("\n") >=0 ):
        ps = [line.strip() for line in fileOrString.split('\n')
                if line.strip() != '' and not line.startswith("#") ];
        if (seperator is None):
            sep = DetermineSeperator(ps[0]);
        else:
            sep = seperator;
        ns = [p.split(sep) for p in ps]
        df1 = pd.DataFrame(ns[1:], columns=ns[0]);
    else:               
        df1 = getDF(fileOrString, debug=False, headers=0, names=None, usecols=usecols)     

    if ( df1 is None or str(type(df1)).find("DataFrame") < 0):
        return df1;
    #df1=df1.convert_objects(convert_numeric=False)
    df2 = df1[columns] if (columns != None ) else df1;

    if (checkForDateTime):
        for i, c in enumerate(df1.columns):
            if (df2.dtypes[i] != object ):
                continue;
            s = df2[c][0]
            if ( len(s) < 8): continue;
            try:
                dateutil.parser.parse(s);
            except:
                continue; 
            print ("Trying to convert to datetime:"+ c);
            df2[c] =  pd.to_datetime(df2[c])  
            
    if debug:
        print ("Printing 5 of %d rows", df2.shape[0]);
        print (df2[:5]);
        
    return df2;


def normalizeData(df):
    df1 = df.select_dtypes(exclude=[object])
    vals = df1.values
    cols = df1.columns

    d = preprocessing.scale(vals)
    df2  = pd.DataFrame(d, columns=cols)
    return df2;

#def DisplayInit():
#    style= '''   
#        <style>   
#        .container { width: 100% !important; }
#            div.cell{
#                width:100%;
#                margin-left:0%;
#                margin-right:auto;
#            }
#        </style>    ''';
#    display(HTML(style))
    
#def pcolumns(df):
#    c1 = [ k for k in df.columns]
#    print (c1)
#    return c1
    
#DisplayInit()

