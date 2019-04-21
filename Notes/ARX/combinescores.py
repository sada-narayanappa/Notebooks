#!/usr/local/bin/python 
import pandas as pd
import numpy as np;
import datetime, sys
import xml.etree.ElementTree as ET
from operator import itemgetter
import json

from operator import itemgetter
def getScoreTable(sdf, idx=0, maxRows=1000):
    cols=sdf.columns[1:];
    row=sdf.ix[0][1:]
    ret = [[cols[i], row[i]] for i in range(1,len(cols)-1) if row[i] > 0 ]
    #for ri in ret:
    #    ri.append(len(ri)) 
    ret=sorted(ret, key=itemgetter(1), reverse=True)
    return ret[0:maxRows]

def LoadBrknXML(brkInvFile, xmlTag = "RCResult"):
    with open(brkInvFile, "r") as f:
        xmlText = f.read()

    root = ET.fromstring(xmlText)
    iters= root.findall(path=xmlTag) if (xmlTag) else root

    # This will write columns as follows
    # period, anomalyScore, brokenCount, [[uname, yname, residual, significance]], []
    ret = []
    for i,e in enumerate(iters):
        (p_, a_, b_) = int(e.findtext('period')), float(e.findtext('anomalyScore')), int(e.findtext('brokenCount'))
        bi = ""
        iters1= e.findall('brokenInvariant')
        ret1=[]
        for i1,e1 in enumerate(iters1):
            (u, y, res, sig) =  e1.findtext('uName'),     e1.findtext('yName'),   \
                                e1.findtext('residual'), e1.findtext('significance');
            ret1.append((u, y, res, sig))
            
        ret1=sorted(ret1,key=itemgetter(3))

        bi = json.dumps(ret1)
        ret.append((p_, 0, a_, b_, bi,0))
        
    return ret;

def main(argv):
    if (len(argv) <= 1):
        print ("*ERROR* need prefix - expecting files in format _BrokenInv.csv,_Score.csv")
        return 
    
    pref=argv[1] #"ab1.csv_siat"

    brkn=LoadBrknXML(f"/NEC/SIAT-OLD/SIAT-OLD/benchmarks/{pref}_BrokenInv.csv")
    scoreDF=pd.read_csv(f"/NEC/SIAT-OLD/SIAT-OLD/benchmarks/{pref}_Score.csv")

    brk=pd.DataFrame(brkn)
    brk.columns='period,time,score,brknCount,brknPairs,brknSensors'.split(',')
    brk.time=scoreDF.Time
    bs = [getScoreTable(scoreDF,i) for i in range(scoreDF.shape[0])]
    brk.brknSensors = bs

    brk.to_csv(f"/NEC/SIAT-OLD/SIAT-OLD/benchmarks/{pref}.out", index=False)

    
def inJupyter():
    try:get_ipython; return True
    except: return False;

if __name__ == '__main__':
    if (not inJupyter()):
        t1 = datetime.datetime.now()
        main(sys.argv)
        t2 = datetime.datetime.now()
        print(f"All Done in {str(t2-t1)} ***")
    else:
        pass #d