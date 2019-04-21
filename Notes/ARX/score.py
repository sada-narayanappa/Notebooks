#!/usr/local/bin/python 

import logging as log
import sys
import os
import pandas as pd
import numpy as np;
import datetime
import getopt
from collections import defaultdict
from Jupytils.DataFrameUtils import LoadDataSet
from Jupytils.jcommon import Map
import json, math
import cmp
from numba import jit, autojit

#GLOBALS 
#--------------------------------------------------------------------------
def Usage(msg=''):
    print(msg)
    print(''' Usage:
Ex: python score.py -i data/test.inv.xml -a data/test.csv 
    
score.py -i <model-file> -a <anomaly-file> [-p <prefix for output files>]
        [-s <starttime date or index>] [-e <endtime date or index>]
        
-P : append to results file if it exists        
-S : just Anomaly score

-s : Optional start time to start from the anomaly file ex: "5/6/2017 10:20"
-e : Optional end time to stop in the anomaly file

Anomaly file must have a first column is the index''')

def deleteFile(f):
    try: os.remove(f) 
    except: pass;

@jit(nopython=True, cache=True)
def predict3(x, y, n,m,k, theta, t):
    yh = theta[-1];
    for i in range(n):
        yh += y[t-n+i] * theta[n-1-i]
    for i in range(m+1):
        yh += x[t-k-i] * theta[i+n]
    
    rs = (y[t] - yh)    
    return yh ,rs, None

@jit(nopython=True, cache=True)
def predict31(x, y, n,m,k, theta, t):
    ret = []
    yh = theta[-1];
    for i in range(n):
        yh += y[t-n+i] * theta[n-1-i]
        ret.append(y[t-n+i])
    for i in range(m+1):
        yh += x[t-k-i] * theta[i+n]
        ret.append(x[t-k-i])
    
    ret.append(1.0)
    
    rs = (y[t] - yh)    
    return yh ,rs, ret

@jit(nopython=True, cache=True)
def predictArray(x, y, n,m,k, theta, t):
    ret = []
    for i in range(n-1,-1,-1):
        ret.append(y[t-n+i])
    for i in range(m+1):
        ret.append(x[t-k-i])
    ret.append(1.0)
    return np.array(ret)

@jit(nopython=True, cache=True)
def predict4(x, y, n,m,k, theta, t):
    ar = predictArray(x, y, n,m,k, theta, t)
    yh= np.sum(ar*theta)
    rs= y[t] - yh
    return yh ,rs, ar

def FindScore(t, invDFA,invDFC,anomDFA,anomDFC, sumFitness, edgeCounts, logd=True,ret=None):
    if (not ret): 
        ret = Map({});
    ret.brkns = Map({});
        
    ret.t = t;
    ret.anomScore = 0;
    ret.brknCount= 0;
    ret.SumSignif = 0;
    ret.brkSignif = 0;
    ret.scores = []
    ret.sumFitness = sumFitness; #sum(invDF.fitness.astype(float));

    for i in range(len(invDFA)):
        r = invDFA[i]
        n = r[invDFC.n]; m = r[invDFC.m]; k = r[invDFC.k]
        uName     = r[invDFC.uName]; 
        yName     = r[invDFC.yName]; 
        theta1    = r[invDFC.theta1]
        threshold = r[invDFC.threshold]
        fitness   = r[invDFC.fitness]
        
        brkn= Map({})
        brkn.brknFitSum = 0;
        brkn.fitness    = fitness;
        brkn.threshold  = threshold;
        brkn['uName']   = uName; brkn['yName'] = yName; #brkn.nmk = (r.n, r.m, r.k)
        
        u=anomDFA[:, anomDFC[uName] ]
        y=anomDFA[:, anomDFC[yName] ]
        
        brkn.yh,  brkn.res, ar = predict3(u, y, n,m,k, theta1, t)
        brkn.yh1, brkn.res1,ar = predict4(u, y, n,m,k, theta1, t)
        
        if (threshold):
            brkn.sigf = abs(brkn.res)/threshold
        else:
            brkn.sigf = math.inf
            
        ret.SumSignif += brkn.sigf
        brkn.wsigf     = brkn.sigf * fitness
        brkn.broken = " "
        brkn.y = y[t]
        
        broken1= abs(brkn.res)  > 1.1 * threshold
        broken2= abs(brkn.res1) > 1.1 * threshold
                
        if ( broken1 and broken2 ):
            #print(f'{uName}, {yName}, {broken1}, {broken2}')
            ret.brkns[uName+":"+ yName] = brkn
            brkn.broken = "*"
            ret.anomScore += fitness
            ret.brknCount+= 1
            ret.brkSignif += brkn.sigf

    ret.normedScore = ret.anomScore/ret.sumFitness
    ret.suspicions, ret.suspicionsCounts = suspicions(ret, edgeCounts, logd);
    
    if ( logd ):
        print( "Score: {:3d} Anom:{:10.4f} Norm:{:5.2f} brk:{:3d}".format( \
                                 t, ret.anomScore, ret.normedScore, ret.brknCount) )
    return ret;
#--------------------------------------------------------------------------
def suspicions(ret, edgeCounts, logDetailed=True):
    bvars=defaultdict(float)
    bvarsCount=defaultdict(int)
    bfits=0;
    for i,brkn in (ret.brkns.items()):
        #print(i)
        if (logDetailed):
            print("{} {} fit: {:.4f} Sigf: {:10.4f}({:7.4f}) score: {:6.2f} - {:6.2f} = {:8.2f}, thresh: {:12.6f}".\
            format(brkn.broken,i,  brkn.fitness, brkn.sigf, brkn.wsigf, brkn.y, brkn.yh, brkn.res, brkn.threshold* 1.1) )
        if ( brkn.broken.strip()):
            bvars[brkn.uName] += brkn.fitness
            bvars[brkn.yName] += brkn.fitness
            bvarsCount[brkn.uName] += 1
            bvarsCount[brkn.yName] += 1
            #bvars[brkn.uName] += 1
            
    for i, k in sorted(bvars.items()):
        bvars[i] = k/edgeCounts[i]
        if (logDetailed): print(i, bvars[i]), edgeCounts
    return bvars, bvarsCount;

#--------------------------------------------------------------------------
'''
def formatBrokenInvs(ret, t, fdOut=None, header=False):
    sks = sorted(ret.brkns );
    vs=str(t) + "," + ",".join(["{:5.3f}".format(ret.brkns[c].sigf) for c in sks])
    
    if (fdOut is not None):
        if (header):
            ks="time," + ",".join(sks)
            fdOut.write(ks); fdOut.write("\n");

        fdOut.write(vs); fdOut.write("\n");

    return vs;

def formatRankScores(nodes, ret, t, edgeCounts=None, fdOut=None, header=False, ):
    sks = ret.suspicions;
    s = [ str("{:5.3f}".format(sks[c])) for c in sorted(nodes) ]
    vs1=str(t) + "," + ",".join(s) 
    
    #-NOTE: Following two lines can be used for more detailed info
    #vs2=",".join(["{:3.2f}/{}".format(ret.suspicionsCounts[c], edgeCounts[c]) for c in sks])
    #vs3=",".join(["{:3.2f} {:3.2f}/{}".format(ret.suspicions[c], ret.suspicionsCounts[c], edgeCounts[c]) for c in sks])

    if (fdOut is not None):
        if (header):
            ks= ",".join(sorted(nodes))
            fdOut.write(ks); fdOut.write("\n");
            
        fdOut.write(vs1);fdOut.write("\n");
            
    return vs1;
    
def FindResidualsOld(invDF, anomDF, start=None,stop=None, log=False, filePrefix="Res_"):
    maxLag = max(max(invDF.n),max(invDF.m+invDF.k))
    start  = maxLag if start is None or start < maxLag else start
    stop   = len(anomDF) if stop is None or stop >= len(anomDF) else stop
    
    allUVs   = list(invDF.uName.values) + list(invDF.yName.values)
    nodes  = np.unique(allUVs)
    edgeCounts=pd.value_counts(pd.Series(allUVs))
    edges  = np.unique(list(invDF.uName.values + ":"+ invDF.yName.values))
    sumFit = sum(invDF.fitness.astype(float))
    time1 = datetime.datetime.now()
    print("{}, maxLag:{}, #Nodes:{}, #Invariants:{}, sumFit:{}".format( \
                                    time1, maxLag, len(nodes), len(edges), sumFit))
    #ret=Map({})
    print( "Deleting files: ... " )
    deleteFile(filePrefix +"_py_AbIndicator.csv") 
    deleteFile(filePrefix +"_py_BrokenInvs.csv")
    deleteFile(filePrefix +"_py_Scores.csv")
        
    scoreFile = open(filePrefix +"_py_AbIndicator.csv", "a", 128); # Anomaly Scores
    brInvFile = open(filePrefix +"_py_BrokenInvs.csv", "a");      # Broken Invariant File 
    ranksFile = open(filePrefix +"_py_Scores.csv", "a");      # Ranking File

    scoreFile.write("time,AnomScore,NormScore,NumBrkn\n");

    nLines=0;
    r = None
    invDFA = invDF.as_matrix()
    invDFC = Map({})
    for i, c in enumerate(invDF.columns):
        invDFC[c] = i;
        
    anomDFA = anomDF.as_matrix() 
    anomDFC = Map({})
    for i, c in enumerate(anomDF.columns):
        anomDFC[c] = i;
    
    for t in range(start, stop):
        r=DetectAnomaly(t, invDFA, invDFC, anomDFA, anomDFC, sumFit, edgeCounts, log, r, )

        tm=anomDF.ix[t][0]
        scoreFile.write(f"{tm},{r.anomScore},{r.normedScore},{r.brknCount}\n")

        formatBrokenInvs(r, tm, brInvFile, not nLines);
        formatRankScores(nodes, r, tm, edgeCounts, ranksFile, not nLines);
        nLines += 1

        sys.stdout.write(f'\r Completed {nLines} of {stop-start}')
    
    time2 = datetime.datetime.now()
    print("\n\n*** ALL Done {} ... {}".format(time2, str(time2-time1)) )
    
    scoreFile.close()
    brInvFile.close()
    ranksFile.close()
    
    return edgeCounts;
'''

def formatBrokenPairs(bret):
    ret=[]
    for c in  bret.brkns:
        z= Map(bret.brkns[c])
        ret.append((z.uName, z.yName, z.res1, z.sigf))

    from operator import itemgetter
    ret=sorted(ret,key=itemgetter(3), reverse=True)
    return json.dumps(ret)
    

def formatSuspiciousSensors(nodes, r, t, edgeCounts=None ):
    r1=sorted(r.suspicions.items(), key=lambda kv: kv[1], reverse=True)
    r2=[r.suspicionsCounts[c[0]] for c in r1]
    #return json.dumps(ret)
    rd=[]
    for i in range(len(r1)):
        rd.append((r1[i][0], r1[i][1], r2[i])) 
    return json.dumps(rd)


def FindResiduals(invDF, anomDF, start=None,stop=None, log=False, filePrefix="Res_"):
    maxLag = max(max(invDF.n),max(invDF.m+invDF.k))
    start  = maxLag if start is None or start < maxLag else start
    stop   = len(anomDF) if stop is None or stop >= len(anomDF) else stop
    
    allUVs   = list(invDF.uName.values) + list(invDF.yName.values)
    nodes  = np.unique(allUVs)
    edgeCounts=pd.value_counts(pd.Series(allUVs))
    edges  = np.unique(list(invDF.uName.values + ":"+ invDF.yName.values))
    sumFit = sum(invDF.fitness.astype(float))
    time1 = datetime.datetime.now()
    print("{}, maxLag:{}, #Nodes:{}, #Invariants:{}, sumFit:{}".format( \
                                    time1, maxLag, len(nodes), len(edges), sumFit))
    deleteFile(filePrefix +"_py.out")
    scoreFile = open(filePrefix +"_py.out", "a"); # Anomaly Scores
    scoreFile.write("period,time,score,brknCount,brknPairs,brknSensors\n");

    nLines=0;
    r = None
    invDFA = invDF.as_matrix()
    invDFC = Map({})
    for i, c in enumerate(invDF.columns):
        invDFC[c] = i;
        
    anomDFA = anomDF.as_matrix() 
    anomDFC = Map({})
    for i, c in enumerate(anomDF.columns):
        anomDFC[c] = i;
    
    for t in range(start, stop):
        r=FindScore(t, invDFA, invDFC, anomDFA, anomDFC, sumFit, edgeCounts, log, r, )

        tm=anomDF.ix[t][0]
        bp=formatBrokenPairs(r);
        bs=formatSuspiciousSensors(nodes, r, tm, edgeCounts);
        scoreFile.write(f"{t},{tm},{r.anomScore},{r.brknCount},{bp},{bs}\n")

        nLines += 1
        sys.stdout.write(f'\r Completed {nLines} of {stop-start}')
    
    time2 = datetime.datetime.now()
    print("\n\n*** ALL Done {} ... {}".format(time2, str(time2-time1)) )
    
    scoreFile.close()
    return edgeCounts;
#--------------------------------------------------------------------------
def process(opts):
    startTime = opts['-s']
    endTime   = opts['-e']
    invFile   = opts['-i']
    anomFile  = opts['-a']
    filePrefix= opts['-p'] or anomFile.replace(".csv","");

    print(opts)

    if ( not os.path.exists(anomFile) ):
        print("Anomaly file '{}' does not exist: ".format(anomFile));
        return
    if ( not os.path.exists(invFile) ):
        print("Invariant file '{}' does not exist: ".format(invFile));
        return
    
    anomDF = pd.read_csv(anomFile)
    invDF  = cmp.LoadInvFile(invFile)
    print("Creating Score Files ... FilePrefix: ",  filePrefix)
    
    FindResiduals(invDF, anomDF, filePrefix= filePrefix); # Lets ignore start and stop for now
#--------------------------------------------------------------------------
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:a:p:s:e:P:S")
    except getopt.GetoptError as e:
        Usage("Getopt Error:" + str(e));
        sys.exit(2)
    options = defaultdict(str)
    
    for o,a in opts: options[o] = a;
    options['args'] = args;
    
    if ( "-h" in opts or "-i" not in options or "-a" not in options or len(opts) < 1 ): 
        Usage(f"-i and -a Arguments required:  {len(opts)}"); 
        #sys.exit(1);
    
    process(options)
    
def inJupyter():
    try:
        get_ipython
        return True
    except:
        return False

if __name__ == '__main__':
    if (not inJupyter()):
        t1 = datetime.datetime.now()
        main()
        t2 = datetime.datetime.now()
        print(f"All Done in {str(t2-t1)} ***")
    else:
        pass 
        #dfi1, dfi2 = testCompareInvCSV ("/EXTDATA/app1/test/rtest.inv.xml",'/EXTDATA/app1/test/rtest.model')
        #dfi1, dfi2 = TestCompareInvCSV ("data/test.inv.xml", "data/test.csv.inv.csv")