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
import math

#GLOBALS 
G_sargs = defaultdict(str);
#--------------------------------------------------------------------------
def Usage(msg=''):
    print(msg)
    print(''' Usage:
score.py -i <invariant.xml> -o <anomaly-file> [-p <prefix for output files>]
        [-s <starttime date or index>] [-e <endtime date or index>]
        
-P : append to results file if it exists        
-S : just Anomaly score

-s --start: Optional start time to start from the anomaly file ex: "5/6/2017 10:20"
-e --end  : Optional end time to stop in the anomaly file

Anomaly file must have a first column is the index
''')

def deleteFile(f):
    try:
        os.remove(f) 
    except:
        pass;
    
def predict1(x, y, n,m,k, theta, t, log=False):
    yh=0 
    for i in range(n):
        py = y[t-i-1]
        yh += py * theta[i] * -1.0
#        if ( log): print("yh = {} * {} * {}".format(y[t-i-1] , theta[i] , -1.0), end='')
            
    for i in range(m+1):
        px = x[t-m-k+i]
        yh += px * theta[n+m-i]
#        if ( log): print(" + {} * {} ".format(x[t-m-k+i] , theta[n+m-i]), end='')

#    if ( log):  print(" + {} ".format(theta[-1]))

    yh += theta[-1]
    return yh, y[t] - yh
    
# t = at time
def DetectAnomaly(t, invDF, anomDF, sumFitness, edgeCounts, logdetailed=True, ret=Map({}) ):
    if (not ret): 
        ret = Map({});
        ret.brkns = Map({});
        
    ret.t = t;
    ret.anomScore = 0;
    ret.brokenInvs= 0;
    ret.SumSignif = 0;
    ret.brkSignif = 0;
    ret.scores = []
    ret.sumFitness = sumFitness; #sum(invDF.fitness.astype(float));

#    for i in range(len(invDF)):
#        r = invDF.ix[i]
        
    for i, r in invDF.iterrows():
        brkn= Map({})
        brkn.brknFitSum =0;
        brkn.r = r
        u=anomDF[r.uName].values
        y=anomDF[r.yName].values
        brkn['uName'] = r.uName; brkn['yName'] = r.yName; #brkn.nmk = (r.n, r.m, r.k)
        ret.brkns[brkn.uName+":"+brkn.yName] = brkn
        
        brkn.yh, brkn.res = predict1(u, y, r.n, r.m, r.k, r.theta1, t, False)
        if (r.threshold):
            brkn.sigf  = abs(brkn.res)/r.threshold
        else:
            brkn.sigf = math.inf
            
        ret.SumSignif += brkn.sigf
        brkn.wsigf     = brkn.sigf * r.fitness
        brkn.broken = " "
        brkn.y = y[t]
        
        if (abs(brkn.res) > r.threshold * 1.1):
            brkn.broken = "*"
            ret.anomScore += r.fitness
            ret.brokenInvs+= 1
            ret.brkSignif += brkn.sigf

    ret.normedScore = ret.anomScore/ret.sumFitness
    
    ret.suspicions, ret.suspicionsCounts = suspicions(ret, edgeCounts, logdetailed);
    
    if ( logdetailed ):
        print( "Score: {:3d} Anom:{:10.4f} Norm:{:5.2f} brk:{:3d}".format( \
                                 t, ret.anomScore, ret.normedScore, ret.brokenInvs) )
    return ret;

def DetectAnomaly(t, invDFA,invDFC,anomDFA,anomDFC, sumFitness, edgeCounts, logd=True,ret=None, maxLag=5 ):
    if (not ret): 
        ret = Map({});
        ret.brkns = Map({});
        
    ret.t = t;
    ret.anomScore = 0;
    ret.brokenInvs= 0;
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
        
        u=anomDFA[:, anomDFC[uName] ]
        y=anomDFA[:, anomDFC[yName] ]
        brkn['uName'] = uName; brkn['yName'] = yName; #brkn.nmk = (r.n, r.m, r.k)
        ret.brkns[uName+":"+ yName] = brkn
        
        brkn.yh, brkn.res = predict1(u, y, n,m,k, theta1, t, False)
        
        if (threshold):
            brkn.sigf  = abs(brkn.res)/threshold
        else:
            brkn.sigf = math.inf
            
        ret.SumSignif += brkn.sigf
        brkn.wsigf     = brkn.sigf * fitness
        brkn.broken = " "
        brkn.y = y[t]
        
        if (abs(brkn.res) > threshold * 1.1):
            brkn.broken = "*"
            ret.anomScore += fitness
            ret.brokenInvs+= 1
            ret.brkSignif += brkn.sigf

    ret.normedScore = ret.anomScore/ret.sumFitness
    
    ret.suspicions, ret.suspicionsCounts = suspicions(ret, edgeCounts, logd);
    
    if ( logd ):
        print( "Score: {:3d} Anom:{:10.4f} Norm:{:5.2f} brk:{:3d}".format( \
                                 t, ret.anomScore, ret.normedScore, ret.brokenInvs) )
    return ret;

#--------------------------------------------------------------------------
def suspicions(ret, edgeCounts, logDetailed=True):
    bvars=defaultdict(int)
    bvarsCount=defaultdict(int)
    bfits=0;
    for i,brkn in sorted(ret.brkns.items()):
        #print(i)
        if (logDetailed):
            print("{} {} fit: {:.4f} Sigf: {:10.4f}({:7.4f}) score: {:6.2f} - {:6.2f} = {:8.2f}, thresh: {:12.6f}".format(brkn.broken,\
            i,  brkn.fitness, brkn.sigf, brkn.wsigf, brkn.y, brkn.yh, brkn.res, brkn.threshold* 1.1) )
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

    #ret=Map({})
    if ( "-P" not in  G_sargs ):
        print( "Deleting files: ... " )
        deleteFile(filePrefix +"_AnomScores.csv") 
        deleteFile(filePrefix +"_BrokenInvs.csv")
        deleteFile(filePrefix +"_RankScores.csv")
        
    scoreFile = open(filePrefix +"_AnomScores.csv", "a", 128); # Anomaly Scores
    brInvFile = open(filePrefix +"_BrokenInvs.csv", "a");      # Broken Invariant File 
    ranksFile = open(filePrefix +"_RankScores.csv", "a");      # Ranking File

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
        dt=0;
#       try:
        d1 = datetime.datetime.now()
        r=DetectAnomaly(t, invDFA, invDFC, anomDFA, anomDFC, sumFit, \
                        edgeCounts, log, r, maxLag=maxLag)
        d2 = datetime.datetime.now()
        dt = d2- d1;
#        except:
#            print("\n: error occured during scoring");
#            break;
            
        tm=anomDF.ix[t][0]
        o1 = "{},{},{},{}\n".format(tm, r.anomScore, r.normedScore, r.brokenInvs) 
        scoreFile.write(o1)
        #scoreFile.flush()
        
        formatBrokenInvs(r, tm, brInvFile, not nLines);
        formatRankScores(nodes, r, tm, edgeCounts, ranksFile, not nLines);
        nLines += 1
        
        sys.stdout.write('\r Completed {} of {} {}'.format(nLines, stop-start,str(dt) ))
    
    time2 = datetime.datetime.now()
    print("\n\n*** ALL Done {} ... {}".format(time2, str(time2-time1)) )
    
    scoreFile.close()
    brInvFile.close()
    ranksFile.close()
    
    return edgeCounts;
    
#--------------------------------------------------------------------------
def thetaF(ret):
    t = [ float(i.strip()) for i in  ret['theta'].split(',') if i.strip()]
    return t

def process():
    global G_sargs
    startTime = G_sargs['-s']
    endTime   = G_sargs['-e']
    invFile   = G_sargs['-i']
    anomFile  = G_sargs['-a']
    filePrefix= G_sargs['-p'] or anomFile.replace(".csv","");

    print(G_sargs)

    if ( not os.path.exists(anomFile) ):
        print("Anomaly file '{}' does not exist: ".format(anomFile));
        return
    if ( not os.path.exists(invFile) ):
        print("Invariant file '{}' does not exist: ".format(invFile));
        return
    
    anomDF = LoadDataSet(anomFile)
    
    # Load Invariant File and prepare columns
    invDF=LoadDataSet(invFile, xmlTag="Invariant")
    invDF = invDF.infer_objects()
    invDF = invDF['uName yName fitness theta n m k threshold'.split()]
    invDF.fitness  = invDF.fitness.astype(float)
    invDF.threshold= invDF.threshold.astype(float)
    invDF.n = invDF.n.astype(int)
    invDF.m = invDF.m.astype(int)
    invDF.k = invDF.k.astype(int)
    invDF['theta1']= invDF.apply (lambda row: thetaF(row),axis=1)  # Convert to float array
    
    print("Creating Score Files ... FilePrefix: ",  filePrefix)
    
    FindResiduals(invDF, anomDF, filePrefix= filePrefix); # Lets ignore start and stop for now
    return '';

#--------------------------------------------------------------------------
def main():
    global G_sargs
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:a:p:s:e:P:S",["start=","end="])
    except getopt.GetoptError as e:
        Usage("Getopt Error:" + str(e));
        sys.exit(2)
    opth=False;    
    for opt, arg in opts:
        if opt == '-h': opth = True; 
        elif opt in ("--start"):  G_sargs['-s']=arg;
        elif opt in ("--end"):    G_sargs['-e']=arg;
        else:          
            G_sargs[opt]=arg;

    if ( opth ): 
        Usage("Arguments required: "); sys.exit(1);
        
    process()
    
if __name__ == '__main__':
    main()