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
        dfi1=pd.read_csv(file, comment='#')
        
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

def filterOLD(dfi, eith, both):
    ignore = {};
    keep = {};
    print(f"## ==> processing file: {sys.argv[1]}, either:{ei} both:{bo}")
    for i in range(0, len(df)):
        u1,v1,f1,c1,n1,m1,k1,th1,t1=dfi.loc[i]
        kk1 = u1+","+v1
        kk2 = v1+","+u1

        if ( f1 < both):
            ignore[kk1] = f1
            ignore[kk2] = f1
            continue;
            
        r = dfi[(dfi.uName == v1) & ( dfi.yName == u1)]
            
        if (len(r.values) < 1):
            if (f1 > eith):
                print(f"No {k2} : keeping {k1}")
                keep[k1]=f1
                continue;
            
        u2,v2,f2,c2,n2,m2,k2,th2,t2 = r.values[0]

        if ( kk1 in ignore or kk2 in ignore or kk1 in keep or kk2 in keep):
            continue;
        
        if ( f1 < both or f2 < both):
            ignore[kk1]=1; ignore[kk2]=1
            continue;
        if ( f1 < eith and f2 < eith):
            ignore[kk1]=1; ignore[kk2]=1
            continue;
            
        if (f1 > f2):
            #print(f"+ kk1 Keeping {kk1} {kk2} {f1} {f2}")
            keep[kk1]=[f1, f2]
        else:
            #print(f"+ kk2 Keeping {kk2} {kk1} {f1} {f2}")
            keep[kk2]=[f2, f1]
     
    print(f"uName,yName,f1,f2,n #: {len(keep)}")
    c=0
    for i in keep:
        #    u1,v1,f1,c1,n1,m1,k1,th1,t1=dfi.loc[i]
        print(f"{i},{keep[i][0]},{keep[i][1]},{len(keep)}")
        c +=1
        pass;
    print(f"## Done {c} rows")
    
    
def filter(dfi, eith, both):
    keep = {};
    questionables = {};
    print(f"## ==> processing file: {sys.argv[1]}, either:{ei} both:{bo}")
    for i, r in df.iterrows():
        u1,v1,f1,c1,n1,m1,k1,th1,t1=r
        uv = u1+","+v1
        vu = v1+","+u1

        if (vu not in questionables):
            questionables[uv] = f1
            continue;
            
        f2 = questionables[vu]
        del questionables[vu]
        if ( f1 < both or f2 < both):
            pass # Remove uv, vu
        elif ( f1 < ei and f2 < ei):
            pass # Remove uv, vu
        elif ( f1 > f2):
            keep[uv]=[f1,f2];
        else:
            keep[vu]=[f2,f1];
         
    for k,v in questionables.items():
        if (v >= eith):
            keep[k] = [v, -1];
            
    for i, r in df.iterrows():
        u1,v1,f1,c1,n1,m1,k1,th1,t1=r
        uv = u1+","+v1
        if ( uv not in keep):
            df.at[i,'uName'] = "#" + u1

    #print(f"uName,yName,f1,f2,n,idx #: {len(keep)}")
    c=0
    for k,v in keep.items():
        #    u1,v1,f1,c1,n1,m1,k1,th1,t1=dfi.loc[i]
        # print(f"{k},{v[0]},{v[1]},{len(keep)} {c}")
        c +=1
        u,v = k.split(",")
    print(f"## Done {c} rows")    
    
def inJupyter():
    try:get_ipython; return True
    except: return False;

if __name__ == '__main__':
    if (not inJupyter()):
        if (len(sys.argv) < 4 ): 
            print("sys.argv[0] model <either:real> <both:real>; ex: argv[0] 0.7 0.5 ")
            sys.exit(1)
        ei, bo = float(sys.argv[2]), float(sys.argv[3])
        if( bo > ei):
            print("both cannnot be greater than either {bo} > {ei}")
            sys.exit(10)
        df=LoadInvFile(sys.argv[1], needTheta1=False);
        filter(df, ei, bo);
        if ( sys.argv[1].endswith("inv.xml")):
            sfile = sys.argv[1].replace(".inv.xml", f"._f{bo}_{ei}.inv.xml.csv");
        else:
            sfile = sys.argv[1]+ f".{bo}_{ei}.filtered.csv"
        df.to_csv( sfile, index=False)
        print(f"Saved to {sfile} to remove lines run:\n\tgrep -v '#' {sfile} | wc -l ")