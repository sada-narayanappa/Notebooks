import numpy as np
#import pylab as pl
#from matplotlib.colors import ListedColormap
#from sklearn import neighbors, datasets
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import HTML
import os
import glob;
import re;
plt.style.use('fivethirtyeight') # Good looking plots
pd.set_option('display.max_columns', None) # Display any number of columns
import platform
import matplotlib
from sklearn import neighbors, datasets, cluster, preprocessing, decomposition
from sklearn.decomposition import PCA
from Features import prepareDF
from matplotlib import colors


if (platform == "Windows"):
	from win32com.client import Dispatch
	from win32com.client.gencache import EnsureDispatch
	from win32com.client import constants
	from IPython.display import IFrame

np.set_printoptions(precision=2, linewidth=1000)

matplotlib.rcParams['lines.linewidth'] = 1
matplotlib.rcParams['lines.color'] = 'r'

#=======================================================================================
def Excel2Html(file, overwrite=True, show=True, leaveItOpen = True, 
                   width="100%", length="400px"):
    f = file.replace("/", "\\");
    #xl = Dispatch('Excel.Application')
    xl = EnsureDispatch ("Excel.Application")
    cwd=os.getcwd() + "\\" + f
    ext = cwd.split(".")[-1]
    nef=cwd.replace(ext, "html")
    nhtml = file.replace(ext, "html")

    fileOpenNow = False;
#    if (not os.path.exists(cwd)):
#        display("File " + cwd + " not found");
#        return;
#    else:
#       try:
#           f = open(cwd, "r+");
#           f.close();
#       except:
#           fileOpenNow =True;
    print (cwd, nef, nhtml, fileOpenNow)
           
    wb=xl.Workbooks.Open(cwd)
    xl.Visible = True #-- optional
    if os.path.exists(nef):
        os.remove(nef)
    wb.SaveAs(nef, constants.xlHtml)
    wb.Close()
#    if (fileOpenNow or leaveItOpen):
    wb=xl.Workbooks.Open(cwd)
    del xl;
    if (show):
        display(IFrame(nhtml, width, length))    

#=======================================================================================
# EXAMPLE USE
# graph a function 
# graphFunction(lambda x: x ** 3, 0,10)
#
# graphFunction('x', 0,1, "r", "$x$", "$z$", "", '.') #, 'o', "sada")
# graphFunction('1/(1+np.exp(-x))', -6,6, "g", "$x$", "$z$", "",'.',"$z=\\frac{1}{1+e^{-\\theta^{T} x}}$") #, 'o', "sada")
# graphFunction(lambda x: log(x), 0,6, "b", "$x$", "$z$", "",'.',"$log(x)$") #, 'o', "sada")

def graphFunction(formula, xmin, xmax, c=None, xlabel= None, ylabel=None, title=None, marker=None, 
		  label=None, legend=True, legendLoc=2):
    x = np.linspace(xmin, xmax, 100)
    if ( callable(formula)):
        y = np.apply_along_axis(formula, 0, x)
        #print ("Evaluating Function", str(formula), y)
    elif (type(formula) == str):
        y = eval(formula)
        #print ("Evaluating", formula, y)
    else:
        y = formula

    #return 
    label = label if label else str(formula);
    plt.plot(x, y, c=c, marker=marker, label=label, linewidth=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title (title, fontsize = 15)
    fs = 20 if (label.find("$") > 0) else 15
    plt.legend(fontsize=15, loc=legendLoc) if ( legend) else None;

'''
=======================================================================================
Convenient class A to manage and display matrices
=======================================================================================
'''
from IPython.display import display, Math, Latex
import numpy as np
import re
class A:
    print_width = 120;

    def dtype(self):
        return "A"

    def __init__(self, s = "", name =""):
        self.type = "ClassA"

        if ( type(s) is str):
            self.a = A.setM(s)
        elif ( type(s) is ndarray):
            self.a = s
        elif ( type(s) is list ):
            self.a = np.array(s) #[None]
        self.name = name;

    @staticmethod
    def setM( str = None):
        str = str.replace("[", "")
        str = str.replace("]", "")
        ss = str.split(';')
        r = None;
        for i,s in enumerate(ss):
            s = s.strip();
            if ( len(s) <= 0):
                continue;
            a = A.setA(s)
            r = a if i == 0 else np.vstack((r,a))

        return r

    @staticmethod
    def setA(s = None):
        if (s == None):
            return []
        s = s.replace(",", " ")
        a = s.split()
        try:
            a = map(int, a)
        except:
            a = map(float, a)

        a = np.array(a)
        return a[None]

    def p(self):
        np.set_printoptions(precision=2, linewidth=180)
        print ("shape=", self.a.shape, " name:", self.name)
        s = str(self.a)
        s = s.replace('[', '')
        s = s.replace(']', '')

        s =  "[" + s + "]"
        return s

    @staticmethod
    def M(m, name="", call_display=True, showdim=True):
        np.set_printoptions(precision=2, linewidth=180)
        name = name + " =" if name != "" else ""
        dim = "";
        if (showdim):
            dim = " \\times ".join(map(str, (m.shape) )) ;
            dim = " \\times ".join(map(str, (m.shape) )) ;
        s = str(m)
        s = s.replace('[', '')
        s = s.replace(']', '')
        s = s.replace('\n', '\\\\\\\\<NEW-LINE>')
        s = re.sub( '\s+', ' ', s ).strip()
        s = s.replace('<NEW-LINE>', "\n")
        s = re.sub('\n\s+', '', s)
        s = s.replace(' ', ' & ')
        s = name + "\\begin{bmatrix}\n" + s + "\n\\end{bmatrix}" +  dim + "\n"
        #print self.a
        if ( call_display):
            display(Math(s))
        return s;

    @staticmethod
    def display(*M):
        s = ""
        for m in M:
            s+= A.M(m, call_display=False, showdim=False);
        display(Math(s))


    @property
    def l(self):
        np.set_printoptions(precision=2, linewidth=180)
        name = self.name + " =" if self.name != "" else ""
        dim = " \\times ".join(map(str, (self.a.shape) )) ;
        s = str(self.a)
        s = s.replace('[', '')
        s = s.replace(']', '')
        s = s.replace('\n', '\\\\\\\\<NEW-LINE>')
        s = re.sub( '\s+', ' ', s ).strip()
        s = s.replace('<NEW-LINE>', "\n")
        s = re.sub('\n\s+', '', s)
        s = s.replace(' ', ' & ')
        s = name + "\\begin{bmatrix}\n" + s + "\n\\end{bmatrix}" +  dim + "\n"
        #print self.a
        display(Math(s))
        return s;

    def __str__(self):
        return self.p()

    def __add__(self, o):
        if (type(o) is int):
            r = self.a + 3
        else:
            r = np.add(self.a, o.a)
        ret = A(r)
        return ret


'''
# Example Usage
a= A(('1 7 2 3. 5 5 6 6 7 8.' * 3 +";") * 4, "A")
#a.a = a.a.T # create transpose
a.l          # display latex version of the matrix
'''


'''
=======================================================================================
Print analysis on the Pandas dataframe
=======================================================================================
'''
def searchDF(df, s="", cols=[], maxRows=10):
    rows=[]    
    for i, r in df.iterrows():
        for j, c in r.iteritems():
            if (len(cols) > 0 and not c in cols):
                break;
            if ( str(c).find(s) >= 0):
                rows.append(i)
                break;
        maxRows = maxRows -1;
        if (maxRows ==0):
            break
    df1 = df.iloc[rows];
    return df1;
    

def colTypesDF(df):
    dfTypes=pd.DataFrame(df.dtypes)
    dfTypes=dfTypes.transpose()
    cols=[];
    for c in dfTypes.columns:
        nc = (str(c) + "\n\t(" + str(dfTypes.iloc[0][c])+")")
        cols.append(nc);
    return cols;

def getFileName(df,idx):
    c = df.columns[idx];
    c=re.sub(r'[\s+?\']', '', c)
    prefix = hex(id(df)) + "-"+c;
    figName1 = "temp/"+ prefix + ".png";    
    
    if not os.path.exists("temp"):
        os.mkdir("temp")
    
    return figName1;    

def isAnyDataPoint(t):  
    for k in t:
        if not np.isnan(k):
            return True;
    return False;
    
def createIcon(df,idx):
    t = df.iloc[:,idx]

    scale=7;
    figName1 = getFileName(df,idx)

    if (os.path.exists(figName1)):
        return figName1;

    k = "hist"
    if (df.dtypes[idx] == object or str(df.dtypes[idx]).find("date") >= 0 ):
        u=len(t.unique());
        if (u>100):
            return None;
        else:
            t = t.value_counts();
            k = "bar";
    if not isAnyDataPoint(t):  
       return None;
               
    #print ("N= " , df.columns[idx]);
    ax=t.plot(kind=k, figsize=(1*scale, 0.5*scale), grid=True);
    ax.get_xaxis().set_visible(True)
    ax.get_yaxis().set_visible(True)
    ax.set_frame_on(False)
    fig = ax.get_figure()
    fig.savefig(figName1,  transparent=True);
        
    plt.close();
        
    return figName1;

def getIcons(df,h):
    h1="<tr><td></td>";
    idx=0;
    for c in df.columns:
        try:
           fig = createIcon(df,idx);
        except:
           print ("Error while getting icon for ", c);
        if ( fig):
            #fig = "/files/" + fig;
            h1 = h1 + "<td><a class='thumbnail' href='#thumb'><img src='" + fig;
            h1 = h1 + "' border=0 style='{margins: 0;}' width=64 height=64 ";
            h1 = h1 + "/> <span><img src='"+ fig + "' /><br /></span></a>";
#            h1 = h1 +"onmouseover='this.width=500;' onmouseout='this.width=64' >" 
                        
            h1 = h1 + "</td>";  
        else:
            h1 = h1 + "<td></td>";
            
        idx=idx + 1;          
    h1 = h1 + "</tr>\n";
    idx = h.find("<tr");
    ret = h[:idx] +h1+h[idx:];
    return ret;    

def addDescribe(df,h):
    df1d= df.describe(include='all')
    hd = df1d.to_html();
    idx1 = hd.find("<tbody>") + 8;
    idx2 = hd.find("</tbody>");
    rep = hd[idx1:idx2] 
    rep = rep.replace("<tr>", "<tr bgcolor=#e6e6fa>",4)
    rep = rep.replace("<tr>", "<tr bgcolor=#dddddd>")
    rep = rep.replace("nan", "-")
    rep = rep.replace("NaN", "-")
    idx = h.find("<tr");
    ret = h[:idx] +rep+h[idx:];
    ret = ret.replace("<th>", "<th bgcolor=#e6e6fa>",4)
    ret = ret.replace("<th>", "<th bgcolor=#6495ed>")
   
    return ret;
    
def displayDFs(dfs, maxrows = 6, showTypes = True, showIcons=True, 
               search=None, cols=[],  showStats = True):
                   
    if ( type(dfs) !=list and type(dfs) != tuple):
        dfs = [dfs];
        
    otr = "<table>"
    bg1="#efefef";
    bg2="lightblue";
    bg = bg2;
    for i, nd in enumerate(dfs):
        if ( nd is None  ):
            otr += "<td> None</td><td>&nbsp;</td>"
            continue;
            
        bg = bg2 if ( bg == bg1 ) else bg1;
        dim = str(nd.shape[0]) + " rows x " + str(nd.shape[1]) + " columns";
        d = nd[:maxrows] if (not search) else searchDF(nd,search,cols);
            
        if (showTypes):
            cols=colTypesDF(d);
            d.columns = cols
        h = d.to_html();
        
        shIcons = showIcons; 
        if ( type(showIcons) ==list and type(showIcons) != tuple):
            shIcons = showIcons[i];

        if (shIcons and nd.shape[0] > 0):
            h = getIcons(nd,h);
        if (showStats and nd.shape[0] > 0):
            h = addDescribe(nd,h);
            
        otr += "<td bgcolor=" + bg + ">" + dim + "<br>" + h + "</td><td>&nbsp;</td>"
    otr += "</table>"
    display(HTML(otr))

def formatContent(c):       
    c1 = str(c).lower().strip()
    g=[k.lower().strip() for k in "complete, finished, success, Yes".split(",") ]
    r=[k.lower().strip() for k in "error, err, failed, no".split(",") ]
    y=[k.lower().strip() for k in "pending, ongoing, current".split(",") ]
    s = ""
    if (c1 in g):
         s = "bgcolor=lightgreen";
    elif (c1 in r):
         s = "bgcolor=#FFAEAE";
    elif (c1 in y):
         s = "bgcolor=lightyellow";
        
    return "<td " + s + ">" + str(c) + "</td>"
    
#if not os.path.exists("temp"):
#    os.mkdir("temp")

#Lets Clean up before we start
[os.unlink(f) for f in glob.glob("./temp/*.png")]
    
def PCAPlot(dfL, predictColumn, s =10):
    predictColumnIdx = predictColumn+'_idx'
    
    ny = dfL[predictColumn]
    df = prepareDF(dfL, makeCopy=True)
    df = df.drop(predictColumn, axis=1)
    pca= PCA(n_components= 2)
    pca.fit(df)
    nX = pca.transform(df)
    
    le = preprocessing.LabelEncoder()
    labels = le.fit_transform(ny)
    le.classes_
    
    nDf = pd.DataFrame(nX)
    nDf[predictColumn] = ny;
    nDf[predictColumnIdx] = labels;
    
    c="r,g,b,c,m,y,k,w".split(",")
    
    for i,j in enumerate(le.classes_):
        dd = nDf[nDf[predictColumnIdx] == i]
        ll = str(le.classes_[i]);
        lb = ll + ":" + str(i) if ( ll != str(i)) else str(ll);
        #print( i,j )
        #plt.scatter(dd[[0]], dd[[1]], s=40, c=dd[predictColumnIdx].apply(lambda x:c[x]));
        plt.scatter(dd[[0]], dd[[1]], s=40, c=dd[predictColumnIdx].apply(lambda x:c[x%len(c)]), label=lb);
    
    plt.legend();

    return nDf;


#
# x must be in percentages of use - otherwise convert to percentages before this call
#
def plotPercentHist(x, bins=10, rangeI=(0.0,1.0)):
		  h, be = np.histogram(x,bins=bins,range=rangeI, normed=True)
		  plt.bar(be[:-1],h*100,width=be[1])
		  plt.xticks(be);
		  return h,be

def plotBar(x, labels=None, bottom=None):
	plt.bar(x, range(0, len(x)), bottom=bottom);
	if (labels is None):
		labels = range(0,len(x));
	plt.xticks(labels);




