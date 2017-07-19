import os
import os.path
import sys
import re;
import json;
import os
import numpy as np
from IPython.display import HTML
import pandas as pd
from datetime import timedelta;
import datetime;
from random import randint
from collections import defaultdict
from pylab import rcParams

import matplotlib.pyplot as plt
import pandas as pd
import glob;
import re;
import platform
import matplotlib

import dateutil;
import json;
import urllib.request;

sys.path.append("../../PyUtils/")

from  DataFrameUtils import *


#get_ipython().magic('run "../../PyUtils/common.ipynb"')
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
from sklearn.svm import SVR
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.gaussian_process import GaussianProcess
from sklearn import linear_model
from sklearn import cross_validation
from sklearn import ensemble
from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
#from statsmodels.regression.quantile_regression import QuantReg
#import statsmodels.formula.api as smf
from sklearn.model_selection import learning_curve
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import r2_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV

hand="http://www.smartconnectedhealth.org/aura/webroot/db.jsp?q=SELECT%20*%20FROM%20hand"
kimj="http://www.smartconnectedhealth.org/aura/webroot/db.jsp?q=SELECT%20*%20FROM%20kimj";
leeh="http://www.smartconnectedhealth.org/aura/webroot/db.jsp?q=SELECT%20*%20FROM%20leeh";
sony="http://www.smartconnectedhealth.org/aura/webroot/db.jsp?q=SELECT%20*%20FROM%20sony";
kuky="http://www.smartconnectedhealth.org/aura/webroot/db.jsp?q=SELECT%20*%20FROM%20kuky";

fileName = hand;

dfOriginal = LoadDataSet(fileName, checkForDateTime=False);
#displayDFs (dfOriginal, maxrows=3 )
#for c in dfOriginal.columns: print (c,  end=', ')


# In[3]:

df=dfOriginal.copy()

#1. Remove all the rows that does not have any pef values 
df.insert(4, 'npt',0)
df.insert(4, 'pef',0)
df.npt = pd.notnull(df.pef1)*1 + pd.notnull(df.pef2)*1 + pd.notnull(df.pef3)*1
df = df[df.npt > 0]
df=df.reset_index(drop=True)
df.pef = (df.pef1 + df.pef2 + df.pef3)/df.npt

#2. Combine 
df=df.fillna(0)
df.loc[df.timeofday == 'null', 'timeofday'] = "00:00:00"
sdttm = df.dateofmeasure + " " + df.timeofday 
df.insert(0, 'sdttm', sdttm)
df.sdttm = pd.to_datetime(df.sdttm)

df.sort_values(by='sdttm', ascending=True, inplace=True)
drps  = "cname, ampm, timeofday, dateofmeasure, npt, pef1, pef2, pef3, pef, indexpef".split(', ')
df=df.drop(drps, axis=1, errors='ignore')
df=df.reset_index(drop=True)

## <== do the following for SAP PA tool
#pef = df.pef;
#df=df.drop(['pef'], axis=1)
#df.insert(1, 'pef1', pef)

#pef[-10:]=0     # Set last 10 to zero for predictions
#df.insert(1, 'pef', pef)
#df.to_csv("HanD/hand1.csv", sep=';')

df1 = df.copy()
df1 = df1.set_index(df1.sdttm)
df1=df1.drop('sdttm', axis=1, errors='ignore')

s= pd.qcut(df1.pefmax, 10, labels='a1,a2,a3,a4,a5,a6,a7,a8,a9,a10'.split(','))
df1.insert(1,"pefcat", s)
df1.pefcat = df1.pefcat.astype(str)

columns = '''pefmax,so2,co,o3,no2,temperaturec,windspeedms,precipitationpercent,vaporpressurehpa,
dewpointtemperaturec,airpressurehpa,sealevelpressurehpa,groundtemperaturec,tmax,amax,tmin'''
#amin,pmin,tmaxlesstmin,amaxlessamin,pmaxlesspmin'''
cs = [c.strip() for c in columns.split(',')]

df1=df1[cs]
#HTML(df1.to_html())
#displayDFs (df1, maxrows=3 )


# In[4]:

df2 = df.copy()
df2=df2.rename(columns={"sdttm": "Date"})

columns = '''Date,pefmax,so2,co,o3,no2,temperaturec,windspeedms,precipitationpercent,vaporpressurehpa,
dewpointtemperaturec,airpressurehpa,sealevelpressurehpa,groundtemperaturec,tmax,amax,tmin'''
#amin,pmin,tmaxlesstmin,amaxlessamin,pmaxlesspmin'''
cs = [c.strip() for c in columns.split(',')]

df2=df2[cs]
#df2.dtypes


# In[6]:

import pyaf.ForecastEngine as autof
lEngine = autof.cForecastEngine()

lExogenousData = (df2 , cs[2:6]) 
df3=df2["Date,pefmax".split(',')].copy()


# In[ ]:

lEngine.train(df3, 'Date' , 'pefmax', 12 , lExogenousData);


# In[7]:

import pyaf.ForecastEngine as autof
lEngine_Without_Exogenous = autof.cForecastEngine()
lEngine_Without_Exogenous.train(df3 , 'Date' , 'pefmax', 7);


# In[ ]:

ozone_forecast_without_exog = lEngine_Without_Exogenous.forecast(ozone_dataframe, 12);
ozone_forecast_with_exog = lEngine.forecast(ozone_dataframe, 12);


# In[ ]:

get_ipython().magic('matplotlib inline')
ozone_forecast_without_exog.plot.line('Date', ['Ozone' , 'Ozone_Forecast', 
                                             'Ozone_Forecast_Lower_Bound', 
                                             'Ozone_Forecast_Upper_Bound'], grid = True, figsize=(12, 8))
ozone_forecast_with_exog.plot.line('Date', ['Ozone' , 'Ozone_Forecast', 
                                             'Ozone_Forecast_Lower_Bound', 
                                             'Ozone_Forecast_Upper_Bound'], grid = True, figsize=(12, 8))


