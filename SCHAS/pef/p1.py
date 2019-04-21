import numpy as np
import pandas as pd
import datetime;

csvfile_link = "https://raw.githubusercontent.com/antoinecarme/pyaf/master/data/ozone-la-exogenous-2.csv"
exog_dataframe = pd.read_csv(csvfile_link);
exog_dataframe['Date'] = exog_dataframe['Date'].astype(np.datetime64);

#exog_dataframe

encoded_csvfile_link = "https://raw.githubusercontent.com/antoinecarme/pyaf/master/data/ozone_exogenous_encoded.csv"
encoded_ozone_dataframe = pd.read_csv(encoded_csvfile_link);
# print(encoded_ozone_dataframe.columns)
interesting_Columns = ['Date',
                       'Ozone', 'Exog2', 
                       'Exog2',
                       'Exog3', 
                       'Exog4=E', 'Exog4=F', 'Exog4=C', 'Exog4=D', 'Exog4=B',
                       'Exog5=K', 'Exog5=L', 'Exog5=M', 'Exog5=N'];

encoded_ozone_dataframe = encoded_ozone_dataframe[interesting_Columns]
#print(encoded_ozone_dataframe.info())
encoded_ozone_dataframe.head()

import pyaf.ForecastEngine as autof
lEngine = autof.cForecastEngine()

csvfile_link = "https://raw.githubusercontent.com/antoinecarme/TimeSeriesData/master/ozone-la.csv"
ozone_dataframe = pd.read_csv(csvfile_link);
ozone_dataframe['Date'] = ozone_dataframe['Month'].apply(lambda x : datetime.datetime.strptime(x, "%Y-%m"))

ozone_dataframe.info()

lExogenousData = (exog_dataframe , ['Exog2' , 'Exog3' , 'Exog4',  'Exog5']) 
#lEngine.train(ozone_dataframe , 'Date' , 'Ozone', 12 , lExogenousData);
#lEngine.getModelInfo()
lEngine_Without_Exogenous = autof.cForecastEngine()
lEngine_Without_Exogenous.train(ozone_dataframe , 'Date' , 'Ozone', 12);

lEngine_Without_Exogenous.getModelInfo()


ozone_forecast_without_exog = lEngine_Without_Exogenous.forecast(ozone_dataframe, 12);
ozone_forecast_with_exog = lEngine.forecast(ozone_dataframe, 12);

ozone_forecast_without_exog.plot.line('Date', ['Ozone' , 'Ozone_Forecast', 
                                             'Ozone_Forecast_Lower_Bound', 
                                             'Ozone_Forecast_Upper_Bound'], grid = True, figsize=(12, 8))
ozone_forecast_with_exog.plot.line('Date', ['Ozone' , 'Ozone_Forecast', 
                                             'Ozone_Forecast_Lower_Bound', 
                                             'Ozone_Forecast_Upper_Bound'], grid = True, figsize=(12, 8))

