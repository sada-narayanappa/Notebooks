import numpy as np
import pandas as pd
import datetime;
import pyaf.ForecastEngine as autof


csvfile_link = "https://raw.githubusercontent.com/antoinecarme/TimeSeriesData/master/ozone-la.csv"
ozone_dataframe = pd.read_csv(csvfile_link);
ozone_dataframe['Date'] = ozone_dataframe['Month'].apply(lambda x : datetime.datetime.strptime(x, "%Y-%m"))

#ozone_dataframe.info()

lEngine_Without_Exogenous = autof.cForecastEngine()
lEngine_Without_Exogenous.train(ozone_dataframe , 'Date' , 'Ozone', 12);

lEngine_Without_Exogenous.getModelInfo()


#ozone_forecast_without_exog = lEngine_Without_Exogenous.forecast(ozone_dataframe, 12);

#ozone_forecast_without_exog.plot.line('Date', ['Ozone' , 'Ozone_Forecast', 
                                     #'Ozone_Forecast_Lower_Bound', 
                                     #'Ozone_Forecast_Upper_Bound'], grid = True, figsize=(12, 8))

