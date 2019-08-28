#%%
import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base import *
import pandas as pd
import datetime as dt
#%%

#function to convert np.datetime64 to unix time, regardless of the units
def get_unixtime(dt64):
    return dt64.astype('datetime64[s]').astype('int')

# import data and create empty containers
# %%
forecast_path='/Users/alhaddwt/Google Drive/GitLab/wind-power/raw_data/France/Predictions/forecast_2018.txt'
production_path='/Users/alhaddwt/Google Drive/GitLab/wind-power/raw_data/France/Aggregated_Production/production_2018.xls'
normalization_path='/Users/alhaddwt/Google Drive/GitLab/wind-power/raw_data/France/normalization_2019.txt'

d = pd.read_excel(production_path, delimiter= '\s+', header=None,  dtype=str)
n = pd.read_csv(normalization_path, delimiter= '\s+', header=0)
f = pd.read_csv(forecast_path, delimiter= '\s+', dtype={'Date': str})

data=np.zeros((3,f.shape[0]))
forecast=np.zeros((2,f.shape[0]))
production=np.zeros((2,d.shape[0]))

# %%

#prepare forecast, parse it and clean it
# %%
forecast[0,:]=get_unixtime(pd.to_datetime(f['Date']+' ' +f['Heure'] ).values)
forecast[1,:]=f['PrÈvision'].values
# %%

#prepare real production data, parse it and clean it
# %%
D=d.values
production_list=[]
production_time_list=[]
i=0
while i < d.shape[0]:
    if D[i,0].find('DonnÈes de rÈalisation du ')==0:
        date_str=D[i,0][26:]
        while D[i+2,0] !='nan' and i < d.shape[0]:
            time_stamp=pd.to_datetime(date_str + ' ' + D[i+2,0][:5]).timestamp()
            #reverse_time_stamp=np.datetime64(int(datetime_object),'s')
            #production[0,i]=time_stamp
            if D[i+2,11]=='*':
                pass #production[1,i]=float('nan')
            else:
                production_list.append(D[i+2,11])
                production_time_list.append(int(time_stamp))
                #production[1,i]=D[i+2,11]
            i+=1
            if i+2 >= d.shape[0]: break;
    i+=1
production=np.vstack((np.asanyarray(production_time_list, dtype=int),\
np.asanyarray(production_list, dtype=float)))
# %%



#np.datetime64(int(production[0,0]),'s')
#np.datetime64(int(forecast[0,0]),'s')


#it remains to match the forecast with the production data
# then we normalize the data set using the imported normalization data set
