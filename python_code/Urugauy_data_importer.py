#%%
import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from data_cleansing import *
#%%


# %%
#file paths
forecast_A='/Users/alhaddwt/Google Drive/GitLab/wind-power/raw_data/Uruguay/wind/Wind_Data_2018_provider_A.csv'
forecast_B='/Users/alhaddwt/Google Drive/GitLab/wind-power/raw_data/Uruguay/wind/Wind_Data_provider_B.csv'
production_path='/Users/alhaddwt/Google Drive/GitLab/wind-power/raw_data/Uruguay/wind/anuales_gen_int_2018_10m_aggregated.csv'
#normalization_path='/Users/alhaddwt/Google Drive/GitLab/wind-power/raw_data/France/normalization_2019.txt'

#import files
d = pd.read_csv(production_path, delimiter= ',', header=2,  dtype=str)
#n = pd.read_csv(normalization_path, delimiter= '\s+', header=0)
f_A = pd.read_csv(forecast_A, delimiter= ',',  header=0, dtype={'Date': str})
f_B = pd.read_csv(forecast_B, delimiter= ',', header=0, dtype={'Date': str})

#empty containers
data=np.zeros((3,d.shape[0]))
forecast_A=np.zeros((2,f_A.shape[0]))
forecast_B=np.zeros((2,f_B.shape[0]))
production=np.zeros((2,d.shape[0]))
# %%

# %%
#convert string date and time to timestamp in seconds, strings to numerical values
#forecast A
forecast_A[0,:]=get_unixtime(pd.to_datetime(\
f_A['Unnamed: 0'] + ' ' + f_A['Unnamed: 1'], dayfirst=True ).values)
forecast_A[1,:]=pd.to_numeric(f_A['Forecast']).values

#forecast B
forecast_B[0,:]=get_unixtime(pd.to_datetime(\
f_B['Unnamed: 0'] + ' ' + f_B['Unnamed: 1'], dayfirst=True ).values)
forecast_B[1,:]=pd.to_numeric(f_B['Forecast']).values

#production
production[0,:]=get_unixtime(pd.to_datetime(\
d['Fecha'], dayfirst=False ).values)
production[1,:]=pd.to_numeric(d['Eólica']).values
# %%

# %%
# split forecasts into path
forecast_A_split=split_data(forecast_A,72)
forecast_B_split=split_data(forecast_B,72)

# match data to forecast
observ_path_A=np.stack(match_data_to_forecast(\
forecast_A_split, production, 427), axis=0)
observ_path_B=np.stack(match_data_to_forecast(\
forecast_B_split, production, 427), axis=0)

#interpolate forecast to match the resolutin of the data
forecast_A_split_interp=interp_forecast_to_data(forecast_A_split,observ_path_A )
forecast_B_split_interp=interp_forecast_to_data(forecast_B_split,observ_path_B )

#remove nan paths generated from the minmuim data points needed in the matching
forecast_data_A=remove_nan_paths_combine(forecast_A_split_interp,observ_path_A )
forecast_data_B=remove_nan_paths_combine(forecast_B_split_interp,observ_path_B )
# %%

#crude normalization:
forecast_data_A[:,1:,:]=np.abs(forecast_data_A[:,1:,:]/1500)

np.max(forecast_data_A[:,1:,:])
np.min(forecast_data_A[:,1:,:])

plt.plot(forecast_data_A[640,2,:])

# %%
base_path='/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code/data/cleansed/'
np.save(base_path+ 'URG_forecast_data_A_2018',forecast_data_A)
np.save(base_path+ 'URG_forecast_data_B_2018',forecast_data_B)
# %%
