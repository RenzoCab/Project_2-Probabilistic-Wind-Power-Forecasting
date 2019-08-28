import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import interpolate
import datetime as dt
import pandas as pd
import os
##############
#constants
##############

# %%
# new functions to be added

#function to convert np.datetime64 to unix time, regardless of the units
def get_unixtime(dt64):
    return dt64.astype('datetime64[s]').astype('int')

#function: returns datetime64 object from a unix timestamp in seconds
def get_datetime64(time_stamp):
    return np.datetime64(int(time_stamp),'s')
get_datetime64=np.vectorize(get_datetime64)

#this function splits forecast paths into equal chunks as requested
def split_data(array_in,path_length,axis=1):
    number_paths=array_in.shape[1]//path_length
    if isinstance(number_paths, int):
        temp=np.split(array_in,number_paths, axis=axis)
        out=np.stack(temp, axis=0)
    else:
        raise ValueError('non-integer divisions of forecasts')
    return(out)

#this function get the shape of each array in a list
def get_list_array_sizes(list):
    out=[];
    for i in range(0,len(list)):
        out.append(list[i].shape)
    return(out)

# depreciated, functionality included in func: remove_partial_observations
#this function filters paths with observations less than minimuim
#WARNING: Index order may change !!! modifty it to fill with np.nan
def remove_partial_observations(list, min_threshold_observations ):
    temp = get_list_array_sizes(list)
    temp1,temp2=zip(*temp) #temp1 is selector (datetime, data), temp2 is
    #number  of observations
    complete_observations=[list[i] for i,v in enumerate(temp2) \
     if v >= min_threshold_observations]
    return(complete_observations)

#match data to forecast: input forecast(paths,2 ,72), production is time stamp and production
def match_data_to_forecast(forecast, production, min_threshold_observations):
    observ_path=[] #=np.zeros((forecast.shape[0]))
    for j in range(0,forecast.shape[0]): #forecast.shape[0]
        path_start=forecast[j,0,0];
        path_end=forecast[j,0,71]
        observation_index=np.where(np.logical_and(path_start<=production[0,:], production[0,:]<=path_end  ))
        if production[:,observation_index[0]].shape[1]>=min_threshold_observations:
            observ_path.append(production[:,observation_index[0]])
        else:  observ_path.append(np.full( (2,min_threshold_observations)  , np.nan))
    return(observ_path)


#func: interpolate forecast to match data, input must be array(path, selector, #point)
def interp_forecast_to_data(forecast, observed):
    out=np.full(forecast.shape[:2]+ (observed.shape[2],)  , np.nan)
    for i in range(0,forecast.shape[0]):
        x = forecast[i,0,:];
        y = forecast[i,1,:];
        xvals = np.linspace(x[0], x[-1], observed.shape[2]);
        yinterp = np.interp(xvals, x, y);
        out[i,0,:]=xvals
        out[i,1,:]=yinterp
    return(out)

#func: remove nans from data and its corresponding path in forecast
def remove_nan_paths_combine(forecast, observ):
    forecast_list=[]
    observ_list=[]
    datetime_list=[]
    for i in range(0,forecast.shape[0]):
        if( not any(np.isnan(observ[i,1,:])) ):
            if(all(forecast[i,0,:]==observ[i,0,:])):
                forecast_list.append(forecast[i,1,:])
                observ_list.append(observ[i,1,:])
                datetime_list.append(observ[i,0,:])
            elif (not all(np.isnan(observ[i,0,:])) ):
                raise ValueError('datetime in forecast is not\
                 matching the data in path ' + str(i))
    forecast=np.asarray(forecast_list);
    observ=np.asarray(observ_list);
    datetime=np.asarray(datetime_list);
    return(np.stack( (datetime, observ, forecast), axis=1))
# %%
