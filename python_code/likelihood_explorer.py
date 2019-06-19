
import os
os.chdir('/Users/alhaddwt/Waleed Dropbox/Waleed Al-Haddad/GitLab/wind-power/python_code')
from Base import *
#from multi_path_base import model_modified_drift
forecast_with_data=np.load('data/forecast_with_data.npy')
forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='likelihood_explorer'
current_plotting_dir='data/'+ script_name+'_' +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

N=forecast_data_inter.shape[2]
M= 10    #forecast_data_inter.shape[1]
dt=1
dN=1/N

disct_temp = disct(N,dt,M)
p=forecast_data_inter[0,:,:]
V= forecast_data_inter[1,:,:] -forecast_data_inter[0,:,:]

forecast_data_inter.shape
this_model=model_modified_drift(disct_temp,V, forecast= p)

theta_array=np.linspace(1, 30, 30)
alpha_array=np.linspace(0.1, 3, 30)

np.save(current_plotting_dir+'/num_paths', M)
np.save(current_plotting_dir+'/interpolation_points', N)
likelihood_evaluater(theta_array=theta_array,\
    alpha_array=alpha_array,\
    this_model=this_model,path_dir=current_plotting_dir )
