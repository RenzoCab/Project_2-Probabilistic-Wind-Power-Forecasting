
import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base_plus import *
#from multi_path_base import model_modified_drift
forecast_with_data=np.load('data/forecast_with_data.npy')
forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='likelihood_explorer'
current_plotting_dir='data/'+ script_name+'_' +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

N=forecast_data_inter.shape[2]-1
M= 100 #forecast_data_inter.shape[1]
dt=1
dN=1/N
dN
disct_temp = disct(N,dt,M)
p=forecast_data_inter[0,:,:]
V= forecast_data_inter[1,:,:]-forecast_data_inter[0,:,:]

N
dN

this_model=model_modified_drift(disct_temp,V, forecast= p)

theta_array=np.linspace(1, 100, 50)
alpha_array=np.linspace(0.001, 0.15, 50)

np.save(current_plotting_dir+'/num_paths', M)
np.save(current_plotting_dir+'/interpolation_points', N)


likelihood_evaluater(theta_array=theta_array,\
    alpha_array=alpha_array,\
    this_model=this_model,inference="beta_likelihood" , path_dir=current_plotting_dir )
print(current_plotting_dir)

#this_model.likelihood(np.array((10,0.11)))
#this_model.likelihood(np.array((60,0.1)))
