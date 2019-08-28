import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base_plus import *
#from multi_path_base import model_modified_drift
forecast_data_in=np.load('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code/data/cleansed/URG_forecast_data_A_2018.npy')
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='likelihood_explorer'
current_plotting_dir='data/'+ script_name+'_' +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)

N=forecast_data_inter.shape[2];N
M=forecast_data_inter.shape[1]-500;M
dt=1
dN=1/N;dN

disct_temp = disct(N,dt,M)
p=forecast_data_inter[2,:-500,:]  #240
V= forecast_data_inter[2,:-500,:] -forecast_data_inter[1,:-500,:]

#plt.plot(p[105])
#plt.plot(forecast_data_inter[1,105,:])
#plt.plot(V[105,:])
#X=forecast_data_inter[1,:-500,:]

#Z = np.arcsin(2*X - 1)

#now for V
this_model=model_modified_drift(disct_temp,V, forecast= p) #change it

theta_array=np.linspace(5, 20, 50) #np.linspace(1, 30, 30)
alpha_array=np.linspace(0.005, 0.15, 50) #np.linspace(0.1, 3, 30)

np.save(current_plotting_dir+'/num_paths', M)
np.save(current_plotting_dir+'/interpolation_points', N)
likelihood_evaluater(theta_array=theta_array,\
    alpha_array=alpha_array,\
    this_model=this_model,inference="beta_likelihood", path_dir=current_plotting_dir )
print(current_plotting_dir)
