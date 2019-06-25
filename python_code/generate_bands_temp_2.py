import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base import *
forecast_with_data=np.load('data/forecast_with_data.npy')
data=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='gen_CI_'
current_plotting_dir='plots/'+ script_name +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

print(current_plotting_dir)

current_list=list(range(0, data.shape[1])) #=list(range(0, data.shape[1])) #full_set
disct_in= disct(N=7, dt=1, M=1000) #M is number of paths wanted in ensamble
real_in=real(2.02,1.25) # SDE parameters to be generated from

empirical_Confidence_Interval_plots(forecast_data_inter=data,  real_in=real_in,\
            disct_in=disct_in, list_forecast_number=current_list,\
            dir_path=current_plotting_dir)

print(current_plotting_dir)
