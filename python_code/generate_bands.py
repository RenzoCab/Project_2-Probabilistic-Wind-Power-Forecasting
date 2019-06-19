import os
os.chdir('/Users/alhaddwt/Waleed Dropbox/Waleed Al-Haddad/GitLab/wind-power/python_code')
from Base import *
forecast_with_data=np.load('data/forecast_with_data.npy')
data=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='gen_CI_'
current_plotting_dir='plots/'+ script_name +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots


current_list=[31,437,82,623,75,59,68] #=list(range(11, data.shape[1])) #full_set
disct_in= disct(N=72, dt=1, M=1000) #M is number of paths wanted in ensamble
real_in=real(8,1) # SDE parameters to be generated from

empirical_Confidence_Interval_plots(forecast_data_inter=data,  real_in=real_in,\
            disct_in=disct_in, list_forecast_number=current_list,\
            dir_path=current_plotting_dir)
