import os
os.chdir('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code')
from Base import *
forecast_with_data=np.load('data/forecast_with_data.npy')
data=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='sample_paths_'
current_plotting_dir='plots/'+ script_name +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

print(current_plotting_dir)

current_list= list(range(0, data.shape[1])) #=list(range(0, data.shape[1])) #full_set
M=5
N=427
disct_in= disct(N=N, dt=1, M=M) #M is number of paths wanted in ensamble
real_in=real(22.32967705, 0.049445) # SDE parameters to be generated from

path_simulator(forecast_data_inter=data,hours=12,\
    real_in=real_in, disct_in=disct_in,list_forecast_number=current_list,\
    dir_path=current_plotting_dir)

path_simulator(forecast_data_inter=data,hours=24,\
    real_in=real_in, disct_in=disct_in,list_forecast_number=current_list,\
    dir_path=current_plotting_dir)
