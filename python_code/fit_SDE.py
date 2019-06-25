import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base import *
forecast_with_data=np.load('data/forecast_with_data.npy')
forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='fit_SDE'
current_data_dir='data/'+ script_name +current_time
os.mkdir(current_data_dir)

print(current_data_dir)

N=7 #forecast_data_inter.shape[2]
M=forecast_data_inter.shape[1]
dt=1
#forecast_data_inter[0,:,:] #forecast

disct_temp = disct(N,dt,M)
p=forecast_data_inter[0,:M,:]
V= forecast_data_inter[1,:M,:]-forecast_data_inter[0,:M,:]

#answer 9.06 , 0.473
this_model=model_modified_drift(disct_temp,V, forecast= p)
intial_point=np.array((8.5,0.5))
optim=this_model.optimize(intial_point)

file_object  = open(current_data_dir+'/results.out', 'w')
note='Fit for ' + str(M) +' paths' + ' of ' + str(N) +' hours. Initialized at '\
    + str(intial_point)
file_object.write( note +'\n'+str(optim.x) +'\n' +optim.message[0])
file_object.close()

optim
