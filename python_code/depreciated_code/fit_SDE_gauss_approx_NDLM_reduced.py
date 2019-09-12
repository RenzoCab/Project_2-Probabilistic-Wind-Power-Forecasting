import os
os.chdir('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code')
from Base_plus import *
forecast_data_in=np.load('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code/data/cleansed/URG_forecast_data_A_2018.npy')
#forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='fit_SDE'
current_data_dir='data/'+ script_name +current_time
os.mkdir(current_data_dir)

forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)
print(current_data_dir)



#check no nans
np.sum( np.isnan(forecast_data_inter.flatten()))
np.max(forecast_data_inter[1:,:,:])<=1
np.min(forecast_data_inter[1:,:,:])>=0

N=forecast_data_inter.shape[2]
M=forecast_data_inter.shape[1]-240 #to be changed in generalization
dt=1
#forecast_data_inter[0,:,:] #forecast

N
M
dt

disct_temp = disct(N,dt,M)
p=forecast_data_inter[2,:-240,:]
V= forecast_data_inter[2,:-240,:]-forecast_data_inter[1,:-240,:]
X=forecast_data_inter[1,:-240,:]
#Z = np.arcsin(2*X - 1)

#answer 9.06 , 0.473
this_model=model_modified_drift(disct_temp,V, forecast= p)
intial_point= np.array((1,0.8))
optim=this_model.optimize(inference="lamperti_likelihood_SDE_approx" ,method="Nelder-Mead",param_initial=intial_point)


optim


file_object  = open(current_data_dir+'/results.out', 'w')
note='Fit for ' + str(M) +' paths' + ' of ' + str(N) +' hours. Initialized at '\
    + str(intial_point)
file_object.write( note +'\n'+str(optim.x) +'\n' +optim.message[0])
file_object.close()

optim


# co=np.zeros((forecast_data_inter.shape[1]))
# for sel in range(0,forecast_data_inter.shape[1]):
#     LHS=np.sum(np.diff(forecast_data_inter[1,sel,:])**2)
#     RHS=2*np.sum( forecast_data_inter[1,sel,:]* (1-forecast_data_inter[1,sel,:])   )
#     co[sel]=LHS/RHS
#
# plt.plot(co)
#
# np.mean(co)
#
# forecast_data_inter[1,:M,:]
#using NDLMR remains to check initial value of first transition
########
