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

N=72 #forecast_data_inter.shape[2]
M=forecast_data_inter.shape[1]
dt=1
#forecast_data_inter[0,:,:] #forecast

disct_temp = disct(N,dt,M)
p=forecast_data_inter[0,:M,:]
V= forecast_data_inter[1,:M,:]-forecast_data_inter[0,:M,:]
X=forecast_data_inter[1,:M,:]
Z = np.arcsin(2*X - 1)

#answer 9.06 , 0.473
this_model=model_modified_drift(disct_temp,Z, forecast= p)
#intial_point=np.array((1,0.5))
optim=this_model.optimize()

file_object  = open(current_data_dir+'/results.out', 'w')
note='Fit for ' + str(M) +' paths' + ' of ' + str(N) +' hours. Initialized at '\
    + str(intial_point)
file_object.write( note +'\n'+str(optim.x) +'\n' +optim.message[0])
file_object.close()

optim


optim

#Additions, brackets computations, lamparti
#
# co=np.zeros((forecast_data_inter.shape[1]))
# m3X=np.zeros((forecast_data_inter.shape[1]))
# m3Y=np.zeros((forecast_data_inter.shape[1]))
# skew_Y=np.zeros((forecast_data_inter.shape[1]))
# skew_X=np.zeros((forecast_data_inter.shape[1]))
#
#
# for sel in range(0,forecast_data_inter.shape[1]):
#     LHS=np.sum(np.diff(forecast_data_inter[1,sel,:])**2)
#     RHS=2*np.sum( forecast_data_inter[1,sel,:]* (1-forecast_data_inter[1,sel,:]) )
#     V_sel= forecast_data_inter[0,sel,:]-forecast_data_inter[1,sel,:]
#     m3X[sel]=np.power(V_sel, 3).sum() / 72
#     Y= np.arcsin(2*V_sel-1)
#     m3Y[sel]=np.power(Y, 3).sum() / 72
#     co[sel]=LHS/RHS
#     skew_Y[sel]=scipy.stats.skew(Y)
#     skew_X[sel]=scipy.stats.skew(V_sel)
#
#
# sel=500
# plt.plot(forecast_data_inter[0,sel,:])
# plt.plot(forecast_data_inter[1,sel,:], 'b')
#
# plt.plot(forecast_data_inter[0,sel,:]-forecast_data_inter[1,sel,:] )
#
#
# plt.plot(co)

# plt.plot(m3X,'b')
# plt.plot(m3Y, 'r')
#
# plt.plot(skew_X,'k')
#
# plt.plot(skew_Y,'b')
