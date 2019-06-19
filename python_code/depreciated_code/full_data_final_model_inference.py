from multi_path_base import *
forecast_with_data=np.load('forecast_with_data.npy')

###############################
#cleaning data
interpolation_points=73*2
forecast_data_inter= np.zeros((3,forecast_with_data.shape[0], interpolation_points ))
x = np.arange(1,72+1,1)
xnew = np.linspace(0,72,interpolation_points)
counter=0
i=0;j=0
while j< forecast_with_data.shape[0]: #fix redundant indicies

    y1=forecast_with_data[j,:,1]
    y2=forecast_with_data[j,:,2]
    inter_1=interpolate.interp1d(x, y1, fill_value='extrapolate')
    inter_2=interpolate.interp1d(x, y2, fill_value='extrapolate')
    p=inter_1(xnew)
    d=inter_2(xnew)

    if  sum(1*np.isnan(p))==0 and \
        sum(1*np.isnan(d))==0 and \
        np.max(p)<=1 and \
        np.min(p)>=0 and \
        np.max(d)<=1 and \
        np.min(d)>=0:

        forecast_data_inter[0,i,:]= p
        forecast_data_inter[1,i,:]= d
        forecast_data_inter[2,i,0]= forecast_with_data[j,0,0]
        i+=1
    j+=1
n_paths=i
print('paths left', i)
################################
N=forecast_data_inter.shape[2]
#M=forecast_data_inter.shape[1]
M= n_paths
dt=1

#forecast_data_inter[0,:,:] #forecast

disct_temp = disct(N,dt,M)
p=forecast_data_inter[0,:M,:]
V= forecast_data_inter[1,:M,:]-forecast_data_inter[0,:M,:]

#answer 9.06 , 0.473
this_model=model_modified_drift(disct_temp,V, forecast= p)
##################################
optim=this_model.optimize(np.array((9.5,0.5)))
optim
