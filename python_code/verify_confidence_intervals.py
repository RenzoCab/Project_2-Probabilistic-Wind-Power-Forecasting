import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base import *
forecast_with_data=np.load('data/forecast_with_data.npy')
forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)

chosen_folder='plots/gen_CI_19-06-21-14-44-39'


list_forecast_number=np.load(chosen_folder + '/forecast_list.npy')
N=np.load(chosen_folder + '/N.npy')
theta=np.load(chosen_folder + '/theta_0.npy')
alpha=np.load(chosen_folder + '/alpha.npy')
q975=np.load(chosen_folder + '/q975.npy')
q025=np.load(chosen_folder + '/q025.npy')
q95=np.load(chosen_folder + '/q95.npy')
q05=np.load(chosen_folder + '/q05.npy')
q50=np.load(chosen_folder + '/q50.npy')
q75=np.load(chosen_folder + '/q75.npy')
q25=np.load(chosen_folder + '/q25.npy')


bool_list=[]
for k in range(0,len(list_forecast_number)):
    bool_list.append(np.prod(q975[k,1:] >= forecast_data_inter[1,k,:]))
sum(bool_list)/len(list_forecast_number)


point_wise_upper=[]
point_wise_lower=[]
point_wise_total=[]
for k in range(1,N):
    percentage_upper=sum(q95[:,k] >= forecast_data_inter[1,:,k]*1)/len(list_forecast_number)
    percentage_lower=sum(q05[:,k] <= forecast_data_inter[1,:,k]*1)/len(list_forecast_number)
    percentage_total=sum((q95[:,k] >= forecast_data_inter[1,:,k]*1)*(q05[:,k] <= forecast_data_inter[1,:,k]*1))/len(list_forecast_number)
    point_wise_upper.append(percentage_upper)
    point_wise_lower.append(percentage_lower)
    point_wise_total.append(percentage_total)

plt.plot(point_wise_upper)

plt.plot(point_wise_lower)


plt.plot(point_wise_total)



path_point_wise_upper=[]
path_point_wise_lower=[]
path_point_wise_total=[]
q95[1,1:].shape
forecast_data_inter[1,1,:].shape
for k in range(0,len(list_forecast_number)):
    percentage_upper=sum(q95[k,1:] >= forecast_data_inter[1,k,:]*1)/N
    percentage_lower=sum(q05[k,1:] <= forecast_data_inter[1,k,:]*1)/N
    percentage_total=sum((q95[k,1:] >= forecast_data_inter[1,k,:]*1)*(q05[k,1:] <= forecast_data_inter[1,k,:]*1))/N
    path_point_wise_upper.append(percentage_upper)
    path_point_wise_lower.append(percentage_lower)
    path_point_wise_total.append(percentage_total)


plt.plot(path_point_wise_total,'o')

# calculation of the sum of prob assuming independence which we don't have
(9/10)*72

plt.plot(path_point_wise_upper,'o')

plt.plot(path_point_wise_lower,'o')


path_point_wise_total_array=np.asarray(path_point_wise_total)
np.quantile(path_point_wise_total_array,0.9)
