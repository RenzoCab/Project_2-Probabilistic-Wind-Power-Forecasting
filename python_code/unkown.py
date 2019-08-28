# %%
from multi_path_base import *
import pandas as pd
import datetime as dt
from scipy import interpolate
%matplotlib tk
# %%
f = pd.read_csv('sal_utep3_h72_2016_130317.txt', delimiter= '\s+', index_col=False,header=1)
d = pd.read_csv('sal_uter1_2016_130317_iniper.txt', delimiter= '\s+', index_col=False,header=1)
n = pd.read_csv('historico_potnom_uruy.txt', delimiter= '\s+', index_col=False,header=None)

# %%
#setting time and date of the forecasts
p=f.values
length_forecast=72 #[hr]
number_forecasts=int(len(p[:,0])/length_forecast)
#seperate the forecasts
forecasts=p.reshape((number_forecasts,72, 10))

#change format of date and time to a datetime object
forecast_clean=np.empty((forecasts.shape[0],forecasts.shape[1],2))
for i in range(0,number_forecasts):
    for j in range(0,72):
        date_object= dt.datetime.strptime( forecasts[i,j,1] , "%d/%m/%y").date()
        time_object= dt.datetime.strptime(format(forecasts[i,j,6], '04d'), "%H%M").time()
        datetime_object = dt.datetime.combine(date_object, time_object)
        forecast_clean[i,j,0]=datetime_object.timestamp()
        forecast_clean[i,j,1]=forecasts[i,j,8]

#forecasts[1,1,1].timestamp()
#forecasts[1,1,1].isoformat()
#dt.datetime.fromtimestamp(timestamp)
# %%
#setting time and data of the data
data=d.values
data_clean=np.zeros((2,data.shape[0]))
for i in range(0,data.shape[0]):
    date_object= dt.datetime.strptime( data[i,2] , "%d/%m/%Y").date()
    time_object= dt.datetime.strptime(data[i,3], "%H:%M").time()
    datetime_object = dt.datetime.combine(date_object, time_object)

    data_clean[0,i]=datetime_object.timestamp()
    data_clean[1,i]=data[i,4]
#plt.plot(data_clean[1,:])
# %%
#cleaning normalization data

normalization_data=n.values
normalization_clean=np.zeros((3,normalization_data.shape[0]))

for i in range(0,normalization_data.shape[0]):
    date_object= dt.datetime.strptime( normalization_data[i,0] , "%d/%m/%Y").date()
    time_object= dt.datetime.strptime(normalization_data[i,1], "%H:%M").time()
    datetime_object_st = dt.datetime.combine(date_object, time_object)
    normalization_clean[0,i]=datetime_object_st.timestamp()

    date_object= dt.datetime.strptime( normalization_data[i,2] , "%d/%m/%Y").date()
    time_object= dt.datetime.strptime(normalization_data[i,3], "%H:%M").time()
    datetime_object_st = dt.datetime.combine(date_object, time_object)
    normalization_clean[1,i]=datetime_object_st.timestamp()

    normalization_clean[2,i]= float(normalization_data[i,4].replace(',','.'))

# %%
#matching forecast with data
normalization=float('nan')
forecast_with_data= np.empty((forecast_clean.shape[0],\
                              forecast_clean.shape[1],forecast_clean.shape[2]+1))
for j in range(0,number_forecasts):
    for i in range(0, forecast_clean[j,:,0].size ):
        temp=forecast_clean[j,i,0]
        index=np.where(data_clean[0,:] == temp )
        forecast_with_data[j,i,0]= temp
        forecast_with_data[j,i,1]= forecast_clean[j,i,1]
        if data_clean[1, index[0]].size != 0:
            #print(dt.datetime.fromtimestamp(temp))
            #print(data_clean[1, index[0]].size)
            for k in range(0, normalization_clean.shape[1]):
                if (normalization_clean[0,k] <= temp) and ( temp <= normalization_clean[1,k]):
                    normalization= normalization_clean[2,k]
            #print(dt.datetime.fromtimestamp(temp))
            #print(normalization)
            forecast_with_data[j,i,2]= data_clean[1, index[0]][0]/normalization
            forecast_with_data[j,i,1]=forecast_with_data[j,i,1]/(normalization*1e3)
        else: forecast_with_data[j,i,2]= float('nan')
#index[2][:]=1
# %%
#save output
np.save('forecast_with_data', forecast_with_data)
# %%
fig = plt.figure(1)
fig.clf()
for i in range(0,4): #number_forecasts
    ax = fig.add_subplot(2,2,i+1)
    plt.xlim(0, 72)
    plt.ylim(0, 1)
    plt.plot(forecast_with_data[i,:,1], label='forecast')
    plt.plot(forecast_with_data[i,:,2], label='data')
plt.legend()
# %%
fig = plt.figure(1)
fig.clf()
for i in range(0,100): #number_forecasts
    ax = fig.add_subplot(10,10,i+1)
    plt.xlim(0, 72)
    plt.ylim(0, 1)

    #ind=random.randint(1,30)#forecast_with_data.shape[0])

    plt.plot(forecast_with_data[i+2,:,1], label='forecast')
    plt.plot(forecast_with_data[i+2,:,2], label='data')
    plt.xlabel('Time [hr]')
    plt.ylabel('Power')

plt.legend()
plt.savefig('forecast_data.eps')
# %%
~np.isnan(forecast_with_data)
# %%
forecast_with_data[~np.isnan(forecast_with_data).any(axis=2)]
# %%
np.where(forecast_with_data[:,:,1])
# %%
dt=0.001
x = np.arange(1,72+1,1) #forecast_with_data[5,:,0]
fig = plt.figure(2)
fig.clf()
#plt.title('Cubic-spline interpolation of the Forecast')

for i in range(0,4):

    ind=random.randint(1,30)

    y = forecast_with_data[ind,:,1]
    tck = interpolate.splrep(x, y, s=0)
    xnew = np.arange(1,72+1,dt)
    ynew = interpolate.splev(xnew, tck, der=0)

    ax = fig.add_subplot(2,2,i+1)
    plt.xlim(0, 72)
    plt.ylim(0, 1)
    #plt.axis('off')
    plt.grid(True)
    plt.plot(x,y,'x', label='Discrete' )
    plt.plot(xnew,ynew,'-', label='cubic splines')
    plt.xlabel('Time [hr]')
    plt.ylabel('Forecast')

plt.legend()
plt.show()
plt.savefig('forecast_splines.eps')
# %%
dt=0.001
x = np.arange(1,72+1,1) #forecast_with_data[5,:,0]
fig = plt.figure(2)
fig.clf()
for i in range(0,4):

    ind=random.randint(1,30)

    y = forecast_with_data[ind,:,1]
    tck = interpolate.splrep(x, y, s=0)
    xnew = np.arange(1,72+1,dt)
    ynew = interpolate.splev(xnew, tck, der=0)

    ax = fig.add_subplot(7,7,i+1)
    plt.xlim(0, 72)
    plt.ylim(0, 1)
    #plt.axis('off')
    plt.grid(True)
    plt.plot(x,y,'x')
    plt.plot(xnew,ynew,'-')
    #plt.title('Cubic-spline interpolation of the Forecast')
    plt.show()
# %%
dt=0.001
x = np.arange(1,72+1,1) #forecast_with_data[5,:,0]

y = forecast_with_data[5,:,1]
tck = interpolate.splrep(x, y, s=0.5)
xnew = np.arange(1,72+1,dt)
ynew = interpolate.splev(xnew, tck, der=0)

#ax = fig.add_subplot(7,7,i+1)
plt.figure()
plt.xlim(0, 72)
plt.ylim(0, 1)
plt.plot(x,y,'o')
plt.grid(True)
plt.plot(xnew,ynew,'-')
plt.title('Cubic-spline interpolation of the Forecast')
# %%
np.arange(1,72+1,0.01)
# %%
forecast_with_data.shape
# %%
max(2,3)
# %%
