# %%
from multi_path_base import *
#%matplotlib tk
import os
import datetime as dtM
os.chdir('/Users/alhaddwt/Waleed Dropbox/Waleed Al-Haddad/GitLab/wind-power/python_code')
forecast_with_data=np.load('forecast_with_data.npy')
cwd = os.getcwd()
os.chdir(cwd+'/CI_FIG')
print('path changed to ' + cwd+'/CI_FIG')
# %%

#Further cleaning the data
interpolation_points=72
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

# %%
np.isnan(forecast_data_inter).any() #check that there are non nans left
# %%
(np.max(forecast_data_inter[0,:,:])>1).any()
# %%
(np.max(forecast_data_inter[1,:,:])>1).any()
# %%
(np.min(forecast_data_inter[0,:,:])<0).any()
# %%
(np.min(forecast_data_inter[1,:,:])<0).any()
# %%
def theta_adjust(theta, forecast):
    p=forecast
    N=len(p)
    p_dot=np.diff(p)*N

    zero_drift = np.diff(p)*N/theta + p[:-1]
    zero_drift_fixed=np.zeros_like(zero_drift)
    theta_adjusted=np.zeros_like(zero_drift)

    for i in range(0,len(zero_drift)):
        zero_drift_fixed[i] = p_dot[i]/max(theta, abs(p_dot[i])/min(p[i], 1-p[i]) ) + p[i]
        theta_adjusted[i]= max(theta, abs(p_dot[i])/min(p[i], 1-p[i]) )

    if max(zero_drift_fixed)>1:
            print('WARNING: outside range ! bigger than one ')
    if min(zero_drift_fixed)<0 and abs(min(zero_drift_fixed)) > 1e-10:
            print('WARNING: outside range ! less than zero')

    #plt.figure()
    #plt.plot(p)

    #plt.figure()
    #plt.plot(p_dot)

    #plt.figure()
    #plt.plot(zero_drift)

    #plt.figure()
    #plt.plot(zero_drift_fixed)

    #plt.figure()
    #plt.plot(theta_adjusted)

    return(theta_adjusted, zero_drift_fixed )


# %%
p=forecast_data_inter[0,3,:]
theta_adjusted, zero_drift_fixed=theta_adjust(12.3,p)
plt.plot(zero_drift_fixed)
plt.figure()
plt.plot(theta_adjusted)
# %%
def gen_X_normal_euler_DT_modified(X0,disct,real,forecast):
    p=forecast
    N=disct.N
    M=disct.M
    X0=p[0]
    theta=real.mu
    alpha=real.sigma
    X=np.zeros((M,N))
    #X_zero_drift=np.zeros(N)
    dN=1/N
    i=0
    j=0
    while j < M:
        X[j,0]=X0
        a=0
        b=0
        dW=np.sqrt(dN)*np.random.normal(0, 1, N)
        i=0
        while i<N-1:

            b=  2*alpha*theta[i]*p[i]*(1-p[i])*X[j,i]*(1-X[j,i])
            a= - theta[i]*(X[j,i]-p[i]) + (p[i+1]-p[i])/(dN)
            if b<0:
                print('negative b')

            X[j,i+1]= X[j,i] + a*dN+ np.sqrt(b)*dW[i]

            if X[j,i+1] >1:
                X[j,i+1]=X[j,i];
            if X[j,i+1] <0:
                X[j,i+1]=X[j,i];
            #X_zero_drift[i]= p[i]+ (p[i+1]-p[i])/(dN*theta[i]);

            i+=1
        j+=1
    return( X )
# %%
for k in range(1064,1067): #n_paths
    p=forecast_data_inter[0,k,:] #obtain cleansed forecast
    d=forecast_data_inter[1,k,:]
    dt=1
    M_test=200
    N=interpolation_points

    theta_adjusted, zero_drift_fixed=theta_adjust(theta,p)

    #theta_const= np.zeros_like(theta_array)+12.3

    real_1 = real(theta_adjusted, alpha ) #theta_array

    disct_temp = disct(N,dt,M_test)

    X = np.empty((M_test,N));
    X= gen_X_normal_euler_DT_modified(X0=p[0],disct=disct_temp,real=real_1,forecast=p)
    plt.plot(X[0,:])
# %%
#%%capture

N=interpolation_points

q975=np.zeros((N))
q025=np.zeros((N))
q95=np.zeros((N))
q05=np.zeros((N))
q50=np.zeros((N))
q75=np.zeros((N))
q25=np.zeros((N))

xnew = np.linspace(0,72,interpolation_points)


theta= 7.8155
alpha= 1.072

for k in [31,437,82,623,75,59,68]:#range(0,n_paths): #n_paths
    p=forecast_data_inter[0,k,:] #obtain cleansed forecast
    d=forecast_data_inter[1,k,:]
    dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[2,k,0])
    dt=1
    M_test=1000

    fig=plt.figure(2,figsize=(10, 4))
    fig.clf()
    plt.plot(xnew,p, label='forecast')
    plt.plot(xnew,d, label='actual production')

    plt.xlim(1, 73)
    plt.ylim(-0.1, 1.1)
    plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

    plt.xticks( fontsize = 20);
    plt.yticks( fontsize = 20);
    plt.xlabel('Time [hr]',fontsize = 24)
    plt.ylabel('Power',fontsize = 24)
    plt.legend( prop={'size': 15})
    plt.savefig('Forecast_data_'+ str(k)+'.pdf', bbox_inches="tight")


    theta_adjusted, zero_drift_fixed=theta_adjust(theta,p)

    #theta_const= np.zeros_like(theta_array)+12.3

    real_1 = real(theta_adjusted, alpha ) #theta_array

    disct_temp = disct(N,dt,M_test)

    X = np.empty((M_test,N));
    X= gen_X_normal_euler_DT_modified(X0=p[0],disct=disct_temp,real=real_1,forecast=p)
    #plt.plot(X[0,:])

    for i in range(0,N):
        q975[i]=np.quantile(X[:,i], 0.975)
        q025[i]=np.quantile(X[:,i], 0.025)

        q95[i]=np.quantile(X[:,i], 0.95)
        q05[i]=np.quantile(X[:,i], 0.05)

        q50[i]=np.quantile(X[:,i], 0.5)

        q75[i]=np.quantile(X[:,i], 0.75)
        q25[i]=np.quantile(X[:,i], 0.25)

    fig=plt.figure(2,figsize=(10, 4))
    fig.clf()

    plt.xlim(1, 73)
    plt.ylim(-0.1, 1.1)

    plt.fill_between(xnew, q75,q25,color='k',alpha=0.2, label='90% CI', edgecolor=None)
    plt.fill_between(xnew, q95,q05,color='k',alpha=0.3, label='50% CI', edgecolor=None)

    #plt.plot(xnew,np.mean(normal_X_derivative_tracking_Euler, axis=0),'c-', label='Mean')

    #plt.plot(xnew,q50,'y-', label='Median')
    #dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[2,k,0])
    plt.plot(xnew,p, 'r-', label='forecast',linewidth=3)
    #plt.plot(xnew[:-1],zero_drift_fixed, 'y-', label='Zero Drift Line',linewidth=1)
    plt.plot(xnew,d , 'y-', label='actual production',linewidth=3)
    plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

    plt.xticks( fontsize = 20);
    plt.yticks( fontsize = 20);
    plt.xlabel('Time [hr]',fontsize = 24)
    plt.ylabel('Power',fontsize = 24)
    plt.legend( prop={'size': 15})
    plt.savefig('72hr_forecast_CI_'+ str(k)+'.pdf', bbox_inches="tight")

    fig= plt.figure(2,figsize=(10, 4))
    fig.clf()

    plt.xlim(1, 7)
    plt.ylim(-0.1, 1.1)

    plt.fill_between(xnew, q75,q25,color='k',alpha=0.2, label='90% CI', edgecolor=None)
    plt.fill_between(xnew, q95,q05,color='k',alpha=0.3, label='50% CI', edgecolor=None)

    #plt.plot(xnew,np.mean(normal_X_derivative_tracking_Euler, axis=0),'c-', label='Mean')

    #plt.plot(xnew,q50,'y-', label='Median')

    plt.plot(xnew,p, 'r-', label='forecast',linewidth=3)
    plt.plot(xnew,d , 'y-', label='actual production',linewidth=3)

    #plt.plot(xnew[:-1],zero_drift_fixed, 'y-', label='Zero Drift Line',linewidth=1)

    plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24) #,fontsize=10 Forecast Confidence Intervals
    plt.xlabel('Time [hr]',fontsize = 24)
    plt.ylabel('Power',fontsize = 24)
    plt.legend( prop={'size': 15})
    plt.xlim(1, 7)
    plt.ylim(-0.1, 1.1)
    plt.xticks( fontsize = 20);
    plt.yticks( fontsize = 20);
    plt.savefig('6hr_forecast_CI_'+ str(k)+'.pdf', bbox_inches="tight")

# %%
