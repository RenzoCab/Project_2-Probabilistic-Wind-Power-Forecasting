# %%
import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base import *
from scipy.stats import moment,skew,probplot, norm,shapiro, normaltest
forecast_with_data=np.load('data/forecast_with_data.npy')
forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)

def plot_QQ(list_data,labels, colors=['r','k','y', 'm', 'g']):
    plt.clf()
    fig = plt.figure()
    err=[]
    i=0
    for i in range(0,len(list_data)):
        x,y=probplot(list_data[i].flatten(), dist="norm")[0]
        plt.plot(x,y, 'o'+ colors[i], label=labels[i] )
        plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), colors[i])
        plt.xlabel('Theoretical quantiles', fontsize=18)
        plt.ylabel('Ordered Values', fontsize=18)
        plt.legend()
        err.append( np.mean( (np.poly1d(np.polyfit(x, y, 1))(np.unique(x)) - y)**2 ) )
    return(err)
dX.flatten().shape[0]

N=np.random.normal(size=dX.flatten().shape[0])


forecast_data_inter.shape







shapiro(N)

shapiro(dX.flatten())
shapiro(dX_Lamprt.flatten())

shapiro(dX_center.flatten())
shapiro(dX_lamprt_center.flatten())

skew(dX.flatten())
skew(dX_Lamprt.flatten())

skew(dX_center.flatten())
skew(dX_lamprt_center.flatten())


skew(dX_center.flatten())
skew(dX_lamprt_center.flatten())

def take_center(data,fctr):
    upper_limit=np.mean(data)+ fctr*np.std(data)
    lower_limit=np.mean(data)- fctr*np.std(data)
    return(data[(data>lower_limit) & (data<upper_limit) ])

# %%
fctr_in=0.5
dX_center=take_center(data=dX,fctr=fctr_in )
dV_center=take_center(data=dV,fctr=fctr_in )
dX_lamprt_center=take_center(data=dX_Lamprt,fctr=fctr_in )
dV_lamprt_center=take_center(data=dV_Lamprt,fctr=fctr_in )

# %%
V[P<0.1].shape

np.mean(np.diff(V[P<0.3]))
np.mean(np.diff(V[P<0.4]))
np.mean(np.diff(V[P<0.5]))
np.mean(np.diff(V[P<0.6]))

np.median(np.diff(V[P<0.1]))
np.median(np.diff(V[P<0.3]))
np.median(np.diff(V[P<0.4]))
np.median(np.diff(V[P<0.5]))
np.median(np.diff(V[P<0.6]))






plt.hist(np.diff(V.flatten())[P.flatten()<0.35],30,density=True)

np.diff(V[P>0.35]).shape
np.diff(V[P<=0.35]).shape

plt.hist(np.diff(V[P>=0.35]), 30)
plot_QQ([ np.diff(V[ P<0.1])], ['Low Power'])

plot_QQ([ np.diff(V[ P<0.1])], ['Low Power'])




plt.hist(V[P>0.95])

plot_QQ([V[ P>0.95]], ['High Power'])

plt.hist(V[ (0.3<P) & (P<0.6) ])

plot_QQ([ V[(0.3<P) & (P<0.6)] , V[(0.1<P) & (P<0.8)] ], ['center', 'large center'])

dX_Lamprt_p=np.arcsin(2*X-1)/np.sqrt(P*(1-P))


# Best Result take this until

Q=np.arcsin(2*(V+P)-1).flatten()
Q.shape
Q[np.isfinite(Q)].shape
np.diff(Q).shape
plt.hist(Q[P.flatten()>0.4],20,density=True)


probplot(R[P.flatten()>0.4], plot=plt)

R=V.flatten()
R.shape
R[np.isfinite(Q)].shape
np.diff(R).shape
plt.hist(R[P.flatten()<0.4],20,density=True)

probplot(R[P.flatten()<0.4], plot=plt)


R=V.flatten()
R.shape
R[np.isfinite(Q)].shape
np.diff(R).shape
plt.hist(R[P.flatten()>0.4],20,density=True)


R=V.flatten()
R.shape
R[np.isfinite(Q)].shape
np.diff(R).shape
plt.hist(R,20,density=True)

normaltest(R)
skew(R)
moment(R,4)/np.std(R)**4

normaltest(Q)
skew(Q)
moment(Q,4)/np.std(Q)**4

#filter # R:

normaltest(R[P.flatten()<0.4])
skew(R[P.flatten()<0.4])
moment(R[P.flatten()<0.4],4)/np.std(R[P.flatten()<0.4])**4


probplot(R, sparams=10 , dist='pareto', plot=plt )


# no power law hill estimator -> jumps
est=[]
for i in np.linspace(100,60000, 10):
    sample=random.sample(R.tolist(), int(i))
    sample=np.asanyarray(sample)
    est.append(1/np.mean(np.diff(np.log(np.sort(sample[sample>0])))))

plt.plot(est)

est=[]
for i in np.linspace(100,60000, 10):
    sample=random.sample(R.tolist(), int(i))
    sample=np.asanyarray(sample)
    est.append(1/np.mean(np.diff(np.log(np.sort(-1*sample[sample<0])))))

plt.plot(est)


### until here

 dX.shape

dX_center=dX[(dX>lower_dX) & (dX<upper_dX) ]
dX_center.shape
dX.flatten().shape


upper_dV=np.mean(dV)+ fctr*np.std(dV)
lower_dV=np.mean(dV)- fctr*np.std(dV)



to_Lamperti=lambda data: np.arcsin(2*data-1)/np.sqrt(P*(1-P))
from_Lamperti=lambda data: 0.5*(1+np.sin(data) )
# %%

# %%
X=forecast_data_inter[1,:,:]
P=forecast_data_inter[0,:,:]
V=X-P
dX=np.diff(X)
dV=np.diff(V)
dX_Lamprt = np.diff(to_Lamperti(X))
dV_Lamprt = np.diff(to_Lamperti(V+P))
# %%

diffp=np.diff(P[760,:])

plt.plot(diffp)


P[10,:]

P.shape



L=[X, V ]
labels=['X', 'V']
plot_QQ(L,labels=labels)




L=[dX, dV ]
labels=['dX', 'dV']
plot_QQ(L,labels=labels)


L=[dX,dX_Lamprt,dV , dV_Lamprt]
labels=['dX','dX Lamperti' , 'dV','dV Lamprti' ]
plot_QQ(L,labels=labels)

L=[dV , dV_Lamprt]
labels=['dV','dV Lamprti' ]
plot_QQ(L,labels=labels)


dX_Lamprt_p.flatten()


L=[dX_Lamprt_p.flatten()]
labels=['dX_Lamprt corrected' ]
err=plot_QQ(L,labels=labels)




L=[dX,dX_Lamprt, dX_Lamprt_p]
labels=['dX','dX Lamperti', 'dX_Lamprt corrected' ]
err=plot_QQ(L,labels=labels)

L=[dX_center,dX_lamprt_center]
labels=['dX','dX Lamperti' ]
err=plot_QQ(L,labels=labels)




plt.plot(err[0])
plt.plot(err[1])
#plt.xlim(-2, 2)



L=[dX]
labels=['dX']
plot_QQ(L,labels=labels)



plt.hist(np.diff(X.flatten()), 30)


L=[np.diff(X.flatten()) ]
labels=['dX']
plot_QQ(L,labels=labels)




plt.hist( np.diff(np.arcsin(2*X.flatten() -1)), 30 )

L=[np.diff(X.flatten()) ]
labels=['dX']
plot_QQ(L,labels=labels)


L=[np.diff(np.arcsin(2*X -1)).flatten() ]
labels=['dX Lamperti']
plot_QQ(L,labels=labels)
