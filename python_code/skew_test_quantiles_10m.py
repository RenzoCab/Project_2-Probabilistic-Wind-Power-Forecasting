# %%
import os
import sys
os.chdir(sys.path[0])
from Base_plus import *
from scipy.stats import moment,skew,probplot, norm,shapiro, normaltest

forecast_data_in=np.load('./data/cleansed/URG_forecast_data_A_2018.npy')
#forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)

def take_center(data,fctr):
    upper_limit=np.mean(data)+ fctr*np.std(data)
    lower_limit=np.mean(data)- fctr*np.std(data)
    return(data[(data>lower_limit) & (data<upper_limit) ])

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
to_Lamperti=lambda data: np.arcsin(2*data-1)
from_Lamperti=lambda data: 0.5*(1+np.sin(data) )

# %%

# %%
X=forecast_data_inter[1,:-240,:] #fix this later !!! empty paths
P=forecast_data_inter[2,:-240,:]
V=X-P
dX=np.diff(X)
dV=np.diff(V)
dX_Lamprt = np.diff(to_Lamperti(X))
dV_Lamprt = np.diff(to_Lamperti(V+P))
# %%

np.max(X)
np.min(X)
np.max(V)
np.min(V)



plt.plot(V.flatten())



R=V.flatten()
R.shape
R[np.isfinite(R)].shape
np.diff(R).shape
plt.hist(R,20,density=True)
plt.xticks( fontsize = 15);
plt.yticks( fontsize = 20);
plt.xlabel('Forecast Error',fontsize = 20)
# plt.ylabel('Number of Incident',fontsize = 20)
# plt.legend( prop={'size': 20})
plt.title('Full Power Range' ,fontsize=24)
plt.savefig('hist_full.pdf',bbox_inches="tight")


probplot(R, plot=plt)




R[P.flatten()<0.3].shape
plt.hist(R[P.flatten()<0.3],20,density=True)
plt.xticks( fontsize = 15);
plt.yticks( fontsize = 20);
plt.xlabel('Forecast Error',fontsize = 20)
# plt.ylabel('Number of Incident',fontsize = 20)
# plt.legend( prop={'size': 20})
plt.title('Low Power Range' ,fontsize=24)
plt.savefig('hist_low.pdf',bbox_inches="tight")





probplot(R[P.flatten()<0.3], plot=plt)


R[P.flatten()>0.6].shape
plt.hist(R[P.flatten()>0.6],20,density=True)
plt.xticks( fontsize = 15);
plt.yticks( fontsize = 20);
plt.xlabel('Forecast Error',fontsize = 20)
# plt.ylabel('Number of Incident',fontsize = 20)
# plt.legend( prop={'size': 20})
plt.title('High Power Range' ,fontsize=24)
plt.savefig('hist_low.pdf',bbox_inches="tight")




probplot(R[P.flatten()>0.6], plot=plt)

R[(P.flatten()>=0.4)   & (P.flatten()<=0.6)].shape
plt.hist(R[ (P.flatten()>=0.4)   & (P.flatten()<=0.6)],20,density=True)

probplot(R[(P.flatten()>=0.4)   & (P.flatten()<=0.6)], plot=plt)


# AFTER Lamperti

Q=np.arcsin(2*(V+P)-1).flatten()
Q.shape
Q[np.isfinite(Q)].shape
np.diff(Q).shape
plt.hist(Q,20,density=True)

probplot(Q, plot=plt)


Q[P.flatten()<0.35].shape
plt.hist(Q[P.flatten()<0.35],20,density=True)
probplot(Q[P.flatten()<0.3], plot=plt)


Q[P.flatten()>0.35].shape
plt.hist(Q[P.flatten()>0.35],20,density=True)
probplot(Q[P.flatten()>0.6], plot=plt)

Q[(P.flatten()>=0.4)   & (P.flatten()<=0.6)].shape
plt.hist(Q[ (P.flatten()>=0.4)   & (P.flatten()<=0.6)],20,density=True)

probplot(Q[(P.flatten()>=0.4)   & (P.flatten()<=0.6)], plot=plt)

################################################################

# %%
fctr_in=0.5
dX_center=take_center(data=dX,fctr=fctr_in )
dV_center=take_center(data=dV,fctr=fctr_in )
dX_lamprt_center=take_center(data=dX_Lamprt,fctr=fctr_in )
dV_lamprt_center=take_center(data=dV_Lamprt,fctr=fctr_in )
# %%



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

plt.hist(np.diff(V[P>=0.35]), 30, density=True)
plot_QQ([ np.diff(V[ P<0.1])], ['Low Power'])

plot_QQ([ np.diff(V[ P<0.1])], ['Low Power'])




plt.hist(V[P>0.95])

plot_QQ([V[ P>0.95]], ['High Power'])


plot_QQ([ V[(0.3<P) & (P<0.6)]] , ['center'])

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

n=3000

def hill(sample,k):
    sample = sorted(np.array(sample),reverse=True)
    ratio = np.log(sample[:k-2])-np.log(sample[k-1])
    return np.mean(ratio), np.var(ratio)



a, m = 3, 1.  # shape and mode
k=40

s=m/(np.random.uniform(0,1,n*n)**(1./a))

type(s)

est1,var1= hill(s,k)

print(1/a)
print(est1-1./a)
print(var1)
print(np.sqrt((est1-1./a)**2+var1))



est1,var1=hill(R,k)
print(est1)
print(var1)

# %%
est=[]
var=[]
k=1000
for i in np.linspace(2000,60000, 10):
    sample=R[R>0]
    est_temp, var_temp=hill(sample[:int(i)],k);
    est.append(est_temp)
    var.append(var_temp)

# %%

plt.plot(est, '-o')

plt.plot(var, '-o')



est=[]
for i in np.linspace(100,60000, 10):
    sample=random.sample(R.tolist(), int(i))
    sample=np.asanyarray(sample)
    est.append(1/np.mean(np.diff(np.log(sorted(sample[sample>0], reverse=True)))))

plt.plot(est)







### until here

 dX.shape

dX_center=dX[(dX>lower_dX) & (dX<upper_dX) ]
dX_center.shape
dX.flatten().shape


upper_dV=np.mean(dV)+ fctr*np.std(dV)
lower_dV=np.mean(dV)- fctr*np.std(dV)



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
