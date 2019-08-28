import os
os.chdir('/Users/alhaddwt/Google Drive/GitLab/wind-power/python_code')
from Base import *
from scipy.stats import moment,skew,probplot, norm
forecast_with_data=np.load('data/forecast_with_data.npy')
forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)


# %%
def plot_QQ(list_data,labels, colors=['r','k','y', 'm', 'g']):
    plt.clf()
    fig = plt.figure()
    for i in range(0,len(list_data)):
        x,y=probplot(list_data[i])[0]
        plt.plot(x,y, 'o'+ colors[i], label=labels[i] )
        plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), colors[i])
        plt.xlabel('Theoretical quantiles', fontsize=18)
        plt.ylabel('Ordered Values', fontsize=16)
        plt.legend()
    return()
# %%

X=forecast_data_inter[1,:,:]
P=forecast_data_inter[0,:,:]
V=forecast_data_inter[1,:,:] -forecast_data_inter[0,:,:]


Xm3=moment(V, axis=1, moment=3)
Vm3=moment(X, axis=1, moment=3)
Z=np.arcsin(2*X-1)
Zm3=moment(Z, axis=1, moment=3)


plt.hist(np.diff(V.flatten()), 30, density=True)

plt.hist(np.diff(Z).flatten() ,30, density=True, color='r')






mean,std=norm.fit(np.diff(X).flatten())
range_1=np.linspace(-1,1,100)
y=norm.pdf(range_1,mean,std)
plt.plot(range_1,y, 'k')
plt.hist(np.diff(X).flatten() ,30, density=True , color='k')

mean,std=norm.fit(np.diff(Z).flatten())
range_1=np.linspace(-1,1,100)
y=norm.pdf(range_1,mean,std)
plt.plot(range_1,y, 'r')
plt.hist(np.diff(Z).flatten() ,30, density=True, color='r')




















Xm4=moment(V, axis=1, moment=4)
Vm4=moment(X, axis=1, moment=4)

plt.plot(Xm4)

plt.plot(Vm4, 'r')


plt.plot(X[350,:])
plt.ylim(0, 1)
plt.hist(V[350,:], 10)

plt.hist(X[350,:], 10)

plt.hist(Z[350,:], 10)


Xm3[350]

plt.plot(Xm3)

plt.plot(Vm3)


plt.plot(Zm3)






plt.plot(skew(V, axis=1))

plt.plot(skew(np.diff(X), axis=1))
plt.plot(skew(np.diff(Z), axis=1), 'r')


plt.plot(X[350])
plt.ylim(0, 1)


plt.plot(np.diff(X[350]))
plt.ylim(0, 1)

(np.diff(X[350])**3).sum()
Xm3[350]

plt.plot(X[350])
plt.ylim(0, 1)
plt.hist(np.diff(X[350]),10)

x1,x2=probplot(np.diff(X[350]), dist="norm", plot=plt)


plt.plot(X[350])
plt.ylim(0, 1)

fig = plt.figure()
ax = fig.add_subplot(111)
y1,y2=probplot(np.diff(X[350]), dist="norm", plot=ax)
y1,y2=probplot(np.diff(Z[350]), dist="norm", plot=ax)
ax.set_title("-")
ax.get_lines()[0].set_marker('o')
ax.get_lines()[0].set_markerfacecolor('r')
ax.get_lines()[0].set_markersize(6.0)
plt.show()


probplot(np.diff(Z[350]), dist="norm", plot=plt)


fig = plt.figure()
ax = fig.add_subplot(111)
#_,_=probplot(np.diff(X.flatten()), dist="norm", plot=ax)
_,_=probplot(np.diff(Z.flatten()), dist="norm", plot=ax)
ax.set_title("-")
ax.get_lines()[0].set_marker('x')
ax.get_lines()[0].set_markerfacecolor('r')
ax.get_lines()[0].set_markersize(2.0)
ax.get_lines()[1].set_linewidth(8.0)
plt.show()

X.shape
P.shape
q=X + P

np.arcsin( 2*X -1)
np.arcsin( 2*(V + P) -1)
W=np.arcsin( 2*(V + P) -1)

probplot(np.diff(W.flatten()), dist="norm", plot=plt)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = fig.add_subplot(111)
_,_=probplot(np.diff(V.flatten()), dist="norm", plot=ax1)
_,_=probplot(np.diff(W.flatten()), dist="norm", plot=ax2)

ax1.get_lines()[0].set_marker('x')
ax1.get_lines()[0].set_markerfacecolor('r')
ax1.get_lines()[0].set_markersize(2.0)
ax1.get_lines()[0].set_linewidth(8.0)

ax2.get_lines()[0].set_marker('o')
ax2.get_lines()[0].set_markerfacecolor('k')
ax2.get_lines()[0].set_markersize(1.0)
ax2.get_lines()[0].set_linewidth(8.0)

plt.show()



fig = plt.figure()
ax = fig.add_subplot(111)
_,_=probplot(np.diff(X[350]), dist="norm", plot=ax)
_,_=probplot(np.diff(Z[350]), dist="norm", plot=ax)
ax.set_title("-")
ax.get_lines()[0].set_marker('x')
ax.get_lines()[0].set_markerfacecolor('r')
ax.get_lines()[0].set_markersize(2.0)
ax.get_lines()[1].set_linewidth(8.0)
plt.show()

plt.plot(y1[1],'o')

x = np.linspace(-5,5,100)
y =  y2[0]*x + y2[1]
plt.plot(x,y)
plt.plot(y1[1],'o')


y2
_,_=probplot(np.diff(Z[350]), dist="norm", plot=plt)

_,_=probplot(np.diff(X[350]), dist="norm", plot=plt)
plt.plot(X[190])
plt.ylim(0, 1)

_,_=probplot(np.diff(X[190]), dist="norm", plot=plt)

_,_=probplot(np.diff(V[190]), dist="norm", plot=plt)

_,_=probplot(np.diff(Z[190]), dist="norm", plot=plt)


plt.plot(X[604])
plt.ylim(0, 1)

_,_=probplot(np.diff(X[604]), dist="norm", plot=plt)

_,_=probplot(np.diff(V[604]), dist="norm", plot=plt)

_,_=probplot(np.diff(Z[604]), dist="norm", plot=plt)






plt.plot(*probplot(np.diff(Z[604]))[0], 'o')


probplot(np.diff(Z[604]))[0]
np.diff(Z[604])

# %%
def plot_QQ(list_data,labels, colors=['r','k','y', 'm', 'g']):
    plt.clf()
    fig = plt.figure()
    for i in range(0,len(list_data)):
        x,y=probplot(list_data[i])[0]
        plt.plot(x,y, 'o'+ colors[i], label=labels[i] )
        plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), colors[i])
        plt.xlabel('Theoretical quantiles', fontsize=18)
        plt.ylabel('Ordered Values', fontsize=16)
        plt.legend()
    return()
# %%


labels=['one', 'two', 'three'];
L=[np.diff(Z[604]), np.diff(X[604]), np.diff(X[350])];

plot_QQ(list_data=L, labels=labels)



L[0]


# V_centered=np.zeros_like(V)
# for i in range(0,V.shape[0]):
#     V_centered[i,:]=V[i,:]-m[i]
#
# m3=np.power(V_centered, 3).sum(axis=1)/V_centered.shape[1]
# m3

# E3=np.power(V, 3).sum(axis=1)/V.shape[1]
# E2=np.power(V, 2).sum(axis=1)/V.shape[1]
# E1=V.sum(axis=1)//V.shape[1]
# var=np.var(V,axis=1)
# s=E3-3*E1*var-E3**3
