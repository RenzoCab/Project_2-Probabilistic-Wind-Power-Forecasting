# %%
from multi_path_base import *
from scipy import interpolate
forecast_with_data=np.load('forecast_with_data.npy')
# %%
#cleaning data
interpolation_points=73
forecast_data_inter= np.zeros((3,1193, interpolation_points ))
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
N=forecast_data_inter.shape[2]
M= 100    #forecast_data_inter.shape[1]
dt=1
dN=1/interpolation_points

#forecast_data_inter[0,:,:] #forecast

disct_temp = disct(N,dt,M)
p=forecast_data_inter[0,:,:]
V= forecast_data_inter[1,:,:] -forecast_data_inter[0,:,:]

this_model=model_modified_drift(disct_temp,V, forecast= p)
# %%

x = np.linspace(1, 30, 30)
y = np.linspace(0.1, 3, 30)

X, Y = np.meshgrid(x, y)

Z=np.zeros((len(y), len(x)))
for i in range(0,len(x)):
    for j in range(0,len(y)):
        Z[j,i] = this_model.likelihood(param=np.array((X[j,i],Y[j,i])))

max_Z = np.max(np.abs(Z));
max_Z
norm = cm.colors.Normalize(vmax= -1*abs(Z/max_Z).min(), vmin=-1*abs(Z/max_Z ).max())
cmap = cm.Spectral #RdYlGn
contours=plt.contourf(X, Y, Z/max_Z,50, cmap=cmap)
plt.colorbar()
plt.title(' (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')',fontsize=15); #,fontsize=24
plt.xlabel('$\\theta$',fontsize = 24);
plt.ylabel('$\\alpha$',fontsize = 24);
plt.xticks( fontsize = 20)
plt.yticks( fontsize = 20)
plt.scatter(8, 1, c='k', alpha=0.5)
file_name='ISO_'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf', bbox_inches="tight")



Zx, Zy=np.gradient(Z)

#Zx[0][0]
#Zy[0][0]

Z_grad_x=np.zeros((len(y), len(x)))
Z_grad_y=np.zeros((len(y), len(x)))
for i in range(0,len(x)):
    for j in range(0,len(y)):
        Z_grad_x[i,j] = Zx[i][j]
        Z_grad_y[i,j] = Zy[i][j]


cmap = cm.RdGy
max_Z_grad_x = np.max(np.abs(Z_grad_x));
norm = cm.colors.Normalize(vmax= abs(Z_grad_x/max_Z_grad_x).max(), vmin=abs(Z_grad_x/max_Z_grad_x ).min())
plt.clf()
contours=plt.contourf(X, Y, Z_grad_x/max_Z_grad_x, 50,norm=norm, cmap=cmap)
plt.title('Gradient of Beta Log-likelihood wrt $\\theta $ \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.colorbar()
file_name='ISO_Grad_x'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf')


max_Z_grad_y = np.max(np.abs(Z_grad_y));
#norm = cm.colors.Normalize(vmax= abs(Z_grad_y/max_Z_grad_y).max(), vmin=abs(Z_grad_y/max_Z_grad_y ).min())
plt.clf()
contours=plt.contourf(X, Y, Z_grad_y/max_Z_grad_y, 50,norm=norm, cmap=cmap)
plt.title('Gradient of Beta Log-likelihood wrt $\\alpha $ \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.colorbar()
file_name='ISO_Grad_y'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf')


Z_grad_norm=np.sqrt(Z_grad_x**2 + Z_grad_y**2)
cmap = cm.RdGy
max_Z_grad_norm = np.max(np.abs(Z_grad_norm));
norm = cm.colors.Normalize(vmax= abs(Z_grad_norm/max_Z_grad_norm).max(), vmin=abs(Z_grad_norm/max_Z_grad_norm ).min())
plt.clf()
contours=plt.contourf(X, Y, Z_grad_norm/max_Z_grad_norm, 50,norm=norm, cmap=cmap)
plt.title('$L_2$ Norm of the Gradient of Log-likelihood \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.colorbar()
file_name='ISO_Grad_norm'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf')



Zxx,Zxy =np.gradient(Z_grad_x)
Zyx,Zyy =np.gradient(Z_grad_y)
#Zxx[0][0]
#Zyy[0][0]

Z_grad_xx=np.zeros((len(y), len(x)))
Z_grad_yy=np.zeros((len(y), len(x)))
Z_grad_xy=np.zeros((len(y), len(x)))
Z_grad_yx=np.zeros((len(y), len(x)))
for i in range(0,len(x)):
    for j in range(0,len(y)):
        Z_grad_xx[i,j] = Zxx[i][j]
        Z_grad_yy[i,j] = Zyy[i][j]
        Z_grad_xy[i,j] = Zxy[i][j]
        Z_grad_yx[i,j] = Zyx[i][j]



cmap = cm.RdGy
max_Z_grad_xx = np.max(np.abs(Z_grad_xx));
norm = cm.colors.Normalize(vmax= abs(Z_grad_xx/max_Z_grad_xx).max(), vmin=abs(Z_grad_xx/max_Z_grad_xx ).min())
plt.clf()
contours=plt.contourf(X, Y, Z_grad_xx/max_Z_grad_xx, 50,norm=norm, cmap=cmap)
plt.title('Second derivative of Beta Log-likelihood wrt $\\theta $ \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.colorbar()
file_name='ISO_Grad_xx'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf')


max_Z_grad_yy = np.max(np.abs(Z_grad_yy));
#norm = cm.colors.Normalize(vmax= abs(Z_grad_y/max_Z_grad_y).max(), vmin=abs(Z_grad_y/max_Z_grad_y ).min())
plt.clf()
contours=plt.contourf(X, Y, Z_grad_yy/max_Z_grad_yy, 50,norm=norm, cmap=cmap)
plt.title('Second derivative of Beta Log-likelihood wrt $\\alpha $ \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.colorbar()
file_name='ISO_Grad_yy'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf')




Z_hess_norm= Z_grad_xx**2 + Z_grad_yy**2 + Z_grad_xy**2 + Z_grad_yx**2
cmap = cm.RdGy
max_Z_hess_norm = np.max(np.abs(Z_hess_norm));
norm = cm.colors.Normalize(vmax= abs(Z_hess_norm/max_Z_hess_norm).max(), vmin=abs(Z_hess_norm/max_Z_hess_norm ).min())
plt.clf()
contours=plt.contourf(X, Y, Z_hess_norm/max_Z_hess_norm, 50,norm=norm, cmap=cmap)
plt.title('$L_2$ Norm of the Hessian of Log-likelihood \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.colorbar()
file_name='ISO_hess_norm'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf')





from matplotlib.patches import Ellipse
import matplotlib as mpl

H=np.zeros((60,60))
H[:30,:30]= Z_grad_xx
H[30:,30:]= Z_grad_yy
H[:30,30:]= Z_grad_yx
H[30:,:30]= Z_grad_xy

eig_values, eig_vect=LA.eig(H)

H_opt = np.zeros((2,2))
H_opt[0,0]= Z_grad_xx[7,9]
H_opt[1,1]= Z_grad_yy[7,9]
H_opt[1,0]= Z_grad_yx[7,9]
H_opt[0,1]= Z_grad_xy[7,9]

np.save('H_opt_100', H_opt)

eig_values_opt, eig_vect_opt=LA.eig(H_opt)

a_opt=1/np.sqrt(np.max(np.abs(eig_values_opt)))  #radius on the x-axis
b_opt=1/np.sqrt( np.min(np.abs(eig_values_opt))) #radius on the y-axis


mean = [ 8 ,  1]
width = a_opt**2
height = b_opt**2
b_opt
angle =math.acos(sum(eig_vect_opt[1]*np.array((0,1))));angle
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.1)
fig, ax = plt.subplots()

V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
origin = [8], [1] # origin point

plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
ax.add_patch(ell)
ax.autoscale()
plt.show()
