from multi_path_base import *
from matplotlib.patches import Ellipse
import matplotlib as mpl
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

N=forecast_data_inter.shape[2]

dt=1
dN=1/interpolation_points

eig_theta=[]
eig_alpha=[]
M= 100
h=1e-2
for M in [100,1000]:
    disct_temp = disct(N,dt,M)
    p=forecast_data_inter[0,:,:]
    V= forecast_data_inter[1,:,:] -forecast_data_inter[0,:,:]
    this_model=model_modified_drift(disct_temp,V, forecast= p)
    h=1e-3
    x=8
    y=1

    Z_xx=(this_model.likelihood(param=np.array((x-h,y)))- 2*this_model.likelihood(param=np.array((x,y))) + this_model.likelihood(param=np.array((x+h,y))))/h**2
    Z_yy=(this_model.likelihood(param=np.array((x,y-h)))- 2*this_model.likelihood(param=np.array((x,y))) + this_model.likelihood(param=np.array((x,y+h))))/h**2
    Z_xy=(this_model.likelihood(param=np.array((x+h,y+h)))- this_model.likelihood(param=np.array((x+h,y)))- this_model.likelihood(param=np.array((x,y+h))) + this_model.likelihood(param=np.array((x,y))))/h**2

    H_opt = np.zeros((2,2))
    H_opt[0,0]= Z_xx
    H_opt[1,1]= Z_yy
    H_opt[1,0]= Z_xy
    H_opt[0,1]= Z_xy
    eig_values_opt, eig_vect_opt=LA.eig(H_opt)
    eig_theta.append(eig_values_opt[0])
    eig_alpha.append(eig_values_opt[1])

eig_theta

eig_alpha

plt.plot(np.log10([100,1000]), np.log10(eig_theta))
x=np.arange(2,4)
y=np.arange(2,4)
plt.plot(x,y,'--', label='slope 1')
plt.plot(x,2*y,'--', label='slope 2')
plt.plot(np.log10([100,1000]), np.log10(eig_alpha))



def Ellipse (H_opt, mean=[ 8 ,  1]):
    eig_values_opt, eig_vect_opt=LA.eig(H_opt)
    a_opt=1/np.sqrt(np.max(np.abs(eig_values_opt)))  #radius on the x-axis
    b_opt=1/np.sqrt( np.min(np.abs(eig_values_opt))) #radius on the y-axis
    width = a_opt**2
    height = b_opt**2
    angle =math.acos(sum(eig_vect_opt[1]*np.array((0,1))));angle
    return(width,height,angle, eig_vect_opt )

mean=[ x ,  y]

fig, ax = plt.subplots()

width,height,angle, eig_vect_opt=Ellipse(H_opt)
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.1, color='b')
V = np.array([-1*eig_vect_opt[0],-1*eig_vect_opt[1]])
origin = [x], [y] # origin point
plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
ax.add_patch(ell)
ax.autoscale()


###### Convergence testing of finite differences
h=[0.5, 0.3,1e-1,1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7,1e-8, 1e-9]
Z_xx_list=[]
Z_yy_list=[]
Z_xy_list=[]
for i in range(0,11):
    Z_yy=0
    Z_xx=0
    Z_xx=(this_model.likelihood(param=np.array((x-h[i],y)))- 2*this_model.likelihood(param=np.array((x,y))) + this_model.likelihood(param=np.array((x+h[i],y))) )/h[i]**2
    Z_xx_list.append(Z_xx)
    Z_yy=(this_model.likelihood(param=np.array((x,y-h[i])))- 2*this_model.likelihood(param=np.array((x,y))) + this_model.likelihood(param=np.array((x,y+h[i]))))/h[i]**2
    Z_yy_list.append(Z_yy)
    Z_xy=(this_model.likelihood(param=np.array((x+h[i],y+h[i])))- this_model.likelihood(param=np.array((x+h[i],y)))- this_model.likelihood(param=np.array((x,y+h[i]))) + this_model.likelihood(param=np.array((x,y))))/h[i]**2
    Z_xy_list.append(Z_xy)

h_list=[0.5, 0.3,1e-1,1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7,1e-8]

plt.plot(np.log10(h_list),np.log10(abs(np.diff(Z_xx_list))), '-o' )
Z_xx_list
plt.plot(np.log10(h_list),np.log10(abs(np.diff(Z_yy_list))), '-o' )
Z_yy_list
plt.plot(np.log10(h_list),np.log10(abs(np.diff(Z_xy_list))), '-o' )
Z_xy_list


H_opt = np.zeros((2,2))
H_opt[0,0]= Z_xx
H_opt[1,1]= Z_yy
H_opt[1,0]= Z_xy
H_opt[0,1]= Z_xy

eig_values_opt, eig_vect_opt=LA.eig(H_opt)

H2=np.transpose(H_opt)@H_opt
eig_H2,eignvect_H2=LA.eig(H2)

eig_H2

def Ellipse (H_opt, mean=[ 8 ,  1]):
    eig_values_opt, eig_vect_opt=LA.eig(H_opt)
    a_opt=1/np.sqrt(np.max(np.abs(eig_values_opt)))  #radius on the x-axis
    b_opt=1/np.sqrt( np.min(np.abs(eig_values_opt))) #radius on the y-axis
    width = a_opt**2
    height = b_opt**2
    angle =math.acos(sum(eig_vect_opt[1]*np.array((0,1))));angle
    return(width,height,angle, eig_vect_opt )

mean=[ x ,  y]
alpha=[0.2,0.4, 0.6, 0.8]

fig, ax = plt.subplots()
print(width)
print(height)
width,height,angle, eig_vect_opt=Ellipse(H2)
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.1, color='b')
V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
origin = [x], [y] # origin point
plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
ax.add_patch(ell)
ax.autoscale()
