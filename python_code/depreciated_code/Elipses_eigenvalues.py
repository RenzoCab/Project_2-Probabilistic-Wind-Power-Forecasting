from multi_path_base import *
from matplotlib.patches import Ellipse
import matplotlib as mpl

H_opt_10=np.load('H_opt_10.npy')
H_opt_100=np.load('H_opt_100.npy')
H_opt_500=np.load('H_opt_500.npy')
H_opt_1193=np.load('H_opt_1193.npy')

def Ellipse (H_opt, mean=[ 8 ,  1]):
    eig_values_opt, eig_vect_opt=LA.eig(H_opt)
    a_opt=1/np.sqrt(np.max(np.abs(eig_values_opt)))  #radius on the x-axis
    b_opt=1/np.sqrt( np.min(np.abs(eig_values_opt))) #radius on the y-axis
    width = a_opt**2
    height = b_opt**2
    angle =math.acos(sum(eig_vect_opt[1]*np.array((0,1))));angle
    return(width,height,angle, eig_vect_opt )

interpolation_points=73
dN= 1/interpolation_points
mean=[ 8 ,  1]
# alpha=[0.2,0.4, 0.6, 0.8]
# color=['b','r', 'k', 'm']
fig, ax = plt.subplots()
M=[100,500,1193]
j=0
for i in [H_opt_100,H_opt_500,H_opt_1193 ]:
    width,height,angle, eig_vect_opt=Ellipse(i)
    plt.clf()
    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.2, color='k')
    V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
    origin = [8], [1] # origin point
    plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
    ax.add_patch(ell)
    ax.autoscale()
    file_name='ellipse'+str(M[j]) +'_samples_dN='+'{:.1e}'.format(dN)
    plt.savefig( file_name + '.pdf')
    j+=1



mean=[ 8 ,  1]
M=1193
width,height,angle, eig_vect_opt=Ellipse(H_opt_1193)
plt.clf()
fig, ax = plt.subplots()
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.2, color='g')
V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
origin = [8], [1] # origin point
plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=10)
ax.add_patch(ell)
ax.autoscale()
plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
file_name='ellipse'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.title(' (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')',fontsize = 15); #Ellipse of the Hessian of Log-likelihood \n
plt.xlabel('$\\theta$',fontsize = 24);
plt.ylabel('$\\alpha$',fontsize = 24);
plt.xticks( fontsize = 20);
plt.yticks( fontsize = 20);
plt.savefig( file_name + '.pdf', bbox_inches="tight")


mean=[ 8 ,  1]
M=500
width,height,angle, eig_vect_opt=Ellipse(H_opt_500)
plt.clf()
fig, ax = plt.subplots()
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.2, color='k')
V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
origin = [8], [1] # origin point
plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
ax.add_patch(ell)
ax.autoscale()
plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
file_name='ellipse'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.title('Ellipse of the Hessian of Log-likelihood \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.savefig( file_name + '.pdf')


mean=[ 8 ,  1]
M=100
width,height,angle, eig_vect_opt=Ellipse(H_opt_100)
plt.clf()
fig, ax = plt.subplots()
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.2, color='k')
V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
origin = [8], [1] # origin point
plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
ax.add_patch(ell)
ax.autoscale()
plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
file_name='ellipse'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.title('Ellipse of the Hessian of Log-likelihood \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
plt.savefig( file_name + '.pdf')








j=0
Eig_theta=[]
Eig_alpha=[]
a_list=[]
b_list=[]
for i in [H_opt_100,H_opt_500,H_opt_1193 ]:
    eig_values_opt,_=LA.eig(i)
    Eig_theta.append(eig_values_opt[0])
    Eig_alpha.append(eig_values_opt[1])
    a_list.append(1/np.sqrt(np.max(np.abs(eig_values_opt))) ) #radius on the x-axis
    b_list.append(1/np.sqrt( np.min(np.abs(eig_values_opt)))) #radius on the y-axis
Eig_theta
Eig_alpha
M=[100,500,1193]
x=np.arange(2,4)
y=np.arange(2,4)
plt.plot(x,y,'--', label='slope 1')
plt.plot(x,y**1.5,'--', label='slope 1.5')
plt.plot(np.log10(M), np.log10(Eig_theta), '-ro', label='$\\lambda_\\theta$')
plt.plot(np.log10(M), np.log10(Eig_alpha), '-bo', label='$\\lambda_\\alpha$')
plt.legend()
plt.xlabel('$\\log(M)$');
plt.ylabel('$\\log(\\lambda)$');
np.polyfit(np.log10(M), np.log10(Eig_theta), 1)
np.polyfit(np.log10(M), np.log10(Eig_alpha), 1)
file_name='Eigen_conv'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf')



M=[100,500,1193]
x=np.arange(2,4)
y=np.arange(2,4)
plt.plot(x,-0.5*y,'--', label='slope -1/2')
plt.plot(x,-1*y,'--', label='slope -1')
plt.plot(np.log10(M), np.log10(a_list), '-ro', label='$r_\\theta$')
plt.plot(np.log10(M), np.log10(b_list), '-bo', label='$r_\\alpha$')
plt.xticks( fontsize = 20)
plt.yticks( fontsize = 20)
plt.legend()
plt.xlabel('$\\log(M)$',fontsize=24);
plt.ylabel('$\\log(\\lambda)$',fontsize=24);

file_name='ellipse_conv'+str(M) +'_samples_dN='+'{:.1e}'.format(dN)
plt.savefig( file_name + '.pdf', bbox_inches="tight")
