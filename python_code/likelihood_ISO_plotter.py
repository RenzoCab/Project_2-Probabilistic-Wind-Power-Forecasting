import os
os.chdir('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Ellipse
import matplotlib as mpl
from numpy import linalg as LA
import math

def Ellipse (H_opt, mean=[ 8 ,  1]):
    eig_values_opt, eig_vect_opt=LA.eig(H_opt)
    a_opt=1/np.sqrt(np.max(np.abs(eig_values_opt)))  #radius on the x-axis
    b_opt=1/np.sqrt( np.min(np.abs(eig_values_opt))) #radius on the y-axis
    width = a_opt**2
    height = b_opt**2
    angle =math.acos(sum(eig_vect_opt[1]*np.array((0,1))));angle
    return(width,height,angle, eig_vect_opt )

#19-07-26-15-30-16 Beta
#model of lamperti linearizing fit_SDE 19-09-04-15-54-23
#model using approx lamperti 19-09-04-16-11-26
#model of Beta
base_path=current_plotting_dir='data/'
chosen_folder='likelihood_explorer_'+ '19-09-11-21-32-05'
#beta 19-09-08-13-58-19

cd ps
data_path=base_path+chosen_folder

#plotting chosen_folder

os.mkdir('plots/'+chosen_folder)

Z=np.load(data_path+'/value.npy')
X=np.load(data_path+'/theta.npy')
Y=np.load(data_path+'/alpha.npy')
Z_grad_xx=np.load(data_path+'/grad_xx.npy')
Z_grad_yy=np.load(data_path+'/grad_yy.npy')
Z_grad_xy=np.load(data_path+'/grad_xy.npy')
Z_grad_yx=np.load(data_path+'/grad_yx.npy')

M=np.load(data_path+'/num_paths.npy')
N=np.load(data_path+'/interpolation_points.npy')

dN=1/N



#max_Z = np.max(np.abs(Z));
#max_Z
#norm = cm.colors.Normalize(vmax= -1*abs(Z/max_Z).min(), vmin=-1*abs(Z/max_Z ).max())
cmap = cm.RdGy
#contours=plt.contourf(X[5:20,5:20], Y[5:20, 5:20], Z[5:20, 5:20],200, cmap=cmap)
contours=plt.contourf(X, Y, Z,100, cmap=cmap)
plt.colorbar()
plt.scatter(22.32967705,  0.04929445)
plt.title('Beta Log-likelihood \n (' + str(M) +' samples, $\\Delta N=1/$'+str(N)+')'); #,fontsize=24
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
file_name='ISO_'+str(M) +'_inter='+str(N);
plt.savefig( 'plots/'+chosen_folder+'/' +file_name+'.pdf', bbox_inches="tight")
print('plots/'+chosen_folder+'/' +file_name+'.pdf')

# H=np.zeros((60,60))
# H[:30,:30]= Z_grad_xx;
# H[30:,30:]= Z_grad_yy;
# H[:30,30:]= Z_grad_yx;
# H[30:,:30]= Z_grad_xy;
#
# eig_values, eig_vect=LA.eig(H)
#
# H_opt = np.zeros((2,2))
# H_opt[0,0]= Z_grad_xx[7,9]
# H_opt[1,1]= Z_grad_yy[7,9]
# H_opt[1,0]= Z_grad_yx[7,9]
# H_opt[0,1]= Z_grad_xy[7,9]
#
#
# np.save('H_opt_10', H_opt)
#
# eig_values_opt, eig_vect_opt=LA.eig(H_opt)
#
# a_opt=1/np.sqrt(np.max(np.abs(eig_values_opt)))  #radius on the x-axis
# b_opt=1/np.sqrt( np.min(np.abs(eig_values_opt))) #radius on the y-axis
#
#
# mean = [ 8 ,  1]
# width = a_opt**2
# height = b_opt**2
# b_opt
# angle =math.acos(sum(eig_vect_opt[1]*np.array((0,1))));angle
# ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.1)
# fig, ax = plt.subplots()
#
# V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
# origin = [8], [1] # origin point
#
# plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
# ax.add_patch(ell)
# ax.autoscale()
# plt.savefig( 'plots/'+chosen_folder+'/' +file_name+'.pdf', bbox_inches="tight")
# #plt.show()
#
# mean=[ 8 ,  1]
# width,height,angle, eig_vect_opt=Ellipse(H_opt)
# plt.clf()
# fig, ax = plt.subplots()
# ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=0.2, color='k')
# V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
# origin = [8], [1] # origin point
# plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=4)
# ax.add_patch(ell)
# ax.autoscale()
# plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
# plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
# file_name='ellipse'+str(M) +'_interp='+ str(N)
# plt.title('Ellipse of the Hessian of Log-likelihood \n (' + str(M) +' samples, $\\Delta N= $'+'{:.1e}'.format(dN)+')');
# plt.xlabel('$\\theta$');
# plt.ylabel('$\\alpha$');
# plt.savefig( 'plots/'+chosen_folder+'/'+file_name + '.pdf')
