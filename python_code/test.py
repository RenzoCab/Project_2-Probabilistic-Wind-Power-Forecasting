import matplotlib.pyplot as plt
# import autograd.numpy as np

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from matplotlib import animation
from IPython.display import HTML

# from autograd import elementwise_grad, value_and_grad
from scipy.optimize import minimize
from collections import defaultdict
from itertools import zip_longest
from functools import partial

import os
os.chdir('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Ellipse
import matplotlib as mpl
from numpy import linalg as LA
import math



# fig, ax = plt.subplots(figsize=(10, 6))
#
# ax.contour(x, y, z, levels=np.logspace(0, 5, 35), norm=LogNorm(), cmap=plt.cm.jet)
# ax.plot(*minima_, 'r*', markersize=18)
#
# line, = ax.plot([], [], 'b', label='Newton-CG', lw=2)
# point, = ax.plot([], [], 'bo')
#
# ax.set_xlabel('$x$')
# ax.set_ylabel('$y$')
#
# ax.set_xlim((xmin, xmax))
# ax.set_ylim((ymin, ymax))
#
# ax.legend(loc='upper left')



#model of Beta
base_path=current_plotting_dir='data/'
chosen_folder='likelihood_explorer_'+ '19-09-11-21-32-05'

data_path=base_path+chosen_folder

#plotting chosen_folder
Z=np.load(data_path+'/value.npy')
X=np.load(data_path+'/theta.npy')
Y=np.load(data_path+'/alpha.npy')


M=np.load(data_path+'/num_paths.npy')
N=np.load(data_path+'/interpolation_points.npy')


from matplotlib import ticker, cm
%matplotlib qt
fig = plt.figure(1)

cmap = cm.coolwarm
contours=plt.contourf(X, Y, Z ,600, cmap=cmap)
plt.colorbar()
plt.plot(23.95875, 0.34)
plt.title('Beta Log-likelihood \n (' + str(M) +' samples, $\\Delta N=1/$'+str(N)+')');
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
file_name='ISO_'+str(M) +'_inter='+str(N);

plt.savefig( 'plots/'+chosen_folder+'/' +file_name+'.pdf', bbox_inches="tight")
print('plots/'+chosen_folder+'/' +file_name+'.pdf')






ax.contour(x, y, z, levels=np.logspace(0, 5, 35), norm=LogNorm(), cmap=plt.cm.jet)
