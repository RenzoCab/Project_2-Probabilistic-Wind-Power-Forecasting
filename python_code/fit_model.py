import os
import sys
os.chdir(sys.path[0]) #+'/python_code'
from Base_plus import *
import config.loadconfig
import argparse
#import numdifftools as nd
from matplotlib import animation
from IPython.display import HTML
from matplotlib.colors import LogNorm
from matplotlib import ticker
import time
from shutil import copyfile

########################
#Warning control
warnings.filterwarnings('error', '.*invalid value encountered.*',)

parser = argparse.ArgumentParser(description='Likelihood Evaluator v1.0')
parser.add_argument('-filename', help=' Config file name or path')
parser.add_argument('--version', action='version', version='Likelihood Evaluator v1.0')
args = parser.parse_args()

print( ' loading configuration ')
#args.filename
config_file= open(args.filename)  #config/beta_config.JSON


setup=config.loadconfig.Test(config_file)


######################### Likelihood background

#model of Beta
base_path=current_plotting_dir='data/'
chosen_folder='likelihood_explorer_'+ '19-09-11-21-32-05'

data_path=base_path+chosen_folder

#plotting chosen_folder
Z_grid=np.load(data_path+'/value.npy')
X_grid=np.load(data_path+'/theta.npy')
Y_grid=np.load(data_path+'/alpha.npy')


M=np.load(data_path+'/num_paths.npy')
N=np.load(data_path+'/interpolation_points.npy')


# %matplotlib qt


Z_grid=Z_grid-1*np.min(Z_grid) + 1e-16 #to make Z positive, we use an offset for plotting

fig = plt.figure(figsize=(8, 5))
ax = plt.axes(projection='3d', elev=50, azim=-50)
ax.plot_surface(X_grid, Y_grid, Z_grid, norm=LogNorm(), rstride=1, cstride=1, edgecolor='none', alpha=.8, cmap=plt.cm.jet)
# ax.plot(10,0.1, 10, 'r*', markersize=10)
line, = ax.plot([], [], [], 'b', label='Newton-CG', lw=2)
point, = ax.plot([], [], [], 'bo')
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
ax.set_zlabel('Likelihood')

##########################################

# os.chdir(setup.dir_path)

# orig_stdout = sys.stdout
# f = open( setup.logs_file_name, 'w')
# sys.stdout = f

# def customwarn(message, category, filename, lineno, file=None, line=None):
#     sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))
#
# warnings.showwarning = customwarn

warnings.simplefilter("once")

forecast_data_in=np.load(setup.data_path)
#forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='fit_SDE'
current_data_dir= script_name +current_time
os.makedirs(current_data_dir,exist_ok=True)

forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)
print('Data output will be save in ',current_data_dir)

# WARNING: Notice that we redefine the length in the transitions-direction for each vector.
# It was 72 hrs \approx 427 transitions, now (22/11/2019) it is 6 hrs \approx 59 transitions.
N=forecast_data_inter.shape[2] - 368 # 427 - 368 = 59, we do this_model because
# we only want the first 6 hours for each path to avoid overlaping.
# 427 # If we are considering 3 24-hours paths, notice that N < 72*6 because
# of the repetition at the bopundries of each day.
M= setup.num_paths  #forecast_data_inter.shape[1]-973 #to be changed in generalization
dt=1
p=forecast_data_inter[2,:setup.num_paths,:N]
p = p[::2,:]
# V= forecast_data_inter[2,:setup.num_paths,:]-forecast_data_inter[1,:setup.num_paths,:]
V= forecast_data_inter[1,:setup.num_paths,:N] - forecast_data_inter[2,:setup.num_paths,:N]
V = V[::2,:]
X=forecast_data_inter[1,:setup.num_paths,:N]
X = X[::2,:]
M = V.shape[0]

disct_temp = disct(N,dt,M) # We want to follow V = X - P \implies X = V + P


######## ToComplete:
# def get_datetime64(time_stamp):
#     return np.datetime64(int(time_stamp),'s')
# get_datetime64=np.vectorize(get_datetime64)
# print(forecast_data_inter[0,:setup.num_paths,0])
# print(get_datetime64(forecast_data_inter[0,:setup.num_paths,0]))
######## ToComplete:
# D=V
# O=np.sum(np.diff(D, axis= 1 )**2, axis=1)
# K= np.sum(2*(D+p)*(1-D-p) , axis=1)
# R= np.mean(O/K )
# # print( 'HERE ', R, K, O, O/K, D, p)
# print( 'HERE ', R)
# aux1 = 0; aux2 = 0; aux3 = 0
# for j in range(0,M):
#     for i in range(0,N-1):
#
#         aux1 += (D[j,i+1]-D[j,i])**2
#         aux2 += 2*(D[j,i]+p[j,i])*(1-D[j,i]-p[j,i])
#
#     aux3 += aux1 / aux2
# result = aux3 / M
# print('HERE ',result)



if setup.likelihood=='lamperti_likelihood_SDE_approx':
    V = np.arcsin(2*X - 1)
    print(' Data has been Lamperti transformed ')

if setup.likelihood=='lamperti_likelihood_linearized':
    V = np.arcsin(2*X - 1)
    print(' Data has been Lamperti transformed ')

this_model=model_modified_drift(disct_temp,V, forecast= p)

# print(this_model.expected_cross_time()) # We print this to test the function expected_cross_time.


file_object  = open(current_data_dir+'/results.out', 'w')

copyfile('config/beta_config_test.JSON', current_data_dir + '/beta_config_test.JSON' ) # Needs to be changed.

# width,height,angle, eig_vect_opt,Hess, FAIL = this_model.compute_ellipse(inference=setup.likelihood, param=(12, 0.4 ), batch_size=10 , plot=False);


#initialize parameters
current_batch_size=setup.optimization['initial_batch_size'];
intial_point=np.array(( setup.optimization['theta_init'] , setup.optimization['alpha_init']  ))

print('Computing Hessian at intial point')
#compute hessian at initial point

# width,height,angle, eig_vect_opt,Hess, FAIL = this_model.compute_ellipse(inference=setup.likelihood, param=(setup.optimization['theta_init'], setup.optimization['alpha_init']), batch_size=current_batch_size , plot=False);

# THIS IS WRONG, WE NEED TO USE THE SAME BATCH AS IN THE OPTIMIZOR:

# gen_mini_batch(self, batch_size, send_end_to_evaluater=None)
batch = this_model.gen_mini_batch(setup.optimization['initial_batch_size'], send_end_to_evaluater=None)
# likelihood_evaluate(self, selector, param, batch ,batch_size=None)
likelihood_initial_value= this_model.likelihood_evaluate(selector = setup.likelihood, param = (setup.optimization['theta_init'], setup.optimization['alpha_init']), batch = batch, batch_size = setup.optimization['initial_batch_size'])

#creating initial array
# parmeter_convergence=np.array((setup.optimization['theta_init'],setup.optimization['alpha_init'], likelihood_initial_value , setup.optimization['initial_batch_size'], height, width,angle ))
# hessian_convergance=Hess
# eigen_vect_convergence=eig_vect_opt

parmeter_convergence_limited=np.array((setup.optimization['theta_init'],setup.optimization['alpha_init'], likelihood_initial_value , setup.optimization['initial_batch_size'] ))

#printing
print( ' optimizing ' + setup.likelihood + ' using ' + setup.optimizer + ' with ' + str(setup.optimization['initial_batch_size']) + ' samples in the batch.'  )
# print( 'Ellipse ' + ' H: ',height , 'W: ', width )


while current_batch_size <= setup.optimization['max_batch_size']:
    #optimize
    optim=this_model.optimize(inference=setup.likelihood, batch_size=current_batch_size, method=setup.optimizer, param_initial=intial_point, niter=1, temp=0);

    print(optim.message[0])
    print(optim.x)

    result_note='optimization result for a batch of size ' + str(current_batch_size) + ' initialized at ' + str(intial_point) + ' is ' + str(optim.x) + ' with functional value ' + str(optim.fun) + ' with message ' + optim.message[0] + ' results save in ' + current_data_dir;

    print(result_note)

    file_object.write(result_note)

    # print('Computing Hessian at ', intial_point)

    # width,height,angle, eig_vect_opt,Hess, FAIL = this_model.compute_ellipse(inference=setup.likelihood, param=optim.x, batch_size=current_batch_size, plot=False);
    #
    # if FAIL==False:
    #     print( 'Ellipse ' + ' H: ', height , 'W: ', width )
    #     #save parameters
    #     parmeter_convergence=np.vstack((parmeter_convergence,np.hstack((optim.x,optim.fun,current_batch_size, height, width, angle) )));
    #     hessian_convergance=np.vstack((hessian_convergance, Hess))
    #     eigen_vect_convergence=np.vstack((eigen_vect_convergence, eig_vect_opt))
    #
    #     np.save(current_data_dir+'/parmeter_convergence_' + str(current_batch_size), parmeter_convergence);
    #     np.save(current_data_dir+'/eig_vect_opt_' + str(current_batch_size), eigen_vect_convergence);
    #     np.save(current_data_dir+'/Hess_' + str(current_batch_size) , hessian_convergance);
    #     print(' results save in ' + current_data_dir)

    parmeter_convergence_limited=np.vstack((parmeter_convergence_limited,np.hstack((optim.x,optim.fun,current_batch_size))));
    np.save(current_data_dir+'/parmeter_convergence_limited_' + str(current_batch_size), parmeter_convergence_limited);
    print(' results save in ' + current_data_dir)

    # #plot likelihood background
    # fig = plt.figure(figsize=(8, 5))
    # ax = plt.axes(projection='3d', elev=50, azim=-50)
    # ax.plot_surface(X_grid, Y_grid, Z_grid, norm=LogNorm(), rstride=1, cstride=1, edgecolor='none', alpha=.8, cmap=plt.cm.jet)
    # # ax.plot(10,0.1, 10, 'r*', markersize=10)
    # line, = ax.plot([], [], [], 'b', label='different Nelder-Mead instances', lw=2)
    # point, = ax.plot([], [], [], 'bo')
    # plt.xlabel('$\\theta$');
    # plt.ylabel('$\\alpha$');
    # ax.set_zlabel('Likelihood')
    #animate
    # path=parmeter_convergence_limited
    # # path=np.swapaxes(path, 0, 1)
    # anim = animation.FuncAnimation(fig, lambda i:animate(line, point,path,i) , init_func=lambda: init(line, point), frames=path.shape[1], interval=1000, repeat_delay=5, blit=True)
    # # writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    # anim.save(current_data_dir+'/convergance_'+ str(current_batch_size)  +'.mp4')

    #update batch size and intial point
    current_batch_size= current_batch_size*setup.optimization['batch_multiplier']
    intial_point=optim.x
    start_note='startimg fit for a batch of size ' + str(current_batch_size) + ' initialized at ' + str(intial_point);
    file_object.write(start_note)
    print(start_note)

file_object.close()


print('results have been saved in ' + current_data_dir + '/results.out' )


######################################################################################################
