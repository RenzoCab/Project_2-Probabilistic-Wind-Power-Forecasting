import os
import sys
os.chdir(sys.path[0])
from Base_plus import *
import config.loadconfig
import itertools as iter
# import argparse
# parser = argparse.ArgumentParser(description='Likelihood Evaluator v1.0')
# parser.add_argument('-filename', help=' Config file name or path')
# parser.add_argument('--version', action='version', version='Likelihood Evaluator v1.0')
# args = parser.parse_args()

print( ' loading configuration ')

config_file= open('beta_config.JSON')  #config/beta_config.JSON

setup=config.loadconfig.Test(config_file)

#os.chdir(setup.dir_path)

print( ' starting ' + setup.likelihood + ' evaluator ')

# orig_stdout = sys.stdout
# f = open('logs/' + setup.logs_file_name, 'w')
# sys.stdout = f
#
# def customwarn(message, category, filename, lineno, file=None, line=None):
#     sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))
#
# warnings.showwarning = customwarn

warnings.simplefilter("once")


# sys.stdin.close()
# sys.stdin = open(os.devnull)


#from multi_path_base import model_modified_drift
forecast_data_in=np.load(setup.data_path)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='likelihood_explorer'
current_plotting_dir='data/'+ script_name+'_' +current_time
# os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)

N=forecast_data_inter.shape[2]
M=setup.num_paths
dt=1
dN=1/N
disct_temp = disct(N,dt,M)
p=forecast_data_inter[2,:-240,:]  #240
V=forecast_data_inter[2,:-240,:] -forecast_data_inter[1,:-240,:]



def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = iter.tee(iterable)
    next(b, None)
    return zip(a, b)

def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)

def count_iterable(i):
    return sum(1 for e in i)


# theta=10
# M=100
# # random.seed(1)
#
# batch_size=50;
# ##############
#
# list_pairs_V=[]
# list_pairs_p=[]
# list_pairs_theta=[]
#
# len(theta_adjust( theta , p[1,1:] )[0])
# len(V[1,1:-1])
#
# for i in range(M):
#     list_pairs_V.append(pairwise( V[i,1:-1]))
#     list_pairs_p.append(pairwise( p[i,1:-1]))
#     list_pairs_theta.append(pairwise( theta_adjust( theta , p[i,1:] )[0] ) )
#
#
# obs_pairs=iter.chain(*list_pairs_V);
# forecast_pairs=iter.chain(*list_pairs_p);
# theta_pairs=iter.chain(*list_pairs_theta);
#
#
# combined_iter=iter.zip_longest(obs_pairs,forecast_pairs,theta_pairs)
# batch=random_combination(combined_iter,10000)
# ##############
#
#
# isinstance(X, (type(None), bytes))
#
#
# any(list(sum(X, ())) ) ==None
#
# if np.any( sum(X, ()) is None) : print('ALERT ! None type ')
#
# if np.any(X is None): print('None')
#
#
# batch[0]

10000/425

this_model=model_modified_drift(disct_temp,V, forecast= p)


objective_arg=(10,0.1,10000)

optim=this_model.optimize(param_initial=objective_arg,inference="rand_beta_objective" ,method="Nelder-Mead", niter=1)



theta_array=np.linspace(setup.theta['start'], setup.theta['end'], setup.theta['num_steps'])
alpha_array=np.linspace(setup.alpha['start'], setup.alpha['end'], setup.alpha['num_steps'])

# np.save(current_plotting_dir+'/num_paths', M)
# np.save(current_plotting_dir+'/interpolation_points', N)

# if setup.num_cores=="all":
#     likelihood_evaluater_pool(theta_array=theta_array,\
#         alpha_array=alpha_array,\
#         this_model=this_model,inference=setup.likelihood, path_dir=current_plotting_dir )
# else:
#     likelihood_evaluater_pool(theta_array=theta_array,\
#         alpha_array=alpha_array,\
#         this_model=this_model,inference=setup.likelihood, path_dir=current_plotting_dir,num_cores=setup.num_cores )

# sys.stdout = orig_stdout
# f.close()

# print( 'output arrays are successfully saved in ', current_plotting_dir)
