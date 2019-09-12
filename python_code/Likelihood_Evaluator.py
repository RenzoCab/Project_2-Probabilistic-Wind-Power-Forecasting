import os
import sys
os.chdir(sys.path[0])
from Base_plus import *
import config.loadconfig
import argparse
parser = argparse.ArgumentParser(description='Likelihood Evaluator v1.0')
parser.add_argument('-filename', help=' Config file name or path')
parser.add_argument('--version', action='version', version='Likelihood Evaluator v1.0')
args = parser.parse_args()

print( ' loading configuration ')

config_file= open(args.filename)  #config/beta_config.JSON

setup=config.loadconfig.Test(config_file)

#os.chdir(setup.dir_path)

print( ' starting ' + setup.likelihood + ' evaluator ')
orig_stdout = sys.stdout
f = open('logs/' + setup.logs_file_name, 'w')
sys.stdout = f

def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))

warnings.showwarning = customwarn

warnings.simplefilter("once")


# sys.stdin.close()
# sys.stdin = open(os.devnull)


#from multi_path_base import model_modified_drift
forecast_data_in=np.load(setup.data_dir)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='likelihood_explorer'
current_plotting_dir='data/'+ script_name+'_' +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)

N=forecast_data_inter.shape[2]
M=setup.num_paths
dt=1
dN=1/N
disct_temp = disct(N,dt,M)
p=forecast_data_inter[2,:-973,:]  #240
V=forecast_data_inter[2,:-973,:] -forecast_data_inter[1,:-973,:]

this_model=model_modified_drift(disct_temp,V, forecast= p)


theta_array=np.linspace(setup.theta['start'], setup.theta['end'], setup.theta['num_steps'])
alpha_array=np.linspace(setup.alpha['start'], setup.alpha['end'], setup.alpha['num_steps'])

np.save(current_plotting_dir+'/num_paths', M)
np.save(current_plotting_dir+'/interpolation_points', N)

if setup.num_cores=="all":
    likelihood_evaluater_pool(theta_array=theta_array,\
        alpha_array=alpha_array,\
        this_model=this_model,inference=setup.likelihood, path_dir=current_plotting_dir )
else:
    likelihood_evaluater_pool(theta_array=theta_array,\
        alpha_array=alpha_array,\
        this_model=this_model,inference=setup.likelihood, path_dir=current_plotting_dir,num_cores=setup.num_cores )

sys.stdout = orig_stdout
f.close()

print( 'output arrays are successfully saved in ', current_plotting_dir)
