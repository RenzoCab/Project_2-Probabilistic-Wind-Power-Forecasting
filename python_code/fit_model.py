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
#args.filename
config_file= open(args.filename)  #config/beta_config.JSON

setup=config.loadconfig.Test(config_file)

# os.chdir(setup.dir_path)

# orig_stdout = sys.stdout
# f = open( setup.logs_file_name, 'w')
# sys.stdout = f

def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))

warnings.showwarning = customwarn

warnings.simplefilter("once")

from Base_plus import *
forecast_data_in=np.load(setup.data_path)
#forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='fit_SDE'
current_data_dir= script_name +current_time
os.mkdir(current_data_dir)

forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)
print('Data output will be save in ',current_data_dir)


# #check no nans
# np.sum( np.isnan(forecast_data_inter.flatten()))
# np.max(forecast_data_inter[1:,:,:])<=1
# np.min(forecast_data_inter[1:,:,:])>=0

N=forecast_data_inter.shape[2]
M= setup.num_paths  #forecast_data_inter.shape[1]-973 #to be changed in generalization
dt=1

disct_temp = disct(N,dt,M)
p=forecast_data_inter[2,:-973,:]
V= forecast_data_inter[2,:-973,:]-forecast_data_inter[1,:-973,:]
X=forecast_data_inter[1,:-973,:]


if setup.likelihood=='lamperti_likelihood_SDE_approx':
    V = np.arcsin(2*X - 1)
    print(' Data has been Lamperti transformed ')

if setup.likelihood=='lamperti_likelihood_linearized':
    V = np.arcsin(2*X - 1)
    print(' Data has been Lamperti transformed ')


this_model=model_modified_drift(disct_temp,V, forecast= p)
intial_point= np.array(( setup.optimzation_initial_point['theta_init'], setup.optimzation_initial_point['alpha_init']))

print( ' optimizing ' + setup.likelihood + ' using ' + setup.optimizer )

optim=this_model.optimize(inference= setup.likelihood , method=setup.optimizer ,param_initial=intial_point)


print(optim)


file_object  = open(current_data_dir+'/results.out', 'w')
note='Fit for ' + str(M) +' paths' + ' of ' + str(N) +' hours. Initialized at '\
    + str(intial_point)
file_object.write( note +'\n'+str(optim.x) +'\n' +optim.message[0])
file_object.close()

print('results have been saved in ' + current_data_dir + '/results.out' )
