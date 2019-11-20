import os
import sys
os.chdir(sys.path[0])
from Base_plus import *
import config.loadconfig
import argparse
#import numdifftools as nd

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

# def customwarn(message, category, filename, lineno, file=None, line=None):
#     sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))
#
# warnings.showwarning = customwarn

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

file_object  = open(current_data_dir+'/results.out', 'w')

#initialize parameters
current_batch_size=setup.optimization['initial_batch_size'];
intial_point=np.array(( setup.optimization['theta_init'] , setup.optimization['alpha_init']  ))

print('Computing Hessian at intial point')
#compute hessian at initial point
width,height,angle, eig_vect_opt,Hess = this_model.compute_ellipse(inference=setup.likelihood, param=(setup.optimization['theta_init'], setup.optimization['alpha_init']), batch_size=current_batch_size , plot=False);

#creating initial array
parmeter_convergence=np.array((setup.optimization['theta_init'],setup.optimization['alpha_init'],setup.optimization['initial_batch_size'], height, width,angle ))
hessian_convergance=Hess
eigen_vect_convergence=eig_vect_opt


#printing
print( ' optimizing ' + setup.likelihood + ' using ' + setup.optimizer + ' with ' + str(setup.optimization['initial_batch_size']) + ' batches.'  )
print( 'Ellipse ' + ' H: ',height , 'W: ', width )


while current_batch_size <= setup.optimization['max_batch_size']:
    #optimize
    optim=this_model.optimize(inference=setup.likelihood, batch_size=current_batch_size, method=setup.optimizer, param_initial=intial_point);

    print(optim.message[0])
    print(optim.x)

    result_note='optimization result for a batch of size ' + str(current_batch_size) + ' initialized at ' + str(intial_point) + ' is ' + str(optim.x) + ' with message ' + optim.message[0] ;
    file_object.write(result_note)

    print('Computing Hessian at ' + intial_point)

    width,height,angle, eig_vect_opt,Hess = this_model.compute_ellipse(inference=setup.likelihood, param=optim.x, batch_size=current_batch_size, plot=False);

    print( 'Ellipse ' + ' H: ',height , 'W: ', width )

    #save parameters
    parmeter_convergence=np.vstack((parmeter_convergence,np.hstack((optim.x,current_batch_size, height, width, angle) )));
    hessian_convergance=np.vstack((hessian_convergance, Hess))
    eigen_vect_convergence=np.vstack((eigen_vect_convergence, eig_vect_opt))

    np.save('parmeter_convergence_' + str(current_batch_size), parmeter_convergence);
    np.save('eig_vect_opt_' + str(current_batch_size), eigen_vect_convergence);
    np.save('Hess_' + str(current_batch_size) , hessian_convergance);

    #update batch size and intial point
    current_batch_size= current_batch_size*setup.optimization['batch_multiplier']
    intial_point=optim.x
    start_note='startimg fit for a batch of size ' + str(current_batch_size) + ' initialized at ' + str(intial_point);
    file_object.write(start_note)
    print(start_note)

file_object.close()

print('results have been saved in ' + current_data_dir + '/results.out' )
