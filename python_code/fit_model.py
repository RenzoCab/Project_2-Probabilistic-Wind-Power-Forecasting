import os
import sys
os.chdir(sys.path[0]) #+'/python_code'
from Base_plus import *
import config.loadconfig
import argparse # https://docs.python.org/3/library/argparse.html
#import numdifftools as nd
from matplotlib import animation
from IPython.display import HTML
from matplotlib.colors import LogNorm
from matplotlib import ticker
import time
from shutil import copyfile

######################### Warning Control and Parser #########################

warnings.filterwarnings('error', '.*invalid value encountered.*',)

if 1 == 1: # IGNORE THESE LINES IF WE RUN LINE BY LINE:
    ## WE CREATE THE PARSER:
    # For the next lines, we follow https://docs.python.org/3/library/argparse.html.
    # This line creates the parser, which contains description, data, and version:
    parser = argparse.ArgumentParser(description = 'Likelihood Evaluator v1.0')
    # This line is for loading the location of the configuration file (.JSON file with the configuration):
    parser.add_argument('-filename', help = ' Config file name or path')
    # The last configuration that worked is: python fit_model.py -f config/beta_config.JSON.
    # This line is for creating a version command, i.e., python fit_model.py --version:
    parser.add_argument('--version', action = 'version', version = 'Likelihood Evaluator v1.0')
    args = parser.parse_args()

    ## WE LOAD THE FILE .JSON WITH ALL THE SPECIFICATIONS:
    print('Loading Configuration...')
    # We load the filename from args into the variable config_file:
    config_file = open(args.filename)  # filename = config/beta_config.JSON.
    # Finally... We read the file called filename where the configuration is.
    setup = config.loadconfig.Test(config_file)
    # Now, the setup variable contains all the variables in beta_config.JSON.
else: # We run these lines if we are running line by line, and we can not insert filename.
    config_file = open('config/beta_config.JSON')
    setup = config.loadconfig.Test(config_file)

######################### Likelihood Background #########################

# This is for ploting the Likelihood:

# We define many paths:
base_path     = current_plotting_dir = 'data/'
chosen_folder = 'likelihood_explorer_' + '19-09-11-21-32-05'
data_path     = base_path+chosen_folder

# Plotting chosen_folder:
Z_grid = np.load(data_path + '/value.npy')
X_grid = np.load(data_path + '/theta.npy')
Y_grid = np.load(data_path + '/alpha.npy')
M      = np.load(data_path + '/num_paths.npy')
N      = np.load(data_path + '/interpolation_points.npy')

# To make Z positive, we use an offset for plotting
Z_grid = Z_grid-1*np.min(Z_grid) + 1e-16

fig = plt.figure(figsize=(8, 5))
ax  = plt.axes(projection='3d', elev=50, azim=-50)
ax.plot_surface(X_grid, Y_grid, Z_grid, norm=LogNorm(), rstride=1, cstride=1, edgecolor='none', alpha=.8, cmap=plt.cm.jet)
line,  = ax.plot([], [], [], 'b', label='Newton-CG', lw=2)
point, = ax.plot([], [], [], 'bo')
plt.xlabel('$\\theta$');
plt.ylabel('$\\alpha$');
ax.set_zlabel('Likelihood')

######################### Create Folder with Date #########################

warnings.simplefilter("once")

forecast_data_in = np.load(setup.data_path) # We load the data where data_path is pointing.

#forecast_data_inter=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now(); # dtM is the package datetime, and it is inside Base_plus.py.
current_time = now.strftime("%y-%m-%d-%H-%M-%S")
script_name = 'fit_SDE'
current_data_dir = script_name + current_time
os.makedirs(current_data_dir,exist_ok=True) # We create a folder with the actual date and time.

forecast_data_inter = np.swapaxes(forecast_data_in, 0,1) # This transforms a raw into column.
print('Data output will be save in ',current_data_dir)

######################### Create all the paths #########################

# IMPORTANT: The format for the data is array array: example[2][4] = number.

N  = forecast_data_inter.shape[2] - 368 # 427 - 368 = 59.
M  = setup.num_paths
dt = 1/144; # NOTICE!! Before (28/01/2020) dt = 1.
p  = forecast_data_inter[2,:setup.num_paths,:N]
p  = p[::2,:]
V  = forecast_data_inter[1,:setup.num_paths,:N] - forecast_data_inter[2,:setup.num_paths,:N]
V  = V[::2,:]
X  = forecast_data_inter[1,:setup.num_paths,:N]
X  = X[::2,:]
M  = V.shape[0]

# Now we are going to load the new and corrected data (28/01/2020):
scalar_df, array_df    = load_matlab_csv("./data/cleansed/Table_Training_Complete.csv")
Date                   = scalar_df.Date
Time                   = array_df.Time
Forecast               = array_df.Forecast
Forecast_Dot           = array_df.Forecast_Dot
Real_UTE               = array_df.Real_UTE
Real_ADME              = array_df.Real_ADME
Error                  = array_df.Error
Error_Transitions      = array_df.Error_Transitions
Error_Lamp             = array_df.Error_Lamp
Error_Lamp_Transitions = array_df.Error_Lamp_Transitions

M,N = Forecast.shape # M is the number of paths, and N the number of measurements per path.
p   = Forecast.to_numpy()
V   = Error.to_numpy()
X   = Real_ADME.to_numpy()
dt  = Time[2][1]

# We have p, V, and X. They represent all the data that will be used.

# We use the next lines to save the data we are importing. It is needed to ensure our data is 'good'.
if 1 == 0:
    toSave1 = open('data_p.pckl','wb')
    pickle.dump(p,toSave1)
    toSave1.close()
    toSave2 = open('data_X.pckl','wb')
    pickle.dump(X,toSave2)
    toSave2.close()
    toSave3 = open('data_V.pckl','wb')
    pickle.dump(V,toSave3)
    toSave3.close()
    exit()

disct_temp = disct(N,dt,M) # We want to follow V = X - P $\implies$ X = V + P.
# disct is a class from Base_plus. It just contains the number of paths, transitions, and dt.

######################### If we choose Lamperti: #########################

if setup.likelihood == 'lamperti_likelihood_SDE_approx':
    V = np.arcsin(2*X - 1)
    print(' Data has been Lamperti transformed ')

if setup.likelihood == 'lamperti_likelihood_linearized':
    V = np.arcsin(2*X - 1)
    print(' Data has been Lamperti transformed ')

this_model = model_modified_drift(disct_temp,V, forecast = p)
# model_modified_drift is a class from Base_plus.
# It just contains the number of paths, transitions, and dt (inside disct_temp).
# The data V, and the forcast p.

# print(this_model.expected_cross_time()) # We print this to test the function expected_cross_time.

# Notice that current_data_dir has the directory where most outputs should be written.
file_object = open(current_data_dir+'/results.out', 'w')
# In file_object we will write some outputs.
copyfile('config/beta_config_test.JSON', current_data_dir + '/beta_config_test.JSON' )
# We copy the JSON file into the current_data_dir directory. Then, we can always see what we used.

# width,height,angle, eig_vect_opt,Hess, FAIL = this_model.compute_ellipse(inference=setup.likelihood, param=(12, 0.4 ), batch_size=10 , plot=False);

# Remmeber that setup has inside the JSON data.
# optimization is a class in JSON that has many parameters inside (e.g., initial_batch_size).
current_batch_size = setup.optimization['initial_batch_size'];
intial_point       = np.array(( setup.optimization['theta_init'] , setup.optimization['alpha_init']  ))

print('Computing Hessian at intial point')

# width,height,angle, eig_vect_opt,Hess, FAIL = this_model.compute_ellipse(inference=setup.likelihood, param=(setup.optimization['theta_init'], setup.optimization['alpha_init']), batch_size=current_batch_size , plot=False);

# THIS IS WRONG, WE NEED TO USE THE SAME BATCH AS IN THE OPTIMIZOR:

# We use the method gen_mini_batch(self, batch_size, send_end_to_evaluater = None)
batch = this_model.gen_mini_batch(setup.optimization['initial_batch_size'], send_end_to_evaluater = None)
# We use the method likelihood_evaluate(self, selector, param, batch, batch_size = None).
# We just evaluate the likelihood in the initial parameters for some initial batch called 'batch'.
likelihood_initial_value = this_model.likelihood_evaluate(selector = setup.likelihood,
param = (setup.optimization['theta_init'], setup.optimization['alpha_init']), batch = batch,
batch_size = setup.optimization['initial_batch_size'])

#creating initial array
# parmeter_convergence=np.array((setup.optimization['theta_init'],setup.optimization['alpha_init'], likelihood_initial_value , setup.optimization['initial_batch_size'], height, width,angle ))
# hessian_convergance=Hess
# eigen_vect_convergence=eig_vect_opt

parmeter_convergence_limited = np.array((setup.optimization['theta_init'], setup.optimization['alpha_init'],
likelihood_initial_value, setup.optimization['initial_batch_size']))

# Printing:
print('Optimizing ' + setup.likelihood + ' using ' + setup.optimizer + ' with ' + \
str(setup.optimization['initial_batch_size']) + ' samples in the batch.')

# print( 'Ellipse ' + ' H: ',height , 'W: ', width )

####################################################################################
######################### We iterate over the minibatches: #########################
####################################################################################

print('We start the minibatches loop...')
while current_batch_size <= setup.optimization['max_batch_size']:

    optim = this_model.optimize(inference = setup.likelihood, batch_size = current_batch_size,
    method = setup.optimizer, param_initial = intial_point, niter = 1, temp = 0);

    print(optim.message[0])
    print(optim.x)

    result_note = 'Optimization result for a batch of size ' + str(current_batch_size) + \
    ' initialized at ' + str(intial_point) + ' is ' + str(optim.x) + ' with functional value ' + \
    str(optim.fun) + ' with message ' + optim.message[0] + ' results save in ' + current_data_dir;

    print(result_note)

    file_object.write(result_note)

    parmeter_convergence_limited = np.vstack((parmeter_convergence_limited,np.hstack((optim.x,optim.fun,current_batch_size))));
    np.save(current_data_dir + '/parmeter_convergence_limited_' + str(current_batch_size), parmeter_convergence_limited);
    print('Results save in ' + current_data_dir)

    #update batch size and intial point
    current_batch_size = current_batch_size * setup.optimization['batch_multiplier']
    intial_point       = optim.x
    start_note         = 'Startimg fit for a batch of size ' + str(current_batch_size) + ' initialized at ' + str(intial_point);
    file_object.write(start_note)
    print(start_note)

file_object.close()

print('results have been saved in ' + current_data_dir + '/results.out')

######################################################################################################
