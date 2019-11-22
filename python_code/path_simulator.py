import os
import sys
os.chdir(sys.path[0])
from Base_plus import *

import config.loadconfig

warnings.filterwarnings('error', '.*invalid value encountered.*',)

parser = argparse.ArgumentParser(description='Likelihood Evaluator v1.0')
parser.add_argument('-filename', help=' Config file name or path')
parser.add_argument('--version', action='version', version='Likelihood Evaluator v1.0')
args = parser.parse_args()

print( ' loading configuration ')
# args.filename
config_file= open(args.filename)  #config/beta_config.JSON

setup=config.loadconfig.Test(config_file)


forecast_data_in=np.load(setup.data_path)

forecast_data_in.shape

data=np.swapaxes(forecast_data_in, 0,1)

data_unique=np.unique(data[:,:,:], axis=1  )


now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='sample_paths_'
current_plotting_dir='plots/'+ script_name +current_time
# os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots
os.makedirs(current_plotting_dir,exist_ok=True)

print(current_plotting_dir)

current_list= list(range(0, setup.num_paths)) #=list(range(0, data.shape[1])) #full_set
M=1 #number of simulated paths in each plot
N=427
disct_in= disct(N=N, dt=1, M=M) #M is number of paths wanted in ensamble
real_in=real(1.36, 0.020) # SDE parameters to be generated from


path_simulator(forecast_data_inter=data_unique,hours=59,\
    real_in=real_in, disct_in=disct_in,list_forecast_number=current_list,\
    dir_path=current_plotting_dir)
    # Notice: We plot from 0 transitions to "hours" transitions. E,g,. is "hours" = 100, we plot from 0 to 100
    # where this 100 has units transitions.

# path_simulator(forecast_data_inter=data,hours=427,\
#     real_in=real_in, disct_in=disct_in,list_forecast_number=current_list,\
#     dir_path=current_plotting_dir)
