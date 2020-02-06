import os
import sys
os.chdir(sys.path[0])
from Base_plus import *
import config.loadconfig

# This code is extremely messy. I will do it in MATLAB (30/01/2020).

warnings.filterwarnings('error', '.*invalid value encountered.*',)

if 1 == 1:
    parser = argparse.ArgumentParser(description = 'Likelihood Evaluator v1.0')
    parser.add_argument('-filename', help = ' Config file name or path')
    parser.add_argument('--version', action = 'version', version = 'Likelihood Evaluator v1.0')
    args = parser.parse_args()
    config_file      = open(args.filename)  # config/beta_config.JSON.
    setup            = config.loadconfig.Test(config_file)
else:
    config_file = open('config/beta_config.JSON')
    setup = config.loadconfig.Test(config_file)

print(' loading configuration ')

forecast_data_in = np.load(setup.data_path)
forecast_data_in.shape
data             = np.swapaxes(forecast_data_in, 0,1)
data_unique      = np.unique(data[:,:,:], axis=1  )
now              = dtM.datetime.now();

current_time         = now.strftime("%y-%m-%d-%H-%M-%S")
script_name          = 'sample_paths_'
current_plotting_dir = 'plots/'+ script_name +current_time

os.makedirs(current_plotting_dir, exist_ok = True) # Make a time stamped folder to contain all plots.

print(current_plotting_dir)

current_list = list(range(0, setup.num_paths)) # =list(range(0, data.shape[1])) #full_set.
M            = 10 # Number of simulated paths in each plot.
N            = 427
disct_in     = disct(N = N, dt = 1, M = M) # M is number of paths wanted in ensamble.
real_in      = real(4.3, 0.02) # SDE parameters to be generated from.

################## This now is MY modification:
# Now we are going to load the new and corrected data (30/01/2020):
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

disct_in     = disct(N = N, dt = dt, M = 10) # M is number of paths wanted in ensamble.

##################

path_simulator(forecast_data_inter = data_unique, hours = 24,\
    real_in = real_in, disct_in = disct_in, list_forecast_number = current_list,\
    dir_path = current_plotting_dir)
    # Notice: We plot from 0 transitions to "hours" transitions. E,g,. is "hours" = 100, we plot from 0 to 100
    # where this 100 has units transitions.
