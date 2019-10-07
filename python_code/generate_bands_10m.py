import os
import sys
os.chdir(sys.path[0]+'/python_code')
from Base_plus import *
import config.loadconfig

warnings.filterwarnings('error', '.*invalid value encountered.*',)

parser = argparse.ArgumentParser(description='Likelihood Evaluator v1.0')
parser.add_argument('-filename', help=' Config file name or path')
parser.add_argument('--version', action='version', version='Likelihood Evaluator v1.0')
args = parser.parse_args()

print( ' loading configuration ')
#args.filename
config_file= open('config/approx_lamperti_config_test.JSON')  #config/beta_config.JSON


setup=config.loadconfig.Test(config_file)


forecast_data_in=np.load(setup.data_path)

forecast_data_in.shape

data=np.swapaxes(forecast_data_in, 0,1)


# forecast_with_data=np.load('data/forecast_with_data.npy')
# data=data_check_interpolate(forecast_with_data=forecast_with_data)
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='gen_CI_'
current_plotting_dir='plots/'+ script_name +current_time
# os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots
os.makedirs(current_plotting_dir,exist_ok=True)

print(current_plotting_dir)

current_list= list(range(0, 600)) #=list(range(0, data.shape[1])) #full_set
disct_in= disct(N=427, dt=1, M=1000) #M is number of paths wanted in ensamble
real_in=real(12,0.1) # SDE parameters to be generated from

empirical_Confidence_Interval_plots(forecast_data_inter=data,  real_in=real_in,\
            disct_in=disct_in, list_forecast_number=current_list,\
            dir_path=current_plotting_dir, hours=72)

print(current_plotting_dir)
