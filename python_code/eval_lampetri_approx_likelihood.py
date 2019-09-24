import os
os.chdir('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code')
from Base_plus import *

print( ' starting approx lamperti evaluator ')

orig_stdout = sys.stdout
f = open('approx_lamperti.logs', 'w')
sys.stdout = f

def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))

warnings.showwarning = customwarn

warnings.simplefilter("once")


# sys.stdin.close()
# sys.stdin = open(os.devnull)


#from multi_path_base import model_modified_drift
forecast_data_in=np.load('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code/data/cleansed/URG_forecast_data_A_2018.npy')
now = dtM.datetime.now();
current_time=now.strftime("%y-%m-%d-%H-%M-%S")
script_name='likelihood_explorer'
current_plotting_dir='data/'+ script_name+'_' +current_time
os.mkdir(current_plotting_dir) # make a time stamped folder to contain all plots

forecast_data_inter=np.swapaxes(forecast_data_in, 0,1)

N=forecast_data_inter.shape[2];N
M=forecast_data_inter.shape[1]-973;M
dt=1
dN=1/N;dN

disct_temp = disct(N,dt,M)
p=forecast_data_inter[2,:-973,:]  #240
V=forecast_data_inter[2,:-973,:] -forecast_data_inter[1,:-973,:]


this_model=model_modified_drift(disct_temp,V, forecast= p)



theta_array=np.linspace(0.01, 25, 25)
alpha_array=np.linspace(0.01, 1, 25)

np.save(current_plotting_dir+'/num_paths', M)
np.save(current_plotting_dir+'/interpolation_points', N)
likelihood_evaluater_pool(theta_array=theta_array,\
    alpha_array=alpha_array,\
    this_model=this_model,inference="lamperti_likelihood_SDE_approx", path_dir=current_plotting_dir )

sys.stdout = orig_stdout
f.close()

print( 'output arrays are successfully saved in ', current_plotting_dir)
