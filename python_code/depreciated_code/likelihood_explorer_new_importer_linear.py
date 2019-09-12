import os
os.chdir('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code')
from Base_plus import *

orig_stdout = sys.stdout
f = open('linear_lamperti_out.txt', 'w')
sys.stdout = f

def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))

warnings.showwarning = customwarn

warnings.simplefilter("once")

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


Z = np.arcsin(2*(forecast_data_inter[1,:-973,:]) - 1)
#now for V
this_model=model_modified_drift(disct_temp,Z, forecast= p) #change it

#this_model.lamperti_likelihood_linearized(np.array((15,0.1)))



theta_array=np.linspace(0.01, 25, 10)
alpha_array=np.linspace(0.01, 2, 10)


np.save(current_plotting_dir+'/num_paths', M)
np.save(current_plotting_dir+'/interpolation_points', N)
likelihood_evaluater(theta_array=theta_array,\
    alpha_array=alpha_array,\
    this_model=this_model,inference="lamperti_likelihood_linearized", path_dir=current_plotting_dir )
print(current_plotting_dir)


sys.stdout = orig_stdout
f.close()
