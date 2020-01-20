import numpy as np
import matplotlib.pyplot as plt #general plotting
import matplotlib.cm as cm #plotting color map
import matplotlib as mpl #plotting ellipse
import scipy.stats
from scipy.stats import multivariate_normal, norm  #used to evaluate normal and path generation
from mpl_toolkits.mplot3d import Axes3D #3D plots
import operator
from scipy.optimize import minimize #optimization wraper
from scipy.optimize import basinhopping #opimization wrapper
import math
import random
from scipy import interpolate #interpolation , might not be used anymore
from math import pi #constant
from numpy import linalg as LA # linear algebra operations for hessian
import datetime as dtM #processing of datetime data
import os #make directories / change working path
import multiprocessing as mp # prallelizing likelihood_evaluater
import time #time for timing tasks
from scipy.integrate import solve_ivp #RK4 numerical solving odes
import sys #to obtain system current path
import warnings #make each warning appear only once, warning redirection
from tqdm import tqdm #progress bar output in stderr
import pickle # We use this to save and load variables.
import argparse #Command line input
import itertools as iter #iter object tools for batching
#import numdifftools as nd #differentiation / Hessian computation

########################

#Warning controls

warnings.filterwarnings('error', '.*invalid value encountered.*',)

##############
#constants
Nfeval=0 # for printing .
##############
class disct: #discritization class
  def __init__(self, N,dt, M ):
    self.N = N #numeber of samples
    self.dt = dt  # dt of time , not used efficitvely now
    self.M = M #number of paths

class real: #real parameters of the model (for siumulation only)
  def __init__(self, mu,sigma ):
    self.mu = mu #mean reversion (this is theta_0 in paper)
    self.sigma = sigma #path variability (it is alpha in paper)

class optimal(real): #this is for the optimal parameters we found/fitted
  def __init__(self):
    real.__init__(self, 0,0 ) # This copies the class real but sets its values in (0,0)

class model_modified_drift: #this is the MAIN / MOST recent class for model.
    def __init__(self, disct, data,forecast ): # Initializations. You need to add the data inside.
        self.disct = disct # disc class
        self.data  = data # array with all the paths and (maybe forcasts)
        #self.ic = ic
        self.forecast = forecast # All the forcasts
        self.optimal  = optimal() # Initialization (with (0,0))
        self.Nfeval   = 0 # for printing , constant

#####################################

    def expected_cross_time(self):
        # We are assuming that all the paths have the same length.
        V=self.data
        N=self.disct.N
        M=self.disct.M
        counter = 0;
        for j in range(0,M):
            for i in range(0,N-1):

                if V[j,i]*V[j,i+1] < 0:
                    counter += 1;

        return(M*N/counter) # The units are 10 minutes, i.e., if expected_cross_time = 12.3 then the expected time is 123 minutes.

    def likelihood_evaluate(self, selector, param, batch, batch_size = None): # The none is ignored if we put a value.
        # When you optimize, param is the iteration starting point.
        # param is the theta and alpha of the Likelihood. The Likelihood is constacted by theta, alpha and data.
        # The Likelihood is constructed assuming that theta is not theta 0, which is wrong.

        if selector == "beta_likelihood":
            # print('Evaluating using the Beta likelihood')
            print(' still have to do it ')
            return()
        elif selector == 'lamperti_likelihood_SDE_approx':
            # likelihood=self.lamperti_likelihood_SDE_approx
            # print('Evaluating using lamperti_likelihood_SDE_approx')
            print(' still have to do it ')
            return()
        elif selector == 'lamperti_likelihood_linearized':
            # likelihood=self.lamperti_likelihood_linearized
            # print('Evaluating using lamperti_likelihood_linearized')
            print(' still have to do it ')
            return()
        elif selector == 'rand_beta_objective': # This is the best until NOW!!! V with no Lamparti.
            print('Evaluating using rand_beta_objective')
            return(self.rand_beta_objective(param, batch_size, batch))
        elif selector == 'rand_approx_lamperti': # Best for Lamparti!!! V with Lamparti.
            print('Evaluating using rand_approx_lamperti')
            return(self.rand_approx_lamperti(param, batch_size ) )
        else:
            print('Requested likelihood not found! ')

        return()

#####################################

    # def initial_param(self,param,send_end_to_evaluater=None): ## lamparti likelihod Number 1!
    #     #data
    #     D=self.data # Model data (it can be V, Z or even X).
    #     p=self.forecast
    #     #discretization object
    #     N=self.disct.N
    #     M=self.disct.M
    #     #input parameters
    #     theta,alpha, _ = param;

#####################################

    def lamperti_likelihood_SDE_approx(self,param,send_end_to_evaluater=None): ## lamparti likelihod Number 1!
        #data
        X=self.data
        p=self.forecast
        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M
        #input parameters
        theta,alpha, _ = param;
        L_n=0;
        L_m=0;
        eps=np.finfo(float).eps;
        dN=1/N
        dA=0
        counter=0
        for j in range(0,M):
            mean=0
            var=0
            L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            L_n=0
            for i in range(0,N-2): #start from 1 or zero #Implement Tripizoidal scheme

                mean, var =moment_compute_approx_lamperti(dN, X_prev=X[j,i], X_next=X[j,i+1], p_prev=p[j,i],p_next=p[j,i+1], theta_prev=theta_adjusted[i], theta_next=theta_adjusted[i+1],alpha=alpha )

                L_n = L_n  - ( X[j,i+1]- X[j,i] - mean   )**2/(2*var) -0.5*np.log(2*math.pi*var)

        # display information
        #if self.Nfeval%5 == 0:
        # print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        # self.Nfeval += 1
        #print('counter is ', counter)
        print('l= ', -1*L_m)
        if send_end_to_evaluater != None:
            send_end_to_evaluater.send(-1*L_m)
        return(-1*L_m)

    def gen_mini_batch(self, batch_size, send_end_to_evaluater = None): # param
        #data
        X = self.data
        p = self.forecast
        N = self.disct.N
        M = self.disct.M
        #input parameters
        # theta,alpha = param;
        batch_size = int(batch_size)

        # a batch of pairwise consecutive samples
        list_pairs_V=[]
        list_pairs_p=[]
        list_pairs_pdot=[]
        # list_pairs_theta=[]
        for i in range(0,M-1): # the pairs are (y_{i-1},y_i), all the transitions
            # N=len(p[i,1:])
            p_dot = np.diff(p[i,1:])*N # We have to check if actually p_dot(t*) = dp/dt (t*) and not, e.g., t*+1 or t*-1
            list_pairs_V.append(pairwise(X[i,1:-1]))
            list_pairs_p.append(pairwise(p[i,1:-1]))
            list_pairs_pdot.append(pairwise(p_dot))
            # list_pairs_theta.append(pairwise( theta_adjust( theta , p[i,1:] )[0])) # parameters
            # Be sure that, even when the likelihood is contructed assuing theta=theta_0, it is
            # fine to supply theta_=theta_t.
        obs_pairs      = iter.chain(*list_pairs_V); # the star * transform a list into an array
        forecast_pairs = iter.chain(*list_pairs_p);
        pdot_pairs     = iter.chain(*list_pairs_pdot);
        # theta_pairs=iter.chain(*list_pairs_theta);
        # combined_iter=iter.zip_longest(obs_pairs,forecast_pairs,theta_pairs) #fillvalue=np.nan
        combined_iter = iter.zip_longest(obs_pairs,forecast_pairs,pdot_pairs) #fillvalue=np.nan
        # combined_iter = (x_{i-1},x_1,...,,,theta_i) wirh ALL THE DATA
        batch = random_combination(combined_iter, batch_size) # you sample from the r^6 chain.

        # Testing:
        # toSave = open('allBatches.pckl','wb')
        # pickle.dump(batch,toSave)
        # toSave.close()
        # exit()

        return(batch)

    def rand_beta_objective(self, param, batch_size, batch, send_end_to_evaluater = None):

        X  = self.data
        p  = self.forecast
        N  = self.disct.N
        dt = self.disct.dt
        M  = self.disct.M
        # Input parameters:
        theta,alpha = param;
        # batch_size = int(batch_size):
        L_n     = 0;
        L_m     = 0;
        a       = 0;
        b       = 0;
        dN      = 1/N;
        counter = False;

        ##############
        num_tries = 10
        for attempt in range(num_tries): # WE TRY MANY TIMES IN CASE WE have a bad data point that produces an error.
            try:
                for j in range(batch_size):

                    # batch[j] = (obs_pairs,forecast_pairs,pdot_pairs)
                    pdot_prev  = batch[j][2][0]
                    pdot_next  = batch[j][2][1]
                    p_prev     = batch[j][1][0]
                    p_next     = batch[j][1][1]
                    theta_prev = theta_t(theta = theta, p = p_prev, pdot = pdot_prev) # theta_t(theta, p, pdot).
                    theta_next = theta_t(theta = theta, p = p_next, pdot = pdot_next)
                    m_1,m_2    = beta_moment(dN, X_prev = batch[j][0][0], X_next = batch[j][0][1],
                    p_prev = p_prev, p_next = p_next, theta_prev = theta_prev, theta_next = theta_next, alpha = alpha)

                    a = m_1; # Mean.
                    b = m_2 - m_1**2; # Variance.
                    # Mean and variance of the transitions in the SDE.

                    # Sanity checks:
                    if np.isnan(a): print('a is nan')
                    if np.isnan(b): print('b is nan')
                    if not np.isfinite(a): print('a is infinite')
                    if not np.isfinite(b): print('b is infinite')
                    if b == 0: b = 1e-16 # Remmeber b is variance!
                    if b < 0:
                        b = 1e-16
                        counter = True

                    # Shape parameters of the beta distribution:
                    beta_param_alpha = -((1+a) * (a**2+b-1)) / (2*b)
                    beta_param_beta  = ((a-1) * (a**2+b-1)) / (2*b)

                    L_n = L_n + (beta_param_alpha-1) * np.log((batch[j][0][1]+1)/2) + \
                    (beta_param_beta-1) * np.log(1-(batch[j][0][1]+1)/2) - \
                    scipy.special.betaln(beta_param_alpha, beta_param_beta)
                    # L_n is the value of the Likelihood, given some parameters and given some data batch.

                if np.isnan(L_n):
                    raise ValueError('NaN value in likelihood')
                # Display information:
                # if self.Nfeval%5 == 0:
                print('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_n))
                self.Nfeval += 1 # Nfeval just gives you the number of times the Likelihood has been called.
                if send_end_to_evaluater != None: # send_end_to_evaluater is a communication channel in case we want to parallelize.
                    send_end_to_evaluater.send(-1*L_n)
                return(-1*L_n)

            except ValueError:
                print('Trying again... trial number ' + str(attempt))

    def rand_approx_lamperti(self, param,batch_size, send_end_to_evaluater = None):
        #data
        X=self.data
        p=self.forecast
        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M
        #input parameters
        theta,alpha = param;
        batch_size=int(batch_size)
        L_n=0;
        L_m=0;
        a=0;b=0;
        dN=1/N
        counter=False

        ##############
        # a batch of pairwise consecutive samples

        if alpha <= 0 : return(np.inf)

        list_pairs_V=[]
        list_pairs_p=[]
        list_pairs_theta=[]

        for i in range(M):
            list_pairs_V.append(pairwise( X[i,1:-1]))
            list_pairs_p.append(pairwise( p[i,1:-1]))
            list_pairs_theta.append(pairwise( theta_adjust( theta , p[i,1:] )[0]))

        obs_pairs=iter.chain(*list_pairs_V);
        forecast_pairs=iter.chain(*list_pairs_p);
        theta_pairs=iter.chain(*list_pairs_theta);
        combined_iter=iter.zip_longest(obs_pairs,forecast_pairs,theta_pairs) #fillvalue=np.nan
        batch=random_combination(combined_iter,batch_size)
        ##############
        num_tries=10
        for attempt in range(num_tries):
            try:
                L_n=0
                for j in range(batch_size):

                    m_1,m_2=moment_compute_approx_lamperti( dN, X_prev=batch[j][0][0] , X_next=batch[j][0][1], p_prev=batch[j][1][0] ,p_next=batch[j][1][1], theta_prev=batch[j][2][0], theta_next=batch[j][2][1],alpha=alpha )

                    a=m_1; #mean
                    b=m_2- m_1**2; #variance

                    #sanity checks
                    if np.isnan(a): print('a is nan')
                    if np.isnan(b): print('b is nan')
                    if not np.isfinite(a): print('a is infinite')
                    if not np.isfinite(b): print('b is infinite')
                    if b==0: b=1e-16
                    if b<0:
                        b=1e-16
                        counter = True

                    #shape parameters of the beta distribution
                    beta_param_alpha= - ( (1+a)*(a**2 +b -1)  )/(2*b)
                    beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                    if np.isnan(beta_param_alpha): print('beta_param_alpha is nan')
                    if np.isnan(beta_param_beta): print('beta_param_beta is nan')

                    temp=(beta_param_alpha-1 )*np.log(  (batch[j][0][1]+1)/2 ) +\
                     (beta_param_beta-1)*np.log(1-( batch[j][0][1] +1)/2 )-\
                     scipy.special.betaln(beta_param_alpha, beta_param_beta)


                    # print('input is ',dN, batch[j][0][0] , batch[j][0][1], batch[j][1][0] ,batch[j][1][1], batch[j][2][0], batch[j][2][1],alpha  )
                    # print('current= ' , temp )
                    # print('prev_L_n', L_n)
                    # print('new sum ', L_n+ temp)


                    L_n = L_n + (beta_param_alpha-1 )*np.log(  (batch[j][0][1]+1)/2 ) +\
                     (beta_param_beta-1)*np.log(1-( batch[j][0][1] +1)/2 )-\
                     scipy.special.betaln(beta_param_alpha, beta_param_beta)

                    if L_n==np.inf:
                        # print('input was ',dN, batch[j][0][0] , batch[j][0][1], batch[j][1][0] ,batch[j][1][1], batch[j][2][0], batch[j][2][1],alpha  )
                        return( np.inf )
                        # raise ValueError('L_n is not finite')
                    if L_n==-1*np.inf:
                        # print('input was ',dN, batch[j][0][0] , batch[j][0][1], batch[j][1][0] ,batch[j][1][1], batch[j][2][0], batch[j][2][1],alpha  )
                        return( np.inf )
                        # raise ValueError('L_n is not finite')

                if np.isnan(L_n):
                    raise ValueError('NaN value in likelihood')
                # display information
                # if self.Nfeval%5 == 0:
                print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_n ))
                self.Nfeval += 1
                if L_n==np.inf: L_n=-1*np.inf
                if send_end_to_evaluater != None:
                    send_end_to_evaluater.send(-1*L_n)
                return(-1*L_n)

            except Error:
                print('trying again...trial number '+ str(attempt) )


    def compute_ellipse(self, inference , param , batch_size , plot=False):

        if inference=="beta_likelihood":
            likelihood=self.beta_likelihood
            print('Evaluating using the Beta likelihood')

        elif inference=='lamperti_likelihood_SDE_approx':
            likelihood=self.lamperti_likelihood_SDE_approx
            print('Evaluating using lamperti_likelihood_SDE_approx')

        elif inference=='lamperti_likelihood_linearized':
            likelihood=self.lamperti_likelihood_linearized
            print('Evaluating using lamperti_likelihood_linearized')

        elif inference=='rand_beta_objective':
            likelihood=self.rand_beta_objective
            print('Evaluating using rand_beta_objective')

        elif inference=='rand_approx_lamperti':
            likelihood=self.rand_approx_lamperti
            print('Evaluating using rand_approx_lamperti')

        else: print('Requested likelihood not found! ')

        wraper=lambda param: likelihood( (param[0] , param[1]), batch_size )

        for attempt in range(1):
            try:
                Hess=nd.Hessian(wraper)((param[0], param[1])); #copmute hessian

                width,height,angle, eig_vect_opt=Ellipse(Hess, param)

                if plot==True:
                    plt.clf()
                    fig, ax = plt.subplots()
                    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle = 180+angle, alpha=1, color='k')
                    V = np.array([eig_vect_opt[0],eig_vect_opt[1]])
                    origin = [10],[0.1] # origin point
                    #plt.quiver(*origin, V[:,0],V[:,1], color=['r','b'], scale=10)
                    ax.add_patch(ell)
                    ax.autoscale()
                    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
                    plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
                    file_name='ellipse'+str(M) +'_interp='+ str(N)
                    plt.title('Ellipse of the Hessian of Log-likelihood at the optimal \n ( batch size= ' + str(current_batch_size) + ' )' );
                    plt.xlabel('$\\theta$');
                    plt.ylabel('$\\alpha$');
                    plt.savefig( 'plots/'+chosen_folder+'/'+file_name + '.pdf')

                FAIL= False

                return(width,height,angle, eig_vect_opt,Hess, FAIL )

            except:
                print('Value error in obtaining Hessian, we try again' , attempt )
                FAIL= True;
        return(width,height,angle, eig_vect_opt,Hess, FAIL )


    def beta_likelihood(self,param, send_end_to_evaluater=None):
        #data
        X=self.data
        p=self.forecast
        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M
        #input parameters
        theta,alpha,_ = param;
        L_n=0;
        L_m=0;
        a=0;b=0;
        dN=1/N
        counter=False
        for j in range(0,M):
            L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            L_n=0
            for i in range(0,N-2):
                m_1,m_2=beta_moment(dN, X_prev=X[j,i], X_next=X[j,i+1], p_prev=p[j,i],p_next=p[j,i+1], theta_prev=theta_adjusted[i], theta_next=theta_adjusted[i+1],alpha=alpha )

                a=m_1; #mean
                b=m_2- m_1**2; #variance

                #sanity checks
                if np.isnan(a): print('a is nan')
                if np.isnan(b): print('b is nan')
                if not np.isfinite(a): print('a is infinite')
                if not np.isfinite(b): print('b is infinite')
                if b==0: b=1e-16
                if b<0:
                    b=1e-16
                    counter = True

                #shape parameters of the beta distribution
                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                L_n = L_n + (beta_param_alpha-1 )*np.log(  (X[j,i+1]+1)/2 ) +\
                 (beta_param_beta-1)*np.log(1-( X[j,i+1] +1)/2 )-\
                 scipy.special.betaln(beta_param_alpha, beta_param_beta)

        # display information
        #if self.Nfeval%5 == 0:
        #print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        self.Nfeval += 1
        print('l= ', -1*L_m)
        if send_end_to_evaluater != None:
            send_end_to_evaluater.send(-1*L_m)

        return(-1*L_m)

    def beta_likelihood_parallel(self,param, send_end_to_evaluater):
        #data
        X=self.data
        p=self.forecast
        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M
        #input parameters
        theta,alpha,_ = param;
        #print( theta, alpha)
        L_n=0;
        L_m=0;
        a=0;b=0;
        eps=np.finfo(float).eps;
        m_1=0
        m_2=0
        dN=1/N
        counter=False #compute_path_moments
        jobs = []
        pipe_list = []
        for j in range(0,M):
            #L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            counter=False
            recv_end, send_end =  mp.Pipe(False)
            arg_list=[j,N, X, p, theta_adjusted, alpha]
            p_temp = mp.Process(target=compute_path_moments, args=(arg_list,send_end))
            jobs.append(p_temp)
            pipe_list.append(recv_end)
            p_temp.start()

        for proc in jobs:
            proc.join()
        result_list = [x.recv() for x in pipe_list]
        #print(result_list)
        L_m=np.sum(result_list)

        # display information
        #if self.Nfeval%5 == 0:
        #print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        self.Nfeval += 1
        #print(counter)
        send_end_to_evaluater.send(-1*L_m)
        return(-1*L_m) #, m_1_list,m_2_list, m_1_E_list, m_2_E_list


    def lamperti_likelihood_linearized(self,param,send_end_to_evaluater=None):
        #data
        X=self.data
        p=self.forecast
        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M
        #input parameters
        theta,alpha,_ = param;
        L_n=0;
        L_m=0;
        eps=np.finfo(float).eps;
        dN=1/N
        dA=0
        counter=0
        for j in range(0,M):
            mean=0
            var=0
            L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            #print( 'shape ' , theta_adjusted.shape)
            L_n=0
            for i in range(0,N-3): #start from 1 or zero #Implement Tripizoidal scheme

                mean,var= linear_lamperti_moment(dN, X_prev=X[j,i], X_next=X[j,i+1], p_prev=p[j,i],p_next=p[j,i+1], theta_prev=theta_adjusted[i], theta_next=theta_adjusted[i+1],alpha=alpha )
                #theta_adjusted[i]
                # mean = X[j,i]*np.exp((alpha-1)*dN*theta_adjusted[i])\
                #     - dN*theta_adjusted[i]*(1-2*p[j,i])
                # var = np.exp(dN*(alpha-1)*theta_adjusted[i])*dN*np.sqrt(2*alpha*theta_adjusted[i])

                # mean = X[j,i]*np.exp((alpha-1)*(theta_adjusted[i+1] + theta_adjusted[i])*dN/2) - dN*(theta_adjusted[i+1]*(1-2*p[j,i+1])*np.exp((alpha-1)*(theta_adjusted[i+1] + theta_adjusted[i])*dN/2) + theta_adjusted[i]*(1-2*p[j,i]) )/2
                #
                # var= np.sqrt(2*alpha)*dN*(np.sqrt(theta_adjusted[i+1])*np.exp((alpha-1)*(theta_adjusted[i+1] + theta_adjusted[i])*dN/2)  + np.sqrt(theta_adjusted[i]) )/2

                #if np.isnan(mean): print('mean is nan' , mean, k)
                #if np.isnan(var): print('var is nan', var)

                L_n = L_n  - ( X[j,i+1]- X[j,i] - mean   )**2/(2*var) -0.5*np.log(2*math.pi*var)

        # display information
        #if self.Nfeval%5 == 0:
        # print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        # self.Nfeval += 1
        #print('counter is ', counter)
        print('l= ', -1*L_m)
        if send_end_to_evaluater != None:
            send_end_to_evaluater.send(-1*L_m)
        return(-1*L_m)

                                                #ADJUST UPPER BOUND !
    def optimize(self, param_initial, batch_size, inference, method = "Nelder-Mead", niter = 1 , temp = 0, bnds = ((1e-3, None), (1e-3, None))):
        print("Starting method 'opimize'...")
        # mu_initial, sig_initial=param_initial.

        # Just to print:
        if inference == "beta_likelihood":
            likelihood = self.beta_likelihood
            print('Evaluating using the Beta likelihood')
        elif inference == 'lamperti_likelihood_SDE_approx':
            likelihood = self.lamperti_likelihood_SDE_approx
            print('Evaluating using lamperti_likelihood_SDE_approx')
        elif inference == 'lamperti_likelihood_linearized':
            likelihood = self.lamperti_likelihood_linearized
            print('Evaluating using lamperti_likelihood_linearized')
        elif inference == 'rand_beta_objective':
            likelihood = self.rand_beta_objective
            print('Evaluating using rand_beta_objective')
        elif inference == 'rand_approx_lamperti':
            likelihood = self.rand_approx_lamperti
            print('Evaluating using rand_approx_lamperti')
        else: print('Requested likelihood not found!')

        if method == "Nelder-Mead":
            # gen_mini_batch(self, batch_size, send_end_to_evaluater = None).
            batch = self.gen_mini_batch(batch_size) # We create the batches we will use.
            minimizer_kwargs = {"method": "Nelder-Mead", "options":{'disp': True}}
            # "options": {options={'disp': True} 'gtol': myfactr, 'ftol': myfactr, "maxiter": 10}
            # "options": {'xatol': 10e-2 , 'fatol' : 10e-2}
            min_param = scipy.optimize.minimize(fun = likelihood, x0 = param_initial,
            args = (batch_size, batch), method = "Nelder-Mead")
            #min_param=scipy.optimize.minimize(fun=lambda param: likelihood(param, batch_size, batch) , x0=param_initial ,args=minimizer_kwargs)
            #min_param = basinhopping(lambda param: likelihood(param, batch_size, batch) , param_initial,T=temp ,minimizer_kwargs=minimizer_kwargs,niter=niter)
        if method == "L-BFGS-B":
            minimizer_kwargs = {"method": "L-BFGS-B", "options":{ 'gtol': myfactr, 'ftol': myfactr, "maxiter": 10}}
            min_param        = basinhopping(likelihood, param_initial, T = temp, minimizer_kwargs = minimizer_kwargs, niter = niter)
        if method == "SLSQP":
            minimizer_kwargs = {"method": "SLSQP", "bounds": bnds, "options":{'ftol': 1e-10, 'gtol': 1e-10}}
            min_param        = basinhopping(likelihood, param_initial, T = temp, minimizer_kwargs = minimizer_kwargs, niter = niter)

        #, options={ 'gtol': myfactr * np.finfo(float).eps, 'ftol' : myfactr * np.finfo(float).eps });min_param
        #x,f,d=scipy.optimize.fmin_l_bfgs_b(self.likelihood,param_initial, approx_grad=True)

        #x,f,d=scipy.optimize.fmin_l_bfgs_b(self.likelihood,param_initial,bounds=bnds, \
        #    factr=10, pgtol=1e-30, epsilon=1e-30, iprint=100, maxfun=1e6, maxiter=1e6,maxls=25, approx_grad=True)

        return(min_param)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)

class model_modified_drift_p_term:
    def __init__(self, disct, data,forecast ):
        self.disct = disct
        self.data = data
        #self.ic = ic
        self.forecast = forecast
        self.optimal=optimal()
        self.Nfeval=0


    def likelihood(self,param): # OLD VERSION.
        #data
        X=self.data
        p=self.forecast
        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M
        #input parameters
        theta,alpha = param;
        #print( theta, alpha)
        L_n=0;
        L_m=0;
        a=0;b=0;
        eps=np.finfo(float).eps;
        m_1=0
        m_2=0
        dN=1/N
        for j in range(0,M):
            m_1=0
            m_2=0
            L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            for i in range(0,N-1): #start from 1 or zero
                m_1=X[j,i]*np.exp(- dN*theta_adjusted[i])

                m_2=  (X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta_adjusted[i]*p[j,i]*(1-p[j,i])*(1-2*p[j,i] )) + \
                            alpha*theta_adjusted[i]*p[j,i]**2*(1-p[j,i])**2) ) /(  1+ dN*2*(theta_adjusted[i]+alpha*theta_adjusted[i]*p[j,i]*(1-p[j,i]) ))
                a=m_1
                b=m_2 - m_1**2

                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                L_n = L_n +  (beta_param_alpha-1 )*np.log(  (X[j,i+1]+1)/2 ) +\
                 (beta_param_beta-1)*np.log(1-( X[j,i+1] +1)/2 )-\
                 scipy.special.betaln(beta_param_alpha, beta_param_beta)
        # display information
        #if self.Nfeval%5 == 0:
        print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        self.Nfeval += 1

        return(-1*L_m)
###########################################################################
# Path generators / simulators

def gen_path(X0,disct,real):

    N=disct.N
    dt=disct.dt
    M=disct.M

    X = np.zeros((M,N))

    #p=0.5 + 0.5*np.sin(np.linspace(0,2*math.pi,N))

    dW=np.sqrt(dt)*np.random.normal(0, 1, N)


    for j in range(0,M):

        dW=np.sqrt(dt)*np.random.normal(0, 1, N)

        for i in range(1,N):
            X[j,0]=X0
            b= real.sigma     #sig_real*X[i-1]*(1-X[i-1])
            a= real.mu       #mu_real*(0.5 - X[i-1])
            X[j,i]= X[j,i-1]+ a*dt+ b*dW[i-1]

    return(X)

def gen_path_in_box(X0,disct,real):

    N=disct.N
    dt=disct.dt
    M=disct.M

    X = np.zeros((M,N))

    #p=0.5 + 0.5*np.sin(np.linspace(0,2*math.pi,N))

    dW=np.sqrt(dt)*np.random.normal(0, 1, N)


    for j in range(0,M):

        dW=np.sqrt(dt)*np.random.normal(0, 1, N)

        for i in range(1,N):
            X[j,0]=X0
            b= real.sigma*X[j,i-1]*(1-X[j,i-1])
            a= real.mu*(0.5 - X[j,i-1])
            X[j,i]= X[j,i-1]+ a*dt+ b*dW[i-1]

    return(X)

def gen_path_in_box_sine(X0,disct,real):

    N=disct.N
    dt=disct.dt
    M=disct.M

    X = np.zeros((M,N))

    #p=0.5 + 0.5*np.sin(np.linspace(0,2*math.pi,N))

    dW=np.sqrt(dt)*np.random.normal(0, 1, N)

    for j in range(0,M):

        dW=np.sqrt(dt)*np.random.normal(0, 1, N)

        for i in range(1,N):
            X[j,0]=X0
            b= real.sigma*X[j,i-1]*(1-X[j,i-1])
            a= real.mu*( 0.5 + 0.5*np.sin(i/N * 2*math.pi) - X[j,i-1])
            X[j,i]= X[j,i-1]+ a*dt+ b*dW[i-1]

    return(X)

def gen_path_beta_moments(X0,disct,real):

    N=disct.N
    dt=disct.dt
    M=disct.M

    theta=real.mu
    alpha=real.sigma

    X = np.zeros((M,N))

    dW=np.random.normal(0, 1, N) #removed sqrt(dt)

    for j in range(0,M):

        for i in range(1,N):
            X[j,0]=X0

            b =   alpha/(16*(1+ alpha/4)) * (1- np.exp( -2*dt*theta*(1+ alpha/4)  ))
            a =  X[j,i-1]*np.exp(- dt*theta)

            beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)   #lambd*( (a+1)  /2  ) #( (1-a)/b - 1/a )*a**2
            beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)  #lambd*( (1-a) / 2)

            X[j,i]=np.random.beta(beta_param_alpha,beta_param_beta,1)
            X[j,i]= -1 +  2*X[j,i] #a + removed

    return(X)

def gen_path_model_moments(X0,disct,real):

    N=disct.N
    dt=disct.dt
    M=disct.M

    theta=real.mu
    alpha=real.sigma

    X = np.zeros((M,N))

    dW=np.random.normal(0, 1, N) #removed sqrt(dt)

    for j in range(0,M):

        dW=np.random.normal(0, 1, N) #removed dt

        for i in range(1,N):
            X[j,0]=X0

            a =  X[j,i-1]*np.exp(- dt*theta)
            b =   alpha/(16*(1+ alpha/4)) * (1- np.exp( -2*dt*theta*(1+ alpha/4)  ))

            X[j,i]= a +  b*dW[i-1] #a + removed

    return(X)

def gen_path_model_moments_check(X0,disct,real):

    N=disct.N
    dt=disct.dt
    M=disct.M

    theta=real.mu
    alpha=real.sigma

    X = np.zeros((M,N))
    X0=0
    dW=np.random.normal(0, 1, N) #removed sqrt(dt)

    for j in range(0,M):

        dW=dt*np.random.normal(0, 1, N) #removed sqrt(dt)

        for i in range(1,N):
            X[j,0]=X0

            b=np.sqrt( 2*alpha*theta*(1/4)*(X[j,i-1] + 1/2 )*(1/2 - X[j,i-1] )  )

            a =  - X[j,i-1]*theta

            X[j,i]= X[j,i-1] + a*dt +  b*dW[i-1] #check sign of a

    return(X)

def gen_path_beta_flex(X0,disct,real,forecast):
    # TODO: We have to choose the correct normalization constant for the final time (choose correctly T).
    # Now we are using T = 1, which is not correct because we use dN = 1/N and we want that
    # dN = 1 transition.
    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    dN=1/N
    X=np.zeros((M,N))

    #while( (beta_param_alpha <=0 ) or ( beta_param_beta <= 0) ):

    for j in range(0,M):
        m_1=0
        m_2=0
        X[j,0]=X0
        m_1=X[j,0]
        m_2=X[j,0]**2
        #while( (beta_param_alpha <=0 ) or ( beta_param_beta <= 0) ):
        for i in range(0,N-1):

        #m_1= X[j,i]*(1- theta*dN)
            m_1=X[j,i]*np.exp(- dN*theta)

            m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                            X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                                alpha*theta*p[i]**2*(1-p[i])**2)

            #m_2=  (4*X[j,i]**2 -X[j,i-1]**2  + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
            #                X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
            #                    alpha*theta*p[i]**2*(1-p[i])**2) )/3


            #print('m1=',m_1)
            #print('m2=',m_2)

            a=m_1
            b=m_2 - m_1**2

            #print('mean=',a)
            #print('variance=',b)

            beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
            beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

            if beta_param_alpha<=0 or beta_param_beta<=0:
                raise Exception('Beta shape parameters are not positive definite.\
                 The parameters were: {}'.format(beta_param_alpha) +'and: {}'.format(beta_param_beta) )

            #print(beta_param_alpha, beta_param_beta)

            X[j,i+1]=np.random.beta(beta_param_alpha,beta_param_beta,1)

            X[j,i+1]= -1 +  2*X[j,i+1]

    return(X)

def gen_path_beta_robust(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    dN=1/N
    X=np.zeros((M,N))

    max_tries=50 #max(int(M*0.1), 2 )
    try_count = 0
    j=0;
    while j <M:
        m_1=0
        m_2=0
        X[j,0]=X0
        m_1=X[j,0]
        m_2=X[j,0]**2
        try:
             for i in range(0,N-1):
                if p[i]+X[j,i] >1:
                    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                if p[i]+X[j,i] <0:
                    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                m_1=X[j,i]*np.exp(- dN*theta)
                #
                # m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                #                 X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #                     alpha*theta*p[i]**2*(1-p[i])**2)

                m_2=  (X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                            alpha*theta*p[i]**2*(1-p[i])**2) ) /(  1+ dN*2*(theta+alpha*theta*p[i]*(1-p[i]) ))
                a=m_1
                b=m_2 - m_1**2

                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                if beta_param_alpha<=0 or beta_param_beta<=0:
                    raise Exception('Beta shape parameters are not positive definite.\
                     The parameters were: {}'.format(beta_param_alpha) +'and: {}'.format(beta_param_beta) )

                X[j,i+1]=np.random.beta(beta_param_alpha,beta_param_beta,1)

                X[j,i+1]= -1 +  2*X[j,i+1]
        except:
            try_count += 1
            print('Try # : ' +str(try_count))
            #print('Beta shape parameters are not positive definite.\
            # The parameters were: {}'.format(a) +'and: {}'.format(b))
            i=0
            pass
            if try_count >= max_tries:
                raise Exception("Unable to generate after %s tries" % max_tries)
        else:
            j += 1 # increments only if no exception

    return(X)

def gen_path_normal_robust(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    dN=1/N
    X=np.zeros((M,N))

    max_tries=50 #max(int(M*0.1), 2 )
    try_count = 0
    j=0;
    while j <M:
        m_1=0
        m_2=0
        X[j,0]=X0
        m_1=X[j,0]
        m_2=X[j,0]**2
        try:
             for i in range(0,N-1):
                if p[i]+X[j,i] >1:
                    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                if p[i]+X[j,i] <0:
                    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                m_1=X[j,i]*np.exp(- dN*theta)
                #
                # m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                #                 X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #                     alpha*theta*p[i]**2*(1-p[i])**2)

                m_2=  (X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                            alpha*theta*p[i]**2*(1-p[i])**2) ) /(  1+ dN*2*(theta+alpha*theta*p[i]*(1-p[i]) ))
                a=m_1
                b=m_2 - m_1**2
                #if (m_2 - m_1**2)<0:
                # if p[i]+X[j,i] >1:
                #     print(m_1,m_2 - m_1**2,p[i]+X[j,i])
                # if p[i]+X[j,i] <0:
                #     print(m_1,m_2 - m_1**2,p[i]+X[j,i])


                X[j,i+1]=np.random.normal(a,np.sqrt(b),1)


                #X[j,i+1]= -1 +  2*X[j,i+1]
        except:
            try_count += 1
            print('Try # : ' +str(try_count))
            #print('Beta shape parameters are not positive definite.\
            # The parameters were: {}'.format(a) +'and: {}'.format(b))
            i=0
            pass
            if try_count >= max_tries:
                raise Exception("Unable to generate after %s tries" % max_tries)
        else:
            j += 1 # increments only if no exception

    return(X)

def gen_path_normal_Euler(X0,disct,real,forecast):
    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma


    dN=1/N
    X=np.zeros((M,N))

    try_count = 0
    j=0;
    while j <M:
        X[j,0]=X0
        a=0
        b=0
        dW=np.sqrt(dN)*np.random.normal(0, 1, N)

        for i in range(0,N-1):
            X[j,0]=X0
            b=  2*alpha*theta*p[i]*(1-p[i])*(X[j,i]+ p[i])*(1-X[j,i]- p[i])
            a= - theta*X[j,i]

            # if X[j,i] >1:
            #     print('WARNING: X greater than one ')
            #     print(a,b,X[j,i])
            # if X[j,i] <0:
            #     print('WARNING: X less than one ')
            #     print(a,b,X[j,i])

            X[j,i+1]= X[j,i] + a*dN+ np.sqrt(b)*dW[i]

            #a= -theta*(X[j,i] - p[i])
            #b= 2*alpha*theta*p[i]*(1-p[i])*X[j,i]*( 1-X[j,i] )

            #if (m_2 - m_1**2)<0:


            #X[j,i+1]+=np.random.normal(a*dN,np.sqrt(b*dN),1)

        j += 1 # increments only if no exception

    return(X)

def gen_X_beta(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast
    dp=np.gradient(p) #derivative tracking

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    dN=1/N
    X=np.zeros((M,N))

    max_tries=50 #max(int(M*0.1), 2 )
    try_count = 0
    j=0;
    while j <M:
        m_1=0
        m_2=0
        X[j,0]=X0
        m_1=X[j,0]
        m_2=X[j,0]**2
        try:
             for i in range(0,N-1):
                #if p[i]+X[j,i] >1:
                #    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                #if p[i]+X[j,i] <0:
                #    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                m_1=( X[j,i] + dN*(theta*p[i+1] ) )/( 1+ dN*theta )   #X[j,i]*np.exp(- dN*theta)
                #
                # m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                #                 X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #                     alpha*theta*p[i]**2*(1-p[i])**2)

                m_2= ( X[j,i]**2 + 2*dN*(m_1*(-theta*m_1 \
                +alpha*theta*p[i+1]*(1-p[i+1])+ theta*p[i+1])))/( 1 + 2*dN*( alpha*theta*p[i+1]*(1-p[i+1]) ) )

                #(X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #            alpha*theta*p[i]**2*(1-p[i])**2) ) /(  1+ dN*2*(theta+alpha*theta*p[i]*(1-p[i]) ))

                a=m_1
                b=m_2 - m_1**2
                #print(X[j,i])
                #if (m_2 - m_1**2)<0:
                # if X[j,i] >1:
                #     print('WARNING: X greater than one ')
                #     print(m_1,m_2 - m_1**2,X[j,i])
                # if X[j,i] <0:
                #     print('WARNING: X less than one ')
                #     print(m_1,m_2 - m_1**2,X[j,i])

                beta_param_alpha= ( (1-a)/b - 1/a   )*a**2
                beta_param_beta= beta_param_alpha*( 1/a -1 )

                if beta_param_alpha<=0 or beta_param_beta<=0:
                    raise Exception('Beta shape parameters are not positive definite.\
                    /n The parameters were: {}'.format(beta_param_alpha) +'and: {}'.format(beta_param_beta) )

                X[j,i+1]=np.random.beta(beta_param_alpha,beta_param_beta,1)

                #X[j,i+1]= -1 +  2*X[j,i+1]
        except:
            try_count += 1
            print('Try # : ' +str(try_count))
            print('Beta shape parameters are not positive definite.\
            /n The parameters were: {}'.format(beta_param_alpha) +'and: {}'.format(beta_param_beta))
            i=0
            pass
            if try_count >= max_tries:
                raise Exception("Unable to generate after %s tries" % max_tries)
        else:
            j += 1 # increments only if no exception

    return(X)

def gen_X_normal(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast
    dp=np.gradient(p) #derivative tracking

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    dN=1/N
    X=np.zeros((M,N))

    max_tries=50 #max(int(M*0.1), 2 )
    try_count = 0
    j=0;
    while j <M:
        m_1=0
        m_2=0
        X[j,0]=X0
        m_1=X[j,0]
        m_2=X[j,0]**2
        try:
             for i in range(0,N-1):
                #if p[i]+X[j,i] >1:
                #    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                #if p[i]+X[j,i] <0:
                #    X[j,i]=p[i-1]-p[i]+X[j,i-1]

                m_1=( X[j,i] + dN*(theta*p[i+1] )  )/( 1+ dN*theta )   #X[j,i]*np.exp(- dN*theta)


                # m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                #                 X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #                     alpha*theta*p[i]**2*(1-p[i])**2)

                m_2= ( X[j,i]**2 + 2*dN*(m_1*(-theta*m_1  \
                +alpha*theta*p[i+1]*(1-p[i+1])+ theta*p[i+1])))/( 1 + 2*dN*( alpha*theta*p[i+1]*(1-p[i+1]) ) )

                #(X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #            alpha*theta*p[i]**2*(1-p[i])**2) ) /(  1+ dN*2*(theta+alpha*theta*p[i]*(1-p[i]) ))

                a=m_1
                b=m_2 - m_1**2
                #print(X[j,i])
                #if (m_2 - m_1**2)<0:
                # if X[j,i] >1:
                #     print('WARNING: X greater than one ')
                #     print(m_1,m_2 - m_1**2,X[j,i])
                # if X[j,i] <0:
                #     print('WARNING: X less than one ')
                #     print(m_1,m_2 - m_1**2,X[j,i])


                X[j,i+1]=X[j,i+1]= a+ np.sqrt(b)*np.random.normal(0,1,1)

                #X[j,i+1]= -1 +  2*X[j,i+1]
        except:
            try_count += 1
            i=0
            pass
            if try_count >= max_tries:
                raise Exception("Unable to generate after %s tries" % max_tries)
        else:
            j += 1 # increments only if no exception

    return(X)

def gen_X_normal_euler(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma


    dN=1/N
    X=np.zeros((M,N))

    try_count = 0
    j=0;
    while j <M:
        X[j,0]=X0
        a=0
        b=0
        dW=np.sqrt(dN)*np.random.normal(0, 1, N)

        for i in range(0,N-1):
            X[j,0]=X0
            b=  2*alpha*theta*p[i]*(1-p[i])*X[j,i]*(1-X[j,i])
            a= - theta*(X[j,i]-p[i])

            # if X[j,i] >1:
            #     print('WARNING: X greater than one ')
            #     print(a,b,X[j,i])
            # if X[j,i] <0:
            #     print('WARNING: X less than one ')
            #     print(a,b,X[j,i])

            X[j,i+1]= X[j,i] + a*dN+ np.sqrt(b)*dW[i]

            #a= -theta*(X[j,i] - p[i])
            #b= 2*alpha*theta*p[i]*(1-p[i])*X[j,i]*( 1-X[j,i] )

            #if (m_2 - m_1**2)<0:


            #X[j,i+1]+=np.random.normal(a*dN,np.sqrt(b*dN),1)

        j += 1 # increments only if no exception

    return(X)

def gen_X_normal_euler_derivative_tracking(X0,disct,real,forecast):

    p=forecast
    N=disct.N
    dt=disct.dt
    M=disct.M
    X0=p[0]

    theta=real.mu
    alpha=real.sigma
    issue_counter=0
    X=np.zeros((M,N))
    max_tries=100
    dN=1/N
    i=0
    j=0
    while j < M:
        X[j,0]=X0
        a=0
        b=0
        dW=np.sqrt(dN)*np.random.normal(0, 1, N)
        i=0
        while i<N-1:
            try:
                try_count=0
                X[j,0]=X0

                b=  2*alpha*theta*p[i]*(1-p[i])*X[j,i]*(1-X[j,i])
                a= - theta*(X[j,i]-p[i]) + (p[i+1]-p[i])/dN

                X[j,i+1]= X[j,i] + a*dN+ np.sqrt(b)*dW[i]
                #print(X[j,i+1])

                    #remove later
                if X[j,i+1] >1:
                    X[j,i+1]=X[j,i]

                if X[j,i+1] <0:
                    X[j,i+1]=X[j,i]

                if (X[j,i+1]<0) or (X[j,i+1])>1:
                    #print('EXCEPTION !! ')
                    Z0=X[j,i+1]
                    raise Exception('outside: {}'.format(X[j,i+1]) )
                #print('here')
            except:
                #print('issue with ', X[j,i+1])
                interpolation_count=0
                Z_init=X[j,i]
                Z_sol=X[j,i+1]
                while Z_sol<0 or Z_sol>1:
                    #print('Solving Expection')
                    try_count += 1
                    interpolation_count+=1

                    #print('Try # : ' +str(try_count))

                    interpolation_points=10*interpolation_count
                    x = np.array((0,1))
                    y = np.array(( p[i-1] , p[i] ))
                    f = interpolate.interp1d(x, y)
                    xnew=np.linspace(0,1,interpolation_points)
                    p_inside = f(xnew)

                    dN_in=1/(N+ interpolation_points)
                    dW_inner=np.sqrt(dN_in)*np.random.normal(0, 1, interpolation_points)

                    Z=np.zeros((interpolation_points))
                    Z[0]=Z_init
                    for  k in range(0,interpolation_points-1):
                        b=  2*alpha*theta*p_inside[k]*(1-p_inside[k])*Z[k]*(1-Z[k])
                        a= - theta*(Z[k]-p_inside[k]) + (p_inside[k+1]-p_inside[k])/dN_in
                        #print(b)
                        Z[k+1]= Z[k] + a*dN_in+ np.sqrt(b)*dW_inner[k]

                        #print(Z[k+1])
                    if np.isnan(Z[k+1])==False:
                        Z_sol = Z[k+1]
                        #print('achieved ',Z_sol)
                pass
                X[j,i+1]= Z_sol
                #print('achieved ',Z_sol)
                issue_counter+=1
                i+=1
                if try_count >= max_tries:
                    raise Exception("Unable to generate after %s tries" % max_tries)
            else:
                i+=1
        j+=1
    print(issue_counter, ' issues resolved')

    return(X)

def gen_X_normal_euler_DT_modified(X0,disct,real,forecast):
    # TITLE: GENERATE X-SDE from Normal with Euler using Derivative Trancking.
    p=forecast
    N=disct.N
    M=disct.M
    X0=p[0]
    theta=real.mu
    alpha=real.sigma
    X=np.zeros((M,N))
    #X_zero_drift=np.zeros(N)
    dN=1/N
    i=0
    j=0
    while j < M:
        X[j,0]=X0
        a=0
        b=0
        dW=np.sqrt(dN)*np.random.normal(0, 1, N)
        i=0
        NG_var=False
        while i<N-1:

            b=  2*alpha*theta[i]*X[j,i]*(1-X[j,i]) #remove p(1-p)
            a= - theta[i]*(X[j,i]-p[i]) + (p[i+1]-p[i])/(dN)
            if b<0:
                NG_var=True

            X[j,i+1]= X[j,i] + a*dN+ np.sqrt(b)*dW[i]

            if X[j,i+1] >1:
                X[j,i+1]=X[j,i];
            if X[j,i+1] <0:
                X[j,i+1]=X[j,i];
            #X_zero_drift[i]= p[i]+ (p[i+1]-p[i])/(dN*theta[i]);

            # TODO: Generate path in Lapmarti space and then to transform back to X.

            i+=1
        if (NG_var==True): print('negative variance in Generator. ' );
        j+=1

    return( X )

def gen_X_beta_derivative_tracking(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    dN=1/N
    X=np.zeros((M,N))
    max_tries=50 #max(int(M*0.1), 2 )
    try_count = 0
    j=0;
    while j <M:
        m_1=0
        m_2=0
        X[j,0]=X0
        m_1=X[j,0]
        m_2=X[j,0]**2
        try:
             for i in range(0,N-1):

                if X[j,i] >1:
                    X[j,i]=X[j,i-1]

                if X[j,i] <0:
                    X[j,i]=X[j,i-1]

                m_1=( X[j,i] + dN*(theta*p[i+1] ) + (p[i+1]-p[i]  )  )/( 1+ dN*theta )   #X[j,i]*np.exp(- dN*theta)


                # m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                #                 X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #                     alpha*theta*p[i]**2*(1-p[i])**2)

                m_2= ( X[j,i]**2 + 2*dN*(m_1*(-theta*m_1 + (p[i+1]-p[i])/(dN)  \
                    +alpha*theta*p[i+1]*(1-p[i+1])+ theta*p[i+1])))/( 1 + 2*dN*( alpha*theta*p[i+1]*(1-p[i+1]) ) )

                a=m_1
                b=m_2 - m_1**2
                #print(X[j,i])
                #if (m_2 - m_1**2)<0:
                # if X[j,i] >1:
                #     print('WARNING: X greater than one ')
                #     print(m_1,m_2 - m_1**2,X[j,i])
                # if X[j,i] <0:
                #     print('WARNING: X less than one ')
                #     print(m_1,m_2 - m_1**2,X[j,i])

                beta_param_alpha= ( (1-a)/b - 1/a   )*a**2
                beta_param_beta= beta_param_alpha*( 1/a -1 )

                if beta_param_alpha<=0 or beta_param_beta<=0:
                    raise Exception('Beta shape parameters are not positive definite.\
                     The parameters were: {}'.format(p[i+1] - p[i]) +'and: {}'.format(p[i+1] - p[i]) )

                X[j,i+1]=np.random.beta(beta_param_alpha,beta_param_beta,1)

                #X[j,i+1]= -1 +  2*X[j,i+1]
        except:
            try_count += 1
            print('Try # : ' +str(try_count))
            print('Beta shape parameters are not positive definite.\
            /n The parameters were: {}'.format(p[i+1] - p[i]) +'and: {}'.format(X[j,i]))
            i=0
            pass
            if try_count >= max_tries:
                raise Exception("Unable to generate after %s tries" % max_tries)
        else:
            j += 1 # increments only if no exception

    return(X)

def gen_X_normal_derivative_tracking(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    dN=1/N
    X=np.zeros((M,N))

    max_tries=50 #max(int(M*0.1), 2 )
    try_count = 0
    j=0;
    while j <M:
        m_1=0
        m_2=0
        X[j,0]=X0
        m_1=X[j,0]
        m_2=X[j,0]**2

        for i in range(0,N-1):
            #if p[i]+X[j,i] >1:
            #    X[j,i]=p[i-1]-p[i]+X[j,i-1]

            #if p[i]+X[j,i] <0:
            #    X[j,i]=p[i-1]-p[i]+X[j,i-1]

            if X[j,i] >1:
                X[j,i]=X[j,i-1]

            if X[j,i] <0:
                X[j,i]=X[j,i-1]

            m_1=( X[j,i] + dN*(theta*p[i+1] ) + (p[i+1]-p[i]  )  )/( 1+ dN*theta )   #X[j,i]*np.exp(- dN*theta)


            # m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
            #                 X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
            #                     alpha*theta*p[i]**2*(1-p[i])**2)

            m_2= ( X[j,i]**2 + 2*dN*(m_1*(-theta*m_1 + (p[i+1]-p[i])/(dN)  \
                +alpha*theta*p[i+1]*(1-p[i+1])+ theta*p[i+1])))/( 1 + 2*dN*( alpha*theta*p[i+1]*(1-p[i+1]) ) )

            a=m_1
            b=m_2 - m_1**2

            #print(X[j,i])
            #if (m_2 - m_1**2)<0:
            # if X[j,i] >1:
            #     print('WARNING: X greater than one ')
            #     print(m_1,m_2 - m_1**2,X[j,i])
            # if X[j,i] <0:
            #     print('WARNING: X less than zero ')
            #     print(m_1,m_2 - m_1**2,X[j,i])


            X[j,i+1]= a+ np.sqrt(b)*np.random.normal(0,1,1)

                #X[j,i+1]= -1 +  2*X[j,i+1]
        j += 1 # increments only if no exception

    return(X)

def gen_path_beta_robust_interpolate(X0,disct,real,forecast):

    N=disct.N
    dt=disct.dt
    M=disct.M

    p=forecast

    theta=real.mu
    alpha=real.sigma

    beta_param_alpha=0
    beta_param_beta=0

    m_1=0
    m_2=0
    #dN=1/N
    interpolation_points=1000 #add to inputs
    dx=1/interpolation_points
    X=np.zeros((M,interpolation_points))

    max_tries=500 #max(int(M*0.1), 2 ) #lower it !
    try_count = 0
    j=0;
    while j<M:

        X[j,0]=X0

        x = np.arange(1,N+1,1)

        y = forecast[j,:] #X[j,:]
        tck = interpolate.splrep(x, y, s=0)

        #y_data = X[j,:]
        #tck_data = interpolate.splrep(x, y_data, s=0)

        #xnew = np.arange(1,72+1,dx) # 1 or 72 ?
        xnew=np.linspace(1,72,interpolation_points)
        p = interpolate.splev(xnew, tck, der=0)
        #d_inter=interpolate.splev(xnew, tck_data, der=0)

        try:
            for i in range(0,interpolation_points):
                m_1=X[j,i]*np.exp(- dx*theta)

                m_2=  X[j,i]**2 + dx*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                                X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                                    alpha*theta*p[i]**2*(1-p[i])**2)
                a=m_1
                b=m_2 - m_1**2


                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                if beta_param_alpha<=0 or beta_param_beta<=0:
                    raise Exception('Beta shape parameters are not positive definite.\
                     The parameters were: {}'.format(beta_param_alpha) +'and: {}'.format(beta_param_beta) )

                X[j,i+1]=np.random.beta(beta_param_alpha,beta_param_beta,1)

                X[j,i+1]= -1 +  2*X[j,i+1]
        except:
            try_count += 1
            print('Try # : ' +str(try_count))
            i=0
            pass
            if try_count >= max_tries:
                raise Exception("Unable to generate after %s tries" % max_tries)
        else:
            j += 1 # increments only if no exception

    return(X)
#########################################################################
# Functions

def likelihood_evaluater_pool(theta_array,alpha_array,this_model,inference,path_dir, num_cores=mp.cpu_count() ):

    x = theta_array
    y = alpha_array
    X, Y = np.meshgrid(x, y)
    Z=np.zeros((len(x), len(y)))
    grid_list=[]
    for i in range(0,len(y)):
        for j in range(0,len(x)):
            grid_list.append(np.array((x[j],y[i])))

    if inference=="beta_likelihood":
        print('Evaluating the Beta likelihood using ' + str(num_cores) + ' cores ' )

        num_tasks = len(grid_list)
        with mp.Pool(num_cores) as pool:
            result_list= list( tqdm(pool.imap(this_model.beta_likelihood, grid_list ), total=num_tasks))
        pool.close()
        pool.join()
        result_list=np.array(result_list)
        Z=result_list.reshape((len(x), len(y))) #maybe opposite orientation


    if inference=='lamperti_likelihood_SDE_approx':
        print('Evaluating approx. lamperti using ' + str(num_cores) + ' cores ' )

        num_tasks = len(grid_list)
        with mp.Pool(num_cores) as pool:
            result_list= list( tqdm(pool.imap(this_model.lamperti_likelihood_SDE_approx, grid_list ), total=num_tasks))
        pool.close()
        pool.join()
        result_list=np.array(result_list)
        Z=result_list.reshape((len(x), len(y))) #maybe opposite orientation


        ########
    if inference=='lamperti_likelihood_linearized':
        print('Evaluating  linearized lamperti using ' + str(num_cores) + ' cores ' )

        num_tasks = len(grid_list)
        with mp.Pool(num_cores) as pool:
            result_list= list( tqdm(pool.imap(this_model.lamperti_likelihood_linearized, grid_list ), total=num_tasks))
        pool.close()
        pool.join()
        result_list=np.array(result_list)
        Z=result_list.reshape((len(x), len(y))) #maybe opposite orientation
        ########

    Zx,Zy=np.gradient(Z)

    Z_grad_x=np.zeros((len(y), len(x)))
    Z_grad_y=np.zeros((len(y), len(x)))
    for i in range(0,len(x)):
        for j in range(0,len(y)):
            Z_grad_x[i,j] = Zx[i][j]
            Z_grad_y[i,j] = Zy[i][j]

    Zxx,Zxy =np.gradient(Z_grad_x)
    Zyx,Zyy =np.gradient(Z_grad_y)

    Z_grad_xx=np.zeros((len(y), len(x)))
    Z_grad_yy=np.zeros((len(y), len(x)))
    Z_grad_xy=np.zeros((len(y), len(x)))
    Z_grad_yx=np.zeros((len(y), len(x)))

    for i in range(0,len(x)):
        for j in range(0,len(y)):
            Z_grad_xx[i,j] = Zxx[i][j]
            Z_grad_yy[i,j] = Zyy[i][j]
            Z_grad_xy[i,j] = Zxy[i][j]
            Z_grad_yx[i,j] = Zyx[i][j]

    np.save(path_dir+'/theta', X)
    np.save(path_dir+'/alpha', Y)
    np.save(path_dir+'/value', Z)
    np.save(path_dir+'/grad_x', Z_grad_x)
    np.save(path_dir+'/grad_y', Z_grad_y)
    np.save(path_dir+'/grad_xx', Z_grad_xx)
    np.save(path_dir+'/grad_yy', Z_grad_yy)
    np.save(path_dir+'/grad_xy', Z_grad_xy)
    np.save(path_dir+'/grad_yx', Z_grad_yx)
    print('output arrays save successfully in ' + path_dir )
    return()

def likelihood_evaluater(theta_array,alpha_array,this_model,inference,path_dir ):

    x = theta_array
    y = alpha_array
    X, Y = np.meshgrid(x, y)
    Z=np.zeros((len(x), len(y)))

    if inference=="beta_likelihood":
        print('Evaluating using the Beta likelihood')
        jobs = []
        pipe_list = []
        for i in range(0,len(y)):
            for j in range(0,len(x)):
                recv_end, send_end_to_evaluater = mp.Pipe(False)
                p = mp.Process(target=this_model.beta_likelihood, args=(np.array((x[j],y[i])), send_end_to_evaluater))
                jobs.append(p)
                pipe_list.append(recv_end)
                p.start()

        for proc in jobs:
            proc.join()
        result_list = [x.recv() for x in pipe_list]
        result_list=np.array(result_list)
        print(result_list)
        Z=result_list.reshape((len(x), len(y))) #maybe opposite orientation

    if inference=='lamperti_likelihood_SDE_approx':
        print('Evaluating using approx. lamperti ')
        jobs = []
        pipe_list = []
        for i in range(0,len(y)):
            for j in range(0,len(x)):
                recv_end, send_end_to_evaluater = mp.Pipe(False)
                p = mp.Process(target=this_model.lamperti_likelihood_SDE_approx, args=(np.array((x[j],y[i])), send_end_to_evaluater))
                jobs.append(p)
                pipe_list.append(recv_end)
                p.start()

        for proc in jobs:
            proc.join()
        result_list = [x.recv() for x in pipe_list]
        result_list=np.array(result_list)
        #print(result_list)
        Z=result_list.reshape((len(x), len(y))) #maybe opposite orientation
        ########
    if inference=='lamperti_likelihood_linearized':
        print('Evaluating using linearized lamperti')
        jobs = []
        pipe_list = []
        for i in range(0,len(y)):
            for j in range(0,len(x)):
                recv_end, send_end_to_evaluater = mp.Pipe(False)
                p = mp.Process(target=this_model.lamperti_likelihood_linearized, args=(np.array((x[j],y[i])), send_end_to_evaluater))
                jobs.append(p)
                pipe_list.append(recv_end)
                p.start()

        for proc in jobs:
            proc.join()
        result_list = [x.recv() for x in pipe_list]
        result_list=np.array(result_list)
        print(result_list)
        Z=result_list.reshape((len(x), len(y))) #maybe opposite orientation
        ########

    Zx,Zy=np.gradient(Z)

    Z_grad_x=np.zeros((len(y), len(x)))
    Z_grad_y=np.zeros((len(y), len(x)))
    for i in range(0,len(x)):
        for j in range(0,len(y)):
            Z_grad_x[i,j] = Zx[i][j]
            Z_grad_y[i,j] = Zy[i][j]

    Zxx,Zxy =np.gradient(Z_grad_x)
    Zyx,Zyy =np.gradient(Z_grad_y)

    Z_grad_xx=np.zeros((len(y), len(x)))
    Z_grad_yy=np.zeros((len(y), len(x)))
    Z_grad_xy=np.zeros((len(y), len(x)))
    Z_grad_yx=np.zeros((len(y), len(x)))

    for i in range(0,len(x)):
        for j in range(0,len(y)):
            Z_grad_xx[i,j] = Zxx[i][j]
            Z_grad_yy[i,j] = Zyy[i][j]
            Z_grad_xy[i,j] = Zxy[i][j]
            Z_grad_yx[i,j] = Zyx[i][j]

    np.save(path_dir+'/theta', X)
    np.save(path_dir+'/alpha', Y)
    np.save(path_dir+'/value', Z)
    np.save(path_dir+'/grad_x', Z_grad_x)
    np.save(path_dir+'/grad_y', Z_grad_y)
    np.save(path_dir+'/grad_xx', Z_grad_xx)
    np.save(path_dir+'/grad_yy', Z_grad_yy)
    np.save(path_dir+'/grad_xy', Z_grad_xy)
    np.save(path_dir+'/grad_yx', Z_grad_yx)
    print('output arrays save successfully !')
    return()


def data_check_interpolate(forecast_with_data, interpolation_points=72):
    #input array from np.load('forecast_with_data.npy')
    forecast_data_inter= np.zeros((3,forecast_with_data.shape[0], interpolation_points ))
    x = np.arange(1,72+1,1)
    xnew = np.linspace(0,72,interpolation_points) #adds one starting point as intialization
    counter=0
    i=0;j=0
    while j< forecast_with_data.shape[0]: #fix redundant indicies
        y1=forecast_with_data[j,:,1]
        y2=forecast_with_data[j,:,2] #linear interplation
        inter_1=interpolate.interp1d(x, y1, fill_value='extrapolate')
        inter_2=interpolate.interp1d(x, y2, fill_value='extrapolate')
        p=inter_1(xnew)
        d=inter_2(xnew)

        if  sum(1*np.isnan(p))==0 and \
            sum(1*np.isnan(d))==0 and \
            np.max(p)<=1 and \
            np.min(p)>=0 and \
            np.max(d)<=1 and \
            np.min(d)>=0:

            forecast_data_inter[0,i,:]= p
            forecast_data_inter[1,i,:]= d
            forecast_data_inter[2,i,0]= forecast_with_data[j,0,0]
            i+=1
        j+=1
    n_paths=i
    print('paths left', i)
    forecast_data_inter=forecast_data_inter[:,:n_paths,:]
    return(forecast_data_inter)

def theta_adjust(theta, forecast):
    p=forecast
    N=len(p)
    p_dot=np.diff(p)*N

    zero_drift = np.diff(p)*N/theta + p[:-1]
    zero_drift_fixed=np.zeros_like(zero_drift)
    theta_adjusted=np.zeros_like(zero_drift)

    for i in range(0,len(zero_drift)):
        zero_drift_fixed[i] = p_dot[i]/max(theta, abs(p_dot[i])/(min(p[i], 1-p[i]) + 1e-16) ) + p[i]
        theta_adjusted[i]= max(theta, abs(p_dot[i])/(min(p[i], 1-p[i]) + 1e-16 ) )

    return(theta_adjusted, zero_drift_fixed )

def theta_t(theta, p, pdot): # (theta_0, batches). From batches, we only want $\dot{p}$.
    return(max(theta, abs(pdot)/(min(p, 1-p) + 1e-16)))

def empirical_Confidence_Interval_plots(forecast_data_inter,hours,\
    real_in, disct_in,list_forecast_number,\
    dir_path):
    if isinstance(hours, int) == False:
        raise ValueError('hours is not an integer. Please use integer values only.')
    os.makedirs(dir_path + '/' + str(hours) + 'hr',exist_ok = True)
    os.makedirs(dir_path + '/' + 'data_plots', exist_ok = True)
    N    = disct_in.N
    freq = (N+1)*1 # Evert 10 minutes.
    #interpolation_points=disct_in.N
    num_forecasts = len(list_forecast_number) # To follow the custom forecast order.
    q975          = np.zeros((num_forecasts,freq))
    q025          = np.zeros((num_forecasts,freq))
    q95           = np.zeros((num_forecasts,freq))
    q05           = np.zeros((num_forecasts,freq))
    q50           = np.zeros((num_forecasts,freq))
    q75           = np.zeros((num_forecasts,freq))
    q25           = np.zeros((num_forecasts,freq))

    xnew = np.linspace(0,N,freq) #np.linspace(0,N,N+1)
    x = np.linspace(1,N,N)

    theta= real_in.mu #7.8155
    alpha= real_in.sigma #1.072
    j=0
    for k in list_forecast_number: #   #range(0,n_paths): #n_paths
        p=forecast_data_inter[2,k,:N] #obtain cleansed forecast
        inter_1=interpolate.interp1d(x, p, fill_value='extrapolate')
        p=inter_1(xnew)
        d=forecast_data_inter[1,k,:N]
        dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[0,k,0])
        dt=1
        M_test=disct_in.M

        fig=plt.figure(2,figsize=(10, 4))
        fig.clf()
        plt.plot(xnew,p, label='forecast')
        plt.plot(x,d, label='actual production')

        plt.xlim(1, 427)
        plt.ylim(-0.1, 1.1)
        plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

        plt.xticks( fontsize = 20);
        plt.yticks( fontsize = 20);
        plt.xlabel('Time [min]',fontsize = 24)
        plt.ylabel('Power',fontsize = 24)
        plt.legend( prop={'size': 15})
        plt.savefig(dir_path +'/data_plots' + '/' + str(k)+'.pdf', bbox_inches="tight")


        theta_adjusted, zero_drift_fixed=theta_adjust(theta,p)

        real_1 = real(theta_adjusted, alpha ) #theta_array

        disct_temp = disct(freq,dt,M_test)

        X = np.empty((M_test,N+1));
        X= gen_X_normal_euler_DT_modified(X0=p[0],disct=disct_temp,real=real_1,forecast=p)
        #plt.plot(X[0,:])

        for i in range(0,freq):
            q975[k,i]=np.quantile(X[:,i], 0.975)
            q025[k,i]=np.quantile(X[:,i], 0.025)

            q95[k,i]=np.quantile(X[:,i], 0.95)
            q05[k,i]=np.quantile(X[:,i], 0.05)

            q50[k,i]=np.quantile(X[:,i], 0.5)

            q75[k,i]=np.quantile(X[:,i], 0.75)
            q25[k,i]=np.quantile(X[:,i], 0.25)

        #plotting
        fig=plt.figure(2,figsize=(10, 4))
        fig.clf()

        plt.xlim(1, hours)
        plt.ylim(-0.1, 1.1)

        plt.fill_between(xnew, q75[j,:],q25[j,:],color='k',alpha=0.2, label='90% CI', edgecolor=None)
        plt.fill_between(xnew, q95[j,:],q05[j,:],color='k',alpha=0.3, label='50% CI', edgecolor=None)

        #plt.plot(xnew,np.mean(normal_X_derivative_tracking_Euler, axis=0),'c-', label='Mean')

        #plt.plot(xnew,q50,'y-', label='Median')
        #dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[2,k,0])
        plt.plot(xnew,p, 'r-', label='forecast',linewidth=3)
        #plt.plot(xnew[:-1],zero_drift_fixed, 'y-', label='Zero Drift Line',linewidth=1)
        plt.plot(x,d , 'y-', label='actual production',linewidth=3)
        plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

        plt.xticks( fontsize = 20);
        plt.yticks( fontsize = 20);
        plt.xlabel('Time [hr]',fontsize = 24)
        plt.ylabel('Power',fontsize = 24)
        plt.legend( prop={'size': 12})
        plt.savefig(dir_path+'/'+str(hours)+'hr/'+str(k)+'.pdf', bbox_inches="tight")

        j+=1;

    #save output

    file_object  = open(dir_path+'/parameter.info', 'w')
    note='Confidence intervals for ' + str(num_forecasts)+' forecasts of '\
    + str(N) +' hours.'+ '\n' +'Parameters used are theta= '\
    +str(theta) + ' and alpha= ' + str(alpha) + '.Generated by '+str(M_test) \
    + ' simulations.'
    file_object.write(note)
    file_object.close()

    np.save(dir_path + '/forecast_list', list_forecast_number)
    np.save(dir_path + '/N', N)
    np.save(dir_path + '/theta_0', theta)
    np.save(dir_path + '/alpha', alpha)
    np.save(dir_path + '/q975', q975)
    np.save(dir_path + '/q025', q025)
    np.save(dir_path + '/q95', q95)
    np.save(dir_path + '/q05', q05)
    np.save(dir_path + '/q50', q50)
    np.save(dir_path + '/q75', q75)
    np.save(dir_path + '/q25', q25)

def empirical_Confidence_Interval_plots_old(forecast_data_inter,\
    real_in, disct_in,list_forecast_number,\
    dir_path):
    os.mkdir(dir_path + '/72hr')
    os.mkdir(dir_path + '/6hr')
    N = disct_in.N
    #interpolation_points=disct_in.N
    num_forecasts = len(list_forecast_number) # to follow the custom forecast order
    q975          = np.zeros((num_forecasts,N+1))
    q025          = np.zeros((num_forecasts,N+1))
    q95           = np.zeros((num_forecasts,N+1))
    q05           = np.zeros((num_forecasts,N+1))
    q50           = np.zeros((num_forecasts,N+1))
    q75           = np.zeros((num_forecasts,N+1))
    q25           = np.zeros((num_forecasts,N+1))

    xnew = np.linspace(0,N,N+1)
    x = np.linspace(1,N,N)

    theta= real_in.mu #7.8155
    alpha= real_in.sigma #1.072
    j=0
    for k in list_forecast_number: #   #range(0,n_paths): #n_paths
        p=forecast_data_inter[0,k,:N] #obtain cleansed forecast
        inter_1=interpolate.interp1d(x, p, fill_value='extrapolate')
        p=inter_1(xnew)
        d=forecast_data_inter[1,k,:N]
        dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[0,k,0])
        dt=1
        M_test=disct_in.M

        fig=plt.figure(2,figsize=(10, 4))
        fig.clf()
        plt.plot(xnew,p, label='forecast')
        plt.plot(xnew,d, label='actual production')

        plt.xlim(1, 427)
        plt.ylim(-0.1, 1.1)
        plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

        plt.xticks( fontsize = 20);
        plt.yticks( fontsize = 20);
        plt.xlabel('Time [hr]',fontsize = 24)
        plt.ylabel('Power',fontsize = 24)
        plt.legend( prop={'size': 15})
        plt.savefig('Forecast_data_'+ str(k)+'.pdf', bbox_inches="tight")


        theta_adjusted, zero_drift_fixed=theta_adjust(theta,p)

        real_1 = real(theta_adjusted, alpha ) #theta_array

        disct_temp = disct(N+1,dt,M_test)

        X = np.empty((M_test,N+1));
        X= gen_X_normal_euler_DT_modified(X0=p[0],disct=disct_temp,real=real_1,forecast=p)
        #plt.plot(X[0,:])

        for i in range(0,N+1):
            q975[k,i]=np.quantile(X[:,i], 0.975)
            q025[k,i]=np.quantile(X[:,i], 0.025)

            q95[k,i]=np.quantile(X[:,i], 0.95)
            q05[k,i]=np.quantile(X[:,i], 0.05)

            q50[k,i]=np.quantile(X[:,i], 0.5)

            q75[k,i]=np.quantile(X[:,i], 0.75)
            q25[k,i]=np.quantile(X[:,i], 0.25)

        #plotting
        fig=plt.figure(2,figsize=(10, 4))
        fig.clf()

        plt.xlim(1, N)
        plt.ylim(-0.1, 1.1)

        plt.fill_between(xnew, q75[j,:],q25[j,:],color='k',alpha=0.2, label='90% CI', edgecolor=None)
        plt.fill_between(xnew, q95[j,:],q05[j,:],color='k',alpha=0.3, label='50% CI', edgecolor=None)

        #plt.plot(xnew,np.mean(normal_X_derivative_tracking_Euler, axis=0),'c-', label='Mean')

        #plt.plot(xnew,q50,'y-', label='Median')
        #dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[2,k,0])
        plt.plot(xnew,p, 'r-', label='forecast',linewidth=3)
        #plt.plot(xnew[:-1],zero_drift_fixed, 'y-', label='Zero Drift Line',linewidth=1)
        plt.plot(x,d , 'y-', label='actual production',linewidth=3)
        plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

        plt.xticks( fontsize = 20);
        plt.yticks( fontsize = 20);
        plt.xlabel('Time [hr]',fontsize = 24)
        plt.ylabel('Power',fontsize = 24)
        plt.legend( prop={'size': 15})
        plt.savefig(dir_path+'/72hr/'+'72hr_'+str(k)+'.pdf', bbox_inches="tight")

        fig= plt.figure(2,figsize=(10, 4))
        fig.clf()

        plt.xlim(1, 7)
        plt.ylim(-0.1, 1.1)

        plt.fill_between(xnew, q75[j,:],q25[j,:],color='k',alpha=0.2, label='90% CI', edgecolor=None)
        plt.fill_between(xnew, q95[j,:],q05[j,:],color='k',alpha=0.3, label='50% CI', edgecolor=None)

        #plt.plot(xnew,np.mean(normal_X_derivative_tracking_Euler, axis=0),'c-', label='Mean')

        #plt.plot(xnew,q50,'y-', label='Median')

        plt.plot(xnew,p, 'r-', label='forecast',linewidth=3)
        plt.plot(x,d , 'y-', label='actual production',linewidth=3)

        #plt.plot(xnew[:-1],zero_drift_fixed, 'y-', label='Zero Drift Line',linewidth=1)

        plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24) #,fontsize=10 Forecast Confidence Intervals
        plt.xlabel('Time [hr]',fontsize = 24)
        plt.ylabel('Power',fontsize = 24)
        plt.legend( prop={'size': 15})
        plt.xlim(1, 7)
        plt.ylim(-0.1, 1.1)
        plt.xticks( fontsize = 20);
        plt.yticks( fontsize = 20);
        plt.savefig( dir_path+'/6hr/'+'6hr_' +str(k)+'.pdf', bbox_inches="tight")

        j+=1;

    #save output

    file_object  = open(dir_path+'/parameter.info', 'w')
    note='Confidence intervals for ' + str(num_forecasts)+' forecasts of '\
    + str(N) +' hours.'+ '\n' +'Parameters used are theta= '\
    +str(theta) + ' and alpha= ' + str(alpha) + '.Generated by '+str(M_test) \
    + ' simulations.'
    file_object.write(note)
    file_object.close()

    np.save(dir_path + '/forecast_list', list_forecast_number)
    np.save(dir_path + '/N', N)
    np.save(dir_path + '/theta_0', theta)
    np.save(dir_path + '/alpha', alpha)
    np.save(dir_path + '/q975', q975)
    np.save(dir_path + '/q025', q025)
    np.save(dir_path + '/q95', q95)
    np.save(dir_path + '/q05', q05)
    np.save(dir_path + '/q50', q50)
    np.save(dir_path + '/q75', q75)
    np.save(dir_path + '/q25', q25)



def path_simulator(forecast_data_inter,hours,\
    real_in, disct_in,list_forecast_number,\
    dir_path):

    if isinstance(hours, int)==False:
        raise ValueError('hours is not an integer. Please use integer values only.')

    os.makedirs(dir_path + '/'+str(hours)+'hr',exist_ok=True)
    #os.makedirs(dir_path + '/6hr',exist_ok=True)
    N=disct_in.N
    #interpolation_points=disct_in.N
    num_forecasts=len(list_forecast_number) # to follow the custom forecast order

    xnew = np.linspace(0,N,(N+1)) #np.linspace(0,N,N+1)
    x = np.linspace(1,N,N)


    theta= real_in.mu #7.8155
    alpha= real_in.sigma #1.072
    j=0
    for k in list_forecast_number: #   #range(0,n_paths): #n_paths
        p=forecast_data_inter[2,k,:N] #obtain cleansed forecast
        inter_1=interpolate.interp1d(x, p, fill_value='extrapolate')
        p=inter_1(xnew)
        d=forecast_data_inter[1,k,:N]
        dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[0,k,0])
        dt=1
        M_test=disct_in.M

        # THIS: plots the forcast data.
        fig=plt.figure(2,figsize=(10, 4))
        fig.clf()
        # plt.plot(xnew,p, label='forecast')
        plt.plot(x,d, label='actual production')

        # plt.xlim(1, 73)
        # plt.ylim(-0.1, 1.1)
        # plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

        # plt.xticks( fontsize = 20);
        # plt.yticks( fontsize = 20);
        # plt.xlabel('Time [hr]',fontsize = 24)
        # plt.ylabel('Power',fontsize = 24)
        # plt.legend( prop={'size': 15})
        # plt.savefig('Forecast_data_'+ str(k)+'.pdf', bbox_inches="tight")

        theta_adjusted, zero_drift_fixed=theta_adjust(theta,p)

        real_1 = real(theta_adjusted, alpha ) #theta_array

        disct_temp = disct((N+1),dt,M_test)

        X = np.empty((M_test,(N+1)));
        X= gen_X_normal_euler_DT_modified(X0=p[0],disct=disct_temp,real=real_1,forecast=p)

        #plotting
        # fig=plt.figure(2,figsize=(10, 4))
        # fig.clf()

        plt.xlim(1, hours)
        plt.ylim(-0.1, 1.1)

        plt.plot(xnew,p, 'k-', label='forecast',linewidth=5)
        #plt.plot(x,d , 'y-', label='actual production',linewidth=2)

        for i in range(0,len(X)):
            plt.plot(xnew,X[i] , '-',linewidth=1,markersize=4)
        plt.plot(xnew,X[i] , '-', label='simulated production',linewidth=1,markersize=4)

        plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

        plt.xticks( fontsize = 20);
        plt.yticks( fontsize = 20);
        plt.xlabel('Time [min]',fontsize = 24)
        plt.ylabel('Power',fontsize = 24)
        plt.legend( prop={'size': 20})
        plt.savefig(dir_path+'/'+ str(hours) +'hr/'+str(k)+'.pdf', bbox_inches="tight");

    file_object  = open(dir_path+'/parameter.info', 'w')
    note='Simulation of ' + str(num_forecasts)+' forecasts spanning '\
    + str(N) +' hours.'+ '\n' +'Parameters used are theta= '\
    +str(theta) + ' and alpha= ' + str(alpha) + '.Generated by '+str(M_test) \
    + ' simulations.'
    file_object.write(note)
    file_object.close()
    np.save(dir_path + '/forecast_list', list_forecast_number)



def compute_path_moments(j,N,X,p,theta_adjusted, alpha, send_end):
    m_1=0  #fix initial condition as zero error in the first point
    m_2=0
    dN=1/N
    L_n=0
    for i in range(0,N-1): #start from 1 or zero #Implement Tripizoidal scheme
        m_1=X[j,i]*np.exp(- dN*theta_adjusted[i])
        m_2=  (X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta_adjusted[i]*p[j,i]*(1-p[j,i])*(1-2*p[j,i] )) + \
                     alpha*theta_adjusted[i]*p[j,i]**2*(1-p[j,i])**2) ) /(  1+ dN*2*(theta_adjusted[i]+alpha*theta_adjusted[i]*p[j,i]*(1-p[j,i]) ))
        a=m_1;
        b=m_2- m_1**2;

        if np.isnan(a): print('a is nan')
        if np.isnan(b): print('b is nan')
        if not np.isfinite(a): print('a is infinite')
        if not np.isfinite(b): print('b is infinite')

        if b==0: b=1e-16
        if b<0:
            b=1e-16

        beta_param_alpha = - ( (1+a)*(a**2 + b - 1) )/ (2*b)
        beta_param_beta  = ( (a-1)*(a**2 + b - 1) ) / (2*b)

        L_n_current= (beta_param_alpha-1 )*np.log(  (X[j,i+1]+1)/2 ) +\
         (beta_param_beta-1)*np.log(1-( X[j,i+1] +1)/2 )-\
         scipy.special.betaln(beta_param_alpha, beta_param_beta)
        L_n = L_n + L_n_current
    send_end.send(L_n)

# def moment_compute(dN, X_prev, X_next, p_prev, theta_current,       alpha_current,send_end ):
#         m_1=X_prev*np.exp(- dN*theta_current)
#         m_2=  (X_prev**2 + 2*dN*( X_prev*(alpha_current*theta_current*p_prev*(1-p_prev)*(1-2*p_prev )) + \
#                      alpha_current*theta_current*p_prev**2*(1-p_prev)**2) ) /(  1+ dN*2*(theta_current+alpha_current*theta_current*p_prev*(1-p_prev) ))
#         a=m_1;
#         b=m_2- m_1**2;
#
#         if b==0: b=1e-16
#         if b<0:
#             b=1e-16
#
#         beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
#         beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)
#
#         L_n_current=(beta_param_alpha-1 )*np.log(  (X_next+1)/2 ) +\
#          (beta_param_beta-1)*np.log(1-( X_next +1)/2 )-\
#          scipy.special.betaln(beta_param_alpha, beta_param_beta)
#
#         send_end.send(L_n_current)

def moment_compute(dN, X_prev, X_next, p_prev,p_next, theta_prev, theta_next,alpha ):
    start=0 ; end=1
    N_INT=100
    x=np.linspace(start,end, N_INT )
    X_INT=np.interp(x,[start,end], [X_prev,X_next] )
    P_INT=np.interp(x,[start,end], [p_prev,p_next] )
    theta_INT=np.interp(x,[start,end], [theta_prev,theta_next] )
    d_INT= 10/N_INT #change time scale to choose which units ? we choose minutes
    m_1=X_prev
    m_2=X_prev**2
    m_1_prev=0
    I_1_prev=0
    I_1=0
    I_2=0
    for q in range( 0, len(X_INT)-1 ):

        m_1_prev=m_1
        I_1_prev=I_1

        I_1= I_1 + (theta_INT[q] + theta_INT[q+1])*d_INT/2

        m_1= X_prev*np.exp(-I_1)

        I_2 = I_2 + ( ( 2*theta_INT[q+1]*m_1*(1-2*P_INT[q+1])    + 2*theta_INT[q+1]*P_INT[q+1]*(1-P_INT[q+1]) )*np.exp(-2*(1+alpha )*I_1)  + ( 2*theta_INT[q]*m_1_prev*(1-2*P_INT[q])    + 2*theta_INT[q]*P_INT[q]*(1-P_INT[q]) )*np.exp(-2*(1+alpha )*I_1_prev)   )*d_INT/2

        # m_1_prev=m_1
        # m_1 = m_1*np.exp(-1*(theta_INT[q] + theta_INT[q+1])*d_INT/2 )
        #
        # m_2=(m_2 + 0.5*d_INT*( -2*m_2*theta_INT[q]*(1+alpha) + 2*m_1_prev*alpha*theta_INT[q]*(1-2*P_INT[q]) + 2*alpha*theta_INT[q]*P_INT[q]*(1-P_INT[q]) + 2*m_1*alpha*theta_INT[q+1]*(1-2*P_INT[q+1]) + 2*alpha*theta_INT[q+1]*P_INT[q+1]*(1-P_INT[q+1]) ) )/(1-d_INT*theta_INT[q+1]*(1+alpha))

    m_1= X_prev*np.exp(-I_1)
    m_2=X_prev**2*np.exp(-2*(1+alpha)*I_1) + alpha*I_2

    return(m_1,m_2)

def moment_compute_approx_lamperti(dN, X_prev, X_next, p_prev,p_next, theta_prev, theta_next,alpha ):

    p_func = interpolate.interp1d([0,1], [p_prev,p_next] , kind='linear')
    theta_func = interpolate.interp1d([0,1], [theta_prev,theta_next], kind='linear')

    fun=lambda t, m: approx_lamperti_ODE_RHS(t, m, p_func=p_func,theta_func=theta_func, alpha=alpha)
    sol = solve_ivp(fun, [0, 1], [X_prev,X_prev**2-X_prev] ) #, rtol=1e-2, atol=1e-2
    m_1= sol.y[0,-1]
    var= sol.y[1,-1]
    return(m_1,var)

    # start=0 ; end=1
    # N_INT=10
    # x=np.linspace(start,end, N_INT )
    # X_INT=np.interp(x,[start,end], [X_prev,X_next] )
    # P_INT=np.interp(x,[start,end], [p_prev,p_next] )
    # theta_INT=np.interp(x,[start,end], [theta_prev,theta_next] )
    # d_INT= (10/N_INT)
    # I_1=0
    # I_2=0
    # I_3=0
    # I_4=0
    # I_1_prev=0
    # I_3_prev=0
    # I_3_list=[]
    # K=0
    # mean_current=0
    # for q in range( 0, len(X_INT)-1 ):
    #     I_1_prev=I_1
    #     mean_prev=mean_current
    #     I_3_prev=I_3
    #
    #     I_1= I_1 + ( theta_INT[q] + theta_INT[q+1] )*d_INT/2
    #
    #     I_2= I_2 + (theta_INT[q+1]*(2*P_INT[q+1]-1)*np.exp((1-alpha)*I_1) + theta_INT[q]*(2*P_INT[q]-1)*np.exp((1-alpha)*I_1_prev  ) )*d_INT/2
    #
    #     K=np.exp(-1*(1-alpha)*I_1 )*I_2 + np.sin(X_prev)
    #
    #     if K>1:
    #         #print(K-1)
    #         K=1
    #     if K<-1:
    #         #print(K+1)
    #         K=-1
    #
    #     mean_current= np.arcsin(K)
    #
    #
    #     I_3= I_3 + ( theta_INT[q+1]*(2*P_INT[q+1]-1)*np.tan(mean_current)/np.cos(mean_current) + theta_INT[q+1]*(alpha - 1)/np.cos(mean_current)**2 + theta_INT[q]*(2*P_INT[q]-1)*np.tan(mean_prev)/np.cos(mean_prev) + theta_INT[q]*(alpha - 1)/np.cos(mean_prev)**2 )*d_INT/2
    #
    #     I_3_list.append(I_3)

        #print(I_3)


        # if not np.isfinite(theta_INT[q+1]): print('theta_INT[q+1] is infinite', theta_INT[q+1])
        # if not np.isfinite(theta_INT[q]): print('theta_INT[q] is infinite', theta_INT[q])
        # if not np.isfinite(I_3): print('I_3 is infinite', I_3)
        # if not np.isfinite(I_3_prev): print('I_3_prev is infinite', I_3_prev)


        #I_4= I_4 + (theta_INT[q+1]*np.exp(-2*I_3) + theta_INT[q]*np.exp(-2*I_3_prev) )*d_INT/2


        # if not np.isfinite(I_4):
        #     print('I_4 is infinite')
        #     print('theta_INT[q+1]= ' , theta_INT[q+1])
        #     print('mean_current= ' , mean_current )
        #     print('theta_INT[q+1]= ' , theta_INT[q+1])
        #     print('P_INT[q]= ' , P_INT[q])
        #
        # if np.isnan(I_3):
        #     print('I_3 is nan')
        # if np.isnan(I_4):
        #     print('I_4 is nan')

    # for t in range(0,len(X_INT)-1):
    #     I_4= I_4 + (theta_INT[t+1]*np.exp(-2*I_3_list[-t-1]) + theta_INT[t]*np.exp(-2*I_3_list[-t]) )*d_INT/2
    #
    # K=np.exp(-1*(1-alpha)*I_1 )*I_2 + np.sin(X_prev)
    # if K>1: K=1
    # if K<-1: K=-1
    # mean=np.arcsin(K)
    # var=alpha*I_4
    # if np.isnan(var):
    #     print('var is nan, we have I_3 and I_4 as ', I_3, I_4)

    # return(mean,var)





# def moment_compute_linear_lamperti(dN, X_prev, X_next, p_prev,p_next, theta_prev, theta_next,alpha ):
#     start=0 ; end=1
#     N_INT=100
#     x=np.linspace(start,end, N_INT )
#     X_INT=np.interp(x,[start,end], [X_prev,X_next] )
#     P_INT=np.interp(x,[start,end], [p_prev,p_next] )
#     theta_INT=np.interp(x,[start,end], [theta_prev,theta_next] )
#     d_INT= (10/N_INT)
#     I_3=0
#     I_1=0
#     I_2=0
#     I_4=0
#     I_3_prev=0
#     for q in range( 0, len(X_INT)-1 ):
#         I_3_prev=I_3
#         I_1= I_1 +( theta_INT[q] + theta_INT[q+1]  )*d_INT/2
#         for w in range(q, len(X_INT)-1 ):
#             I_3 = I_3 + ( theta_INT[q] + theta_INT[q+1]  )*d_INT/2
#
#         I_2 = I_2 + ( theta_INT[q]*(1-2*P_INT[q])*np.exp( (alpha-1)*I_3_prev ) + theta_INT[q+1]*(1-2*P_INT[q+1])*np.exp( (alpha-1)*I_3 ) )*d_INT/2
#         I_4 = I_4 + ( np.sqrt(2*alpha*theta_INT[q+1] )*np.exp( (alpha-1)*I_3 )   + np.sqrt(2*alpha*theta_INT[q] )*np.exp( (alpha-1)*I_3_prev )  )*d_INT/2
#
#     mean= X_prev*np.exp( (alpha-1)*I_1 ) - I_2
#     var= I_4
#
#     return(mean,var)

####################################################################################################
#helper Functions
#####################################################################################################
def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--', label='slope='+str(slope))

def abline_log(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.loglog(x_vals, y_vals, '--', label='slope='+str(slope))

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

    #############################################################################################333

def beta_ODE_RHS(t, m ,p_func,theta_func, alpha):
    m_1_RHS= -m[0]*theta_func(t)
    m_2_RHS= -2*m[1]*theta_func(t)*(1+alpha) + 2*alpha*theta_func(t)*m[0]*(1-2*p_func(t)) + 2*alpha*theta_func(t)*p_func(t)*( 1-p_func(t) )
    return(m_1_RHS, m_2_RHS)


def approx_lamperti_ODE_RHS(t, m ,p_func,theta_func, alpha):
    m_1_RHS= -m[0]*theta_func(t)*(1-alpha) - theta_func(t)*(1-2*p_func(t))
    var_RHS= 2*(theta_func(t)*(2*p_func(t) -1)*np.sin(m[0])/np.cos(m[0])**2 + theta_func(t)*(alpha -1 )/np.cos(m[0])**2  )*m[1] + 2*theta_func(t)*alpha
    return(m_1_RHS, var_RHS)

def linear_lamperti_ODE_RHS(t, m ,p_func,theta_func, alpha):
    m_1_RHS= m[0]*theta_func(t)*(alpha -1) - theta_func(t)*(1-2*p_func(t))
    var_RHS= (alpha -1)*theta_func(t)*m[1] + np.sqrt(2*theta_func(t)*alpha)
    return(m_1_RHS, var_RHS)


# p_func = interpolate.interp1d([0,1], [0.4,0.3] , kind='linear')
#
# theta_func = interpolate.interp1d([0,1], [10,13] , kind='linear')
#
# beta_ODE_RHS(0.5, m=np.array((1,2)), p_func=p_func , theta_func=theta_func, alpha=0.1  )
#
# fun=lambda t, m: beta_ODE_RHS(t, m, p_func=p_func,theta_func=theta_func, alpha=0.1)
# sol = solve_ivp(fun, [0, 1], [2,8], method='RK45', t_eval=[1])
# m_1= np.asscalar(sol.y[0])
# m_2= np.asscalar(sol.y[1])

def beta_moment(dN, X_prev, X_next, p_prev, p_next, theta_prev, theta_next, alpha):

    # dN = 1/N.
    transition_duration = dN # Before 21/11/2019, we had transition_duration = 1.

    p_func     = interpolate.interp1d([0,transition_duration], [p_prev,p_next],         kind = 'linear', fill_value = 'extrapolate')
    theta_func = interpolate.interp1d([0,transition_duration], [theta_prev,theta_next], kind = 'linear', fill_value = 'extrapolate')

    fun = lambda t, m: beta_ODE_RHS(t, m, p_func = p_func, theta_func = theta_func, alpha = alpha)
    sol = solve_ivp(fun, [0, transition_duration], [X_prev, X_prev**2], rtol = 1e-2, atol = 1e-2)
    m_1 = sol.y[0,-1]
    m_2 = sol.y[1,-1]

    return(m_1,m_2)

def linear_lamperti_moment(dN, X_prev, X_next, p_prev, p_next, theta_prev, theta_next, alpha ):

    p_func = interpolate.interp1d([0,1], [p_prev,p_next], kind = 'linear', fill_value = 'extrapolate')

    theta_func = interpolate.interp1d([0,1], [theta_prev,theta_next], kind = 'linear',fill_value = 'extrapolate')

    fun = lambda t, m: linear_lamperti_ODE_RHS(t, m, p_func = p_func, theta_func = theta_func, alpha = alpha)
    sol = solve_ivp(fun, [0, 1], [X_prev,X_prev**2-X_prev], rtol = 1e-2, atol = 1e-2)
    m_1 = sol.y[0,-1]
    var = sol.y[1,-1]

    return(m_1,var)

def Ellipse (H_opt, mean):
    eig_values_opt, eig_vect_opt=LA.eig(H_opt)
    a_opt=1/np.sqrt(np.max(np.abs(eig_values_opt)))  #radius on the x-axis
    b_opt=1/np.sqrt( np.min(np.abs(eig_values_opt))) #radius on the y-axis
    width = a_opt**2
    height = b_opt**2
    angle =math.acos(sum(eig_vect_opt[1]*np.array((0,1))));angle
    return(width,height,angle, eig_vect_opt )

f= lambda x,y: 0
f= np.vectorize(f,otypes=[np.float])

def init(line, point):
    line.set_data([], [])
    line.set_3d_properties([])
    point.set_data([], [])
    point.set_3d_properties([])
    return(line, point)

def animate(line,point,path,i):
    line.set_data(path[0,:i], path[1,:i])
    line.set_3d_properties(f(*path[::,:i]))
    point.set_data(path[0,i-1:i], path[1,i-1:i])
    point.set_3d_properties(f(*path[::,i-1:i]))
    return(line, point)
