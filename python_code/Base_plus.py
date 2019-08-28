import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats
from scipy.optimize import minimize
from scipy.stats import multivariate_normal, norm
from mpl_toolkits.mplot3d import Axes3D
import operator
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import math
import random
from scipy import interpolate
from math import pi
from numpy import linalg as LA
import datetime as dtM
import os
##############
#constants
Nfeval=0
##############
class disct:
  def __init__(self, N,dt, M ):
    self.N = N
    self.dt = dt
    self.M = M

class real:
  def __init__(self, mu,sigma ):
    self.mu = mu
    self.sigma = sigma

class optimal(real):
  def __init__(self):
    real.__init__(self, 0,0 )

class model_modified_drift:
    def __init__(self, disct, data,forecast ):
        self.disct = disct
        self.data = data
        #self.ic = ic
        self.forecast = forecast
        self.optimal=optimal()
        self.Nfeval=0

    def lamperti_likelihood_SDE_approx(self,param):
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
        counter=0
        for j in range(0,M):
            h=0
            g=0
            mean=0
            var=0
            L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            for i in range(0,N-1): #start from 1 or zero #Implement Tripizoidal scheme
                #theta_adjusted[i]
                g= g + theta_adjusted[i]*(1-2*p[j,i])*h*dN #this h is approx.
                h= h + theta_adjusted[i]*(alpha-1)*dN

                mean= X[j,i]* np.exp(h) - g
                var= var + 2*alpha*theta*np.exp(2*h) *dN

                if np.isnan(mean): print('mean is nan')
                if np.isnan(var): print('var is nan')

                #norm(0, 1).pdf(10)

                L_n = L_n  - ( X[j,i+1]- X[j,i] - mean   )**2/(2*var) -0.5*np.log(2*math.pi*var)

        # display information
        #if self.Nfeval%5 == 0:
        print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        self.Nfeval += 1
        #print('counter is ', counter)

        return(-1*L_m)

    def beta_likelihood(self,param):
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
        counter=False
        for j in range(0,M):
            m_1=0  #fix initial condition as zero error in the first point
            m_2=0
            L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            counter=False
            h=0
            for i in range(0,N-1): #start from 1 or zero #Implement Tripizoidal scheme
                #element=(theta_adjusted[i] + theta_adjusted[i+1])/2 * dN
                #m_1 = X[j,i]*np.exp(-1*element )

                # m_2= X[j,i]**2 * np.exp(-1*(1+alpha)*element) + \
                #     alpha*dN*( ( theta_adjusted[i+1]*m_1*(1-p[j,i+1]- p[j,i+1]**2) +2 )*np.exp(-1*(1+alpha)*element ) +\
                #     theta_adjusted[i]*X[j,i]*(1-p[j,i] - p[j,i]**2 ) +2 )/2

                m_1=X[j,i]*np.exp(- dN*theta_adjusted[i])

                m_2=  (X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta_adjusted[i]*p[j,i]*(1-p[j,i])*(1-2*p[j,i] )) + \
                            alpha*theta_adjusted[i]*p[j,i]**2*(1-p[j,i])**2) ) /(  1+ dN*2*(theta_adjusted[i]+alpha*theta_adjusted[i]*p[j,i]*(1-p[j,i]) ))

                #m_222=X[j,i]**2 * np.exp( -(1+alpha)*theta_adjusted[i]*dN ) + dN*(alpha*theta_adjusted[i]*X[j,i]*(1-p[j,i]-p[j,i]**2) + 2)*np.exp(-(1+alpha)*theta_adjusted[i]*dN )

                a=m_1;
                b=m_2- m_1**2;
                #print(m_2,m_22)
                #if b<0: print('negative variance !! help')

                if np.isnan(a): print('a is nan')
                if np.isnan(b): print('b is nan')
                if not np.isfinite(a): print('a is infinite')
                if not np.isfinite(b): print('b is infinite')

                if b==0: b=1e-16
                if b<0:
                    b=1e-16
                    counter = True

                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                L_n = L_n +  (beta_param_alpha-1 )*np.log(  (X[j,i+1]+1)/2 ) +\
                 (beta_param_beta-1)*np.log(1-( X[j,i+1] +1)/2 )-\
                 scipy.special.betaln(beta_param_alpha, beta_param_beta)
        # display information
        #if self.Nfeval%5 == 0:
        #print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        self.Nfeval += 1
        #print(counter)


        return(-1*L_m)

    def lamperti_likelihood_moment_approx(self,param):
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
        counter=0
        for j in range(0,M):
            h=0
            g=0
            mean=0
            var=0
            L_m = L_m + L_n
            theta_adjusted, _ = theta_adjust( theta , p[j,:] )
            for i in range(0,N-1): #start from 1 or zero #Implement Tripizoidal scheme
                #theta_adjusted[i]
                g= g + theta_adjusted[i]*(2*p[j,i]-1)*h*dN #this h is approx.
                h= h + theta_adjusted[i]*(1-2*alpha)*dN
                arc_in=np.exp(-h)*(g+np.sin(X[j,i]) )

                if (arc_in >1):
                    arc_in=1
                    counter +=1
                if (arc_in <-1):
                    arc_in=-1
                    counter +=1
                if np.isnan(arc_in): print('arc_in is nan')

                mean= np.arcsin(arc_in )

                var= var + 2*alpha*theta*dN
                if np.isnan(mean): print('mean is nan')
                if np.isnan(var): print('var is nan')

                #norm(0, 1).pdf(10)

                L_n = L_n  - ( X[j,i+1]- X[j,i] - mean   )**2/(2*var) -0.5*np.log(2*math.pi*var)

        # display information
        #if self.Nfeval%5 == 0:
        print ('{0:4d}   {1: 3.12f}   {2: 3.12f}   {3: 3.1f}'.format(self.Nfeval, theta, alpha, -1*L_m ))
        self.Nfeval += 1
        print('counter is ', counter)

        return(-1*L_m)

                                                #ADJUST UPPER BOUND !
    def optimize(self, inference="beta_likelihood" ,method="Nelder-Mead", niter=1 , temp=0,  param_initial=np.random.uniform(size=2), bnds = ((1e-3, None), (1e-3, None))):
        print("starting optimization")
        mu_initial,sig_initial=param_initial

        if inference=="beta_likelihood":
            likelihood=self.beta_likelihood
            print('Evaluating using the Beta likelihood')
        if inference=='lamperti_likelihood_SDE_approx':
            likelihood=self.lamperti_likelihood_SDE_approx
        if inference=='lamperti_likelihood_moment_approx':
            likelihood=self.lamperti_likelihood_moment_approx

        if method=="Nelder-Mead":
            minimizer_kwargs = {"method":"Nelder-Mead"} #"options":{ 'gtol': myfactr , 'ftol' : myfactr, "maxiter":10 }
            min_param = basinhopping(likelihood, param_initial,T=temp ,minimizer_kwargs=minimizer_kwargs,niter=niter)
        if method=="L-BFGS-B":
            minimizer_kwargs = {"method":"L-BFGS-B" , "options":{ 'gtol': myfactr , 'ftol' : myfactr, "maxiter":10 }}
            min_param = basinhopping(likelihood, param_initial,T=temp ,minimizer_kwargs=minimizer_kwargs,niter=niter)
        if method=="SLSQP":
            minimizer_kwargs = {"method": "SLSQP", "bounds":bnds , "options":{ 'ftol': 1e-10, 'gtol': 1e-10} }
            min_param = basinhopping(likelihood, param_initial,T=temp, minimizer_kwargs=minimizer_kwargs, niter=niter)

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


    def likelihood(self,param):
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

def likelihood_evaluater(theta_array,alpha_array,this_model,inference,path_dir ):

    x = theta_array
    y = alpha_array
    X, Y = np.meshgrid(x, y)
    Z=np.zeros((len(x), len(y)))

    if inference=="beta_likelihood":
        print('Evaluating using the Beta likelihood')
        for i in range(0,len(y)):
            for j in range(0,len(x)):
                Z[j,i] = this_model.beta_likelihood(param=np.array((X[j,i],Y[j,i])))
    if inference=='lamperti_likelihood_SDE_approx':
        for i in range(0,len(y)):
            for j in range(0,len(x)):
                Z[j,i] = this_model.lamperti_likelihood_SDE_approx(param=np.array((X[j,i],Y[j,i])))
    if inference=='lamperti_likelihood_moment_approx':
        for i in range(0,len(y)):
            for j in range(0,len(x)):
                Z[j,i] = this_model.lamperti_likelihood_moment_approx(param=np.array((X[j,i],Y[j,i])))

    Zx, Zy=np.gradient(Z)

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

def empirical_Confidence_Interval_plots(forecast_data_inter,hours,\
    real_in, disct_in,list_forecast_number,\
    dir_path):
    if isinstance(hours, int)==False:
        raise ValueError('hours is not an integer. Please use integer values only.')
    os.makedirs(dir_path + '/'+str(hours)+'hr',exist_ok=True)
    N=disct_in.N
    freq=(N+1)*1 # evert 10 minutes
    #interpolation_points=disct_in.N
    num_forecasts=len(list_forecast_number) # to follow the custom forecast order
    q975=np.zeros((num_forecasts,freq))
    q025=np.zeros((num_forecasts,freq))
    q95=np.zeros((num_forecasts,freq))
    q05=np.zeros((num_forecasts,freq))
    q50=np.zeros((num_forecasts,freq))
    q75=np.zeros((num_forecasts,freq))
    q25=np.zeros((num_forecasts,freq))

    xnew = np.linspace(0,N,freq) #np.linspace(0,N,N+1)
    x = np.linspace(1,N,N)

    theta= real_in.mu #7.8155
    alpha= real_in.sigma #1.072
    j=0
    for k in list_forecast_number: #   #range(0,n_paths): #n_paths
        p=forecast_data_inter[0,k,:N] #obtain cleansed forecast
        inter_1=interpolate.interp1d(x, p, fill_value='extrapolate')
        p=inter_1(xnew)
        d=forecast_data_inter[1,k,:N]
        dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[2,k,0])
        dt=1
        M_test=disct_in.M

        # fig=plt.figure(2,figsize=(10, 4))
        # fig.clf()
        # plt.plot(xnew,p, label='forecast')
        # plt.plot(xnew,d, label='actual production')
        #
        # plt.xlim(1, 73)
        # plt.ylim(-0.1, 1.1)
        # plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24
        #
        # plt.xticks( fontsize = 20);
        # plt.yticks( fontsize = 20);
        # plt.xlabel('Time [hr]',fontsize = 24)
        # plt.ylabel('Power',fontsize = 24)
        # plt.legend( prop={'size': 15})
        # plt.savefig('Forecast_data_'+ str(k)+'.pdf', bbox_inches="tight")


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
    N=disct_in.N
    #interpolation_points=disct_in.N
    num_forecasts=len(list_forecast_number) # to follow the custom forecast order
    q975=np.zeros((num_forecasts,N+1))
    q025=np.zeros((num_forecasts,N+1))
    q95=np.zeros((num_forecasts,N+1))
    q05=np.zeros((num_forecasts,N+1))
    q50=np.zeros((num_forecasts,N+1))
    q75=np.zeros((num_forecasts,N+1))
    q25=np.zeros((num_forecasts,N+1))

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
        dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[2,k,0])
        dt=1
        M_test=disct_in.M

        # fig=plt.figure(2,figsize=(10, 4))
        # fig.clf()
        # plt.plot(xnew,p, label='forecast')
        # plt.plot(xnew,d, label='actual production')
        #
        # plt.xlim(1, 73)
        # plt.ylim(-0.1, 1.1)
        # plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24
        #
        # plt.xticks( fontsize = 20);
        # plt.yticks( fontsize = 20);
        # plt.xlabel('Time [hr]',fontsize = 24)
        # plt.ylabel('Power',fontsize = 24)
        # plt.legend( prop={'size': 15})
        # plt.savefig('Forecast_data_'+ str(k)+'.pdf', bbox_inches="tight")


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

    xnew = np.linspace(0,N,(N+1)*6) #np.linspace(0,N,N+1)
    x = np.linspace(1,N,N)

    theta= real_in.mu #7.8155
    alpha= real_in.sigma #1.072
    j=0
    for k in list_forecast_number: #   #range(0,n_paths): #n_paths
        p=forecast_data_inter[0,k,:N] #obtain cleansed forecast
        inter_1=interpolate.interp1d(x, p, fill_value='extrapolate')
        p=inter_1(xnew)
        #d=forecast_data_inter[1,k,:N]
        dt_object = dtM.datetime.fromtimestamp(forecast_data_inter[2,k,0])
        dt=1
        M_test=disct_in.M

        # fig=plt.figure(2,figsize=(10, 4))
        # fig.clf()
        # plt.plot(xnew,p, label='forecast')
        # plt.plot(xnew,d, label='actual production')
        #
        # plt.xlim(1, 73)
        # plt.ylim(-0.1, 1.1)
        # plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24
        #
        # plt.xticks( fontsize = 20);
        # plt.yticks( fontsize = 20);
        # plt.xlabel('Time [hr]',fontsize = 24)
        # plt.ylabel('Power',fontsize = 24)
        # plt.legend( prop={'size': 15})
        # plt.savefig('Forecast_data_'+ str(k)+'.pdf', bbox_inches="tight")

        theta_adjusted, zero_drift_fixed=theta_adjust(theta,p)

        real_1 = real(theta_adjusted, alpha ) #theta_array

        disct_temp = disct((N+1)*6,dt,M_test)

        X = np.empty((M_test,(N+1)*6));
        X= gen_X_normal_euler_DT_modified(X0=p[0],disct=disct_temp,real=real_1,forecast=p)

        #plotting
        fig=plt.figure(2,figsize=(10, 4))
        fig.clf()

        plt.xlim(1, hours)
        plt.ylim(-0.1, 1.1)

        plt.plot(xnew,p, 'k-', label='forecast',linewidth=5)
        #plt.plot(x,d , 'y-', label='actual production',linewidth=2)

        for i in range(0,len(X)):
            plt.plot(xnew,X[i] , '-o',linewidth=1,markersize=3)
        plt.plot(xnew,X[i] , '-o', label='simulated production',linewidth=1,markersize=3)

        plt.title('{:%d, %b %Y (%H:%M)}'.format(dt_object),fontsize=24)#,fontsize=24

        plt.xticks( fontsize = 20);
        plt.yticks( fontsize = 20);
        plt.xlabel('Time [hr]',fontsize = 24)
        plt.ylabel('Power',fontsize = 24)
        plt.legend( prop={'size': 20})
        plt.savefig(dir_path+'/'+ str(hours) +'hr/'+str(k)+'.pdf', bbox_inches="tight")

    file_object  = open(dir_path+'/parameter.info', 'w')
    note='Simulation of ' + str(num_forecasts)+' forecasts spanning '\
    + str(N) +' hours.'+ '\n' +'Parameters used are theta= '\
    +str(theta) + ' and alpha= ' + str(alpha) + '.Generated by '+str(M_test) \
    + ' simulations.'
    file_object.write(note)
    file_object.close()
    np.save(dir_path + '/forecast_list', list_forecast_number)


def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--', label='slope='+str(slope))
