import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats
from scipy.optimize import minimize
from scipy.stats import multivariate_normal
from mpl_toolkits.mplot3d import Axes3D
import operator
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import math
import random
from scipy import interpolate

class disct:
  def __init__(self, N, dt, M ):
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


class model:
    def __init__(self, disct, data,ic=0 ):
        self.disct = disct
        self.data = data
        self.ic = ic
        self.optimal=optimal()

    def likelihood(self,param):
        #data
        X=self.data

        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        #input parameters
        a,b = param;

        L_n=0
        L_m=0

        for j in range(0,M):
            L_m = L_m + L_n
            for i in range(1,N):
                L_n = L_n  - (X[j,i]- X[j,i-1] - a*dt )**2/(2*dt*b**2) -0.5*np.log(2*math.pi*b**2*dt)

        return(-1*L_m)

    def renew(self,real):

        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        X = np.zeros((M,N))

        #p=0.5 + 0.5*np.sin(np.linspace(0,2*math.pi,N))

        dW=np.sqrt(dt)*np.random.normal(0, 1, N)

        X0=self.ic

        for j in range(0,M):

            dW=np.sqrt(dt)*np.random.normal(0, 1, N)

            for i in range(1,N):
                X[j,0]=X0
                b= real.sigma      #sig_real*X[i-1]*(1-X[i-1])
                a= real.mu       #mu_real*(0.5 - X[i-1])
                X[j,i]= X[j,i-1]+ a*dt+ b*dW[i-1]
        self.data = X

        return()

    def optimize(self, param_initial=np.random.uniform(size=2)*10, bnds = ((None, None), (1e-5, None))):

        mu_initial,sig_initial=param_initial

        min_param=minimize(self.likelihood,param_initial,bounds=bnds);min_param

        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)

class model_in_box:
    def __init__(self, disct, data,ic=0 ):
        self.disct = disct
        self.data = data
        self.ic = ic
        self.optimal=optimal()

    def likelihood(self,param):
        #data
        X=self.data

        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        #input parameters
        mu,sigma = param;

        L_n=0
        L_m=0

        for j in range(0,M):
            L_m = L_m + L_n
            for i in range(1,N):
                b= sigma*X[j,i-1]*(1-X[j,i-1])
                a= mu*(0.5 - X[j,i-1])
                L_n = L_n  - (X[j,i]- X[j,i-1] - a*dt )**2/(2*dt*b**2) -0.5*np.log(2*math.pi*b**2*dt)

        return(-1*L_m)

    def renew(self,real):

        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        X = np.zeros((M,N))

        #p=0.5 + 0.5*np.sin(np.linspace(0,2*math.pi,N))

        dW=np.sqrt(dt)*np.random.normal(0, 1, N)

        X0=self.ic

        for j in range(0,M):

            dW=np.sqrt(dt)*np.random.normal(0, 1, N)

            for i in range(1,N):
                X[j,0]=X0
                b= real.sigma*X[j,i-1]*(1-X[j,i-1])
                a= real.mu*(0.5 - X[j,i-1])
                X[j,i]= X[j,i-1]+ a*dt+ b*dW[i-1]
            self.data = X

        return()

    def optimize(self, param_initial=np.random.uniform(size=2), bnds = ((None, None), (1e-5, None))):

        mu_initial,sig_initial=param_initial

        min_param=minimize(self.likelihood,param_initial,bounds=bnds);min_param

        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)

class model_in_box_sine:
    def __init__(self, disct, data,ic=0 ):
        self.disct = disct
        self.data = data
        self.ic = ic
        self.optimal=optimal()

    def likelihood(self,param):
        #data
        X=self.data

        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        #input parameters
        mu,sigma = param;

        L_n=0
        L_m=0

        for j in range(0,M):
            L_m = L_m + L_n
            for i in range(1,N):
                b= sigma*X[j,i-1]*(1-X[j,i-1])
                a= mu*( 0.5 + 0.5*np.sin(i/N * 2*math.pi) - X[j,i-1])
                L_n = L_n  - (X[j,i]- X[j,i-1] - a*dt )**2/(2*dt*b**2) -0.5*np.log(2*math.pi*b**2*dt)

        return(-1*L_m)

    def renew(self,real):

        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        X = np.zeros((M,N))

        #p=0.5 + 0.5*np.sin(np.linspace(0,2*math.pi,N))

        dW=np.sqrt(dt)*np.random.normal(0, 1, N)

        X0=self.ic

        for j in range(0,M):

            dW=np.sqrt(dt)*np.random.normal(0, 1, N)

            for i in range(1,N):
                X[j,0]=X0
                b= real.sigma*X[j,i-1]*(1-X[j,i-1])
                a= real.mu*( 0.5 + 0.5*np.sin(i/N * 2*math.pi) - X[j,i-1])
                X[j,i]= X[j,i-1]+ a*dt+ b*dW[i-1]
            self.data = X

        return()

    def optimize(self, param_initial=np.random.uniform(size=2), bnds = ((None, None), (1e-5, None))):

        mu_initial,sig_initial=param_initial

        min_param=minimize(self.likelihood,param_initial,bounds=bnds);min_param

        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)

class model_moments:
    def __init__(self, disct, data,ic=0 ):
        self.disct = disct
        self.data = data
        self.ic = ic
        self.optimal=optimal()

    def likelihood(self,param):
        #data
        X=self.data

        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        #input parameters
        theta,alpha = param;

        L_n=0;
        L_m=0;
        a=0;b=0;
        for j in range(0,M):
            L_m = L_m + L_n
            for i in range(1,N):
                b =   alpha/(16*(1+ alpha/4)) * (1- np.exp( -2*dt*theta*(1+ alpha/4)  )) #variance not std
                #b=np.sqrt( 2*alpha*theta*(1/4)*(X[j,i-1] + 1/2 )*(1-X[j,i-1] + 1/2 )  )
                a =  X[j,i-1]*np.exp(- dt*theta)   #mu*( 0.5 + 0.5*np.sin(i/N * 2*math.pi) - X[j,i-1])
                L_n = L_n  - ( X[j,i] - a   )**2/(2*b**2) -0.5*np.log(2*math.pi*b**2) #  removed dt from both terms

        return(-1*L_m)

    def renew(self,real):

        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        theta=real.mu
        alpha=real.sigma

        X = np.zeros((M,N))
        X0=0
        dW=np.random.normal(0, 1, N) #removed sqrt(dt)

        for j in range(0,M):

            dW=dt*np.random.normal(0, 1, N) #removed sqrt(dt)

            for i in range(1,N):
                X[j,0]=X0

                b = np.sqrt( 2*alpha*theta*(1/4)*(X[j,i-1] + 1/2 )*(1-X[j,i-1] + 1/2 )  )

                X[j,i]= (X[j,i-1] +  b*dW[i-1] )/(1+  theta*dt )

            self.data = X

        return()

    def optimize(self, param_initial=np.random.uniform(size=2), bnds = ((None, None), (1e-5, None))):

        mu_initial,sig_initial=param_initial
        #L-BFGS-B
        myfactr = 1
        min_param=minimize(self.likelihood,param_initial,bounds=bnds, method='L-BFGS-B', options={ 'ftol' : myfactr * np.finfo(float).eps });min_param

        #minimizer_kwargs = {"method": "SLSQP", "bounds":bnds} #, "options":{ 'ftol': 1e-10, 'gtol': 1e-10} } #'gtol': 1e-9
        #min_param = basinhopping(self.likelihood, x0=param_initial, minimizer_kwargs=minimizer_kwargs) # 'ftol': 1e-25,'gtol': 1e-25
        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x, min_param.message)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)

class model_beta_moments:
    def __init__(self, disct, data,ic=0 ):
        self.disct = disct
        self.data = data
        self.ic = ic
        self.optimal=optimal()

    def likelihood(self,param):
        #data
        X=self.data

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
        for j in range(0,M):
            L_m = L_m + L_n
            for i in range(1,N):

                b =   alpha/(16*(1+ alpha/4)) * (1- np.exp( -2*dt*theta*(1+ alpha/4)  )) #variance not std
                a =  X[j,i-1]*np.exp(- dt*theta)

                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)   #lambd*( (a+1)  /2  ) #( (1-a)/b - 1/a )*a**2
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)  #lambd*( (1-a) / 2)

                #if beta_param_alpha <0 or beta_param_beta <0:
                #    print('WARNING: Negative Shape parameters !')
                #print(scipy.special.beta(beta_param_alpha, beta_param_beta))

                L_n = L_n +  (beta_param_alpha-1 )*np.log(  (X[j,i]+1)/2 ) + (beta_param_beta-1)*np.log(1-( X[j,i] +1)/2 ) -scipy.special.betaln(beta_param_alpha, beta_param_beta)

                #- np.log( scipy.special.beta(beta_param_alpha, beta_param_beta) ) #-np.log(2)

                #L_n = L_n  - ( X[j,i] - a   )**2/(2*b**2) -0.5*np.log(2*math.pi*b**2)

        return(-1*L_m)

    def renew(self,real):

        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        theta=real.mu
        alpha=real.sigma

        X = np.zeros((M,N))
        X0=0
        dW=np.random.normal(0, 1, N) #removed sqrt(dt)

        for j in range(0,M):

            dW=dt*np.random.normal(0, 1, N) #removed sqrt(dt)

            for i in range(1,N):
                X[j,0]=X0

                b = np.sqrt( 2*alpha*theta*(1/4)*(X[j,i-1] + 1/2 )*(1-X[j,i-1] + 1/2 )  )

                X[j,i]= (X[j,i-1] +  b*dW[i-1] )/(1+  theta*dt )

            self.data = X

        return()

    def optimize(self, param_initial=np.random.uniform(size=2), bnds = ((0.1, None), (0.1, None))):

        mu_initial,sig_initial=param_initial
        #L-BFGS-B
        myfactr = 1
        min_param=minimize(self.likelihood,param_initial,bounds=bnds, method='L-BFGS-B', options={ 'ftol' : myfactr * np.finfo(float).eps });min_param
        #minimizer_kwargs = {"method": "SLSQP", "bounds":bnds} #, "options":{ 'ftol': 1e-10, 'gtol': 1e-10} } #'gtol': 1e-9
        #min_param = basinhopping(self.likelihood, x0=param_initial, minimizer_kwargs=minimizer_kwargs) # 'ftol': 1e-25,'gtol': 1e-25
        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x, min_param.message)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)

class model_beta_flex:
    def __init__(self, disct, data,forecast ):
        self.disct = disct
        self.data = data
        #self.ic = ic
        self.forecast = forecast
        self.optimal=optimal()

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
            for i in range(0,N-1): #start from 1 or zero

                m_1=X[j,i]*np.exp(- dN*theta)
                #
                # m_2=  X[j,i]**2 + dN*2*(-X[j,i]**2*(theta+alpha*theta*p[i]*(1-p[i]) ) + \
                #                 X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
                #                     alpha*theta*p[i]**2*(1-p[i])**2)

                m_2=  (X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta*p[j,i]*(1-p[j,i])*(1-2*p[j,i] )) + \
                            alpha*theta*p[j,i]**2*(1-p[j,i])**2) ) /(  1+ dN*2*(theta+alpha*theta*p[j,i]*(1-p[j,i]) ))
                a=m_1
                b=m_2 - m_1**2

                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                L_n = L_n +  (beta_param_alpha-1 )*np.log(  (X[j,i+1]+1)/2 ) +\
                 (beta_param_beta-1)*np.log(1-( X[j,i+1] +1)/2 )-\
                 scipy.special.betaln(beta_param_alpha, beta_param_beta)

        return(-1*L_m)

    def grad(self,alpha, beta ):

        V=self.data
        N=self.disct.N
        M=self.disct.M
        #input parameters

        s1=np.sum(np.log((V+1)/2), axis=None)
        s2=np.sum(np.log( 1-  (V+1)/2), axis=None)

        digamma_1=-M*N*scipy.special.digamma(alpha)
        digamma_2=-M*N*scipy.special.digamma(beta)
        digamma_common=-M*N*scipy.special.digamma(alpha+beta)

        grad_eval=np.array(( s1 + digamma_1 +digamma_common  , s2  + digamma_2 +digamma_common ))

        return(grad_eval) #not useful because, needs grad in terms of theta and alpha not beta parameter


    def renew(self,real):

        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        theta=real.mu
        alpha=real.sigma

        X = np.zeros((M,N))
        X0=0
        dW=np.random.normal(0, 1, N) #removed sqrt(dt)

        for j in range(0,M):

            dW=dt*np.random.normal(0, 1, N) #removed sqrt(dt)

            for i in range(1,N):
                X[j,0]=X0

                b = np.sqrt( 2*alpha*theta*(1/4)*(X[j,i-1] + 1/2 )*(1-X[j,i-1] + 1/2 )  )

                X[j,i]= (X[j,i-1] +  b*dW[i-1] )/(1+  theta*dt )

            self.data = X

        return()
                                                #ADJUST UPPER BOUND !
    def optimize(self, param_initial=np.random.uniform(size=2), bnds = ((1e-2, None), (1e-2, None))):

        mu_initial,sig_initial=param_initial
        #L-BFGS-B
        #myfactr = 1

        min_param=minimize(self.likelihood,param_initial,bounds=bnds, method='L-BFGS-B', options={ 'gtol': 1e-19, 'ftol' : 1e-19 });min_param
        #, options={ 'gtol': myfactr * np.finfo(float).eps, 'ftol' : myfactr * np.finfo(float).eps });min_param

        #minimizer_kwargs = {"method": "SLSQP", "bounds":bnds} #, "options":{ 'ftol': 1e-10, 'gtol': 1e-10} } #'gtol': 1e-9
        #min_param = basinhopping(self.likelihood, x0=param_initial, minimizer_kwargs=minimizer_kwargs) # 'ftol': 1e-25,'gtol': 1e-25
        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x, min_param.message)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)

class model_beta_data:
    def __init__(self, disct, data,forecast ):
        self.disct = disct
        self.data = data
        #self.shape_param_alpha=np.empty_like(data)
        #self.shape_param_beta=np.empty_like(data)
        #self.ic = ic
        self.forecast = forecast
        self.optimal=optimal()

    def compute_shape_parameters(self, param): #needs shape parameters to run
        #data
        X=self.data
        p=self.forecast
        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M
        theta,alpha = param;
        #print( theta, alpha)
        a=0;b=0;
        eps=np.finfo(float).eps;
        m_1=0
        m_2=0
        dN=1/N #set N=72


        for j in range(0,M):
            m_1=0
            m_2=0
            interpolation_points=200 #add to inputs
            dx=1/interpolation_points
            x = np.arange(1,N+1,1)
            y = X[j,:]
            tck = interpolate.splrep(x, y, s=0)
            xnew = np.arange(1,72+1,dx)
            p_inter = interpolate.splev(xnew, tck, der=0)

            for i in range(0,N-1): #start from 1 or zero
                m_1=X[j,i]*np.exp(- dN*theta)

                m_2=X[j,i]**2
                for k in range(i*interpolation_points, interpolation_points*(i+1)):
                    p=p_inter[k]
                    m_2 = m_2 + dN*dx*2*(-X[j,i]**2*(theta+alpha*theta*p*(1-p) ) + \
                            X[j,i]*(alpha*theta*p*(1-p)*(1-2*p )) + \
                                alpha*theta*p**2*(1-p)**2)

                #print('m1=',m_1)
                #print('m2=',m_2)
                a=m_1
                b=m_2 - m_1**2

                #print('mean=',a)
                #print('variance=',b)

                self.shape_param_alpha[j,i]= - ( (1+a)*(a**2 +b -1)    )/(2*b )
                self.shape_param_beta[j,i]= ( (a-1)*(a**2 + b -1)  )  /(2*b)

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
        dN=1/N #set N=72
        for j in range(0,M):
            m_1=0
            m_2=0
            L_m = L_m + L_n

            interpolation_points=10 #add to inputs
            dx=1/interpolation_points
            x = np.arange(1,N+1,1)
            y = self.forecast[j,:] #X[j,:] ###### CHECK !
            tck = interpolate.splrep(x, y, s=0)
            xnew = np.arange(1,72+1,dx)
            p_inter = interpolate.splev(xnew, tck, der=0)

            for i in range(0,N-1): #start from 1 or zero
                m_1=X[j,i]*np.exp(- dN*theta)
                m_2=X[j,i]**2
                for k in range(i*interpolation_points, interpolation_points*(i+1)):
                    p=p_inter[k]
                    # m_2 = m_2 + dN*dx*2*(-X[j,i]**2*(theta+alpha*theta*p*(1-p) ) + \
                    #         X[j,i]*(alpha*theta*p*(1-p)*(1-2*p )) + \
                    #             alpha*theta*p**2*(1-p)**2)
                    m_2=  (X[j,i]**2 + 2*dN*( X[j,i]*(alpha*theta*p[i]*(1-p[i])*(1-2*p[i] )) + \
            alpha*theta*p[i]**2*(1-p[i])**2) ) /(  1+ dN*2*(theta+alpha*theta*p[i]*(1-p[i]) ))
                #print('m1=',m_1)
                #print('m2=',m_2)
                a=m_1
                b=m_2 - m_1**2

                #print('mean=',a)
                #print('variance=',b)

                beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b )
                beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

                #print(beta_param_alpha, beta_param_beta)

                #if beta_param_alpha <0 or beta_param_beta <0:
                #    print('WARNING: Negative Shape parameters !')
                #print(scipy.special.beta(beta_param_alpha, beta_param_beta))

                L_n = L_n +  (beta_param_alpha-1 )*np.log(  (X[j,i+1]+1)/2 ) +\
                 (beta_param_beta-1)*np.log(1-( X[j,i+1] +1)/2 )-\
                 scipy.special.betaln(beta_param_alpha, beta_param_beta)

        return(-1*L_m)

    def optimize(self, param_initial=np.random.uniform(size=2), bnds = ((1e-3, None), (1e-3, None))):

        mu_initial,sig_initial=param_initial
        #L-BFGS-B
        #myfactr = 1
        min_param=minimize(self.likelihood,param_initial,bounds=bnds, method='L-BFGS-B', options={ 'gtol': 1e-16, 'ftol' : 1e-16 });min_param
        #, options={ 'ftol' : myfactr * np.finfo(float).eps }
        #minimizer_kwargs = {"method": "SLSQP", "bounds":bnds} #, "options":{ 'ftol': 1e-10, 'gtol': 1e-10} } #'gtol': 1e-9
        #min_param = basinhopping(self.likelihood, x0=param_initial, minimizer_kwargs=minimizer_kwargs) # 'ftol': 1e-25,'gtol': 1e-25
        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x, min_param.message)



class model_moments_check:
    def __init__(self, disct, data,ic=0 ):
        self.disct = disct
        self.data = data
        self.ic = ic
        self.optimal=optimal()

    def likelihood(self,param):
        #data
        X=self.data

        #discretization object
        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        #input parameters
        theta,alpha = param;

        L_n=0;
        L_m=0;
        a=0;b=0;
        for j in range(0,M):
            L_m = L_m + L_n
            for i in range(1,N):
                #b =   alpha/(16*(1+ alpha/4)) * (1- np.exp( -2*dt*theta*(1+ alpha/4)  )) #variance not std
                b=np.sqrt( 2*alpha*theta*(1/4)*(X[j,i] + 1/2 )*(1/2 - X[j,i] )  )
                a =  - X[j,i]*theta  #mu*( 0.5 + 0.5*np.sin(i/N * 2*math.pi) - X[j,i-1])
                L_n = L_n  - ( X[j,i] - X[j,i-1] + a*dt  )**2/(2*b**2*dt) -0.5*np.log(2*math.pi*b**2*dt) #  -  X[j,i-1] removed
        return(-1*L_m)

    def renew(self,real):

        N=self.disct.N
        dt=self.disct.dt
        M=self.disct.M

        theta=real.mu
        alpha=real.sigma

        X = np.zeros((M,N))
        X0=0
        dW=np.random.normal(0, 1, N) #removed sqrt(dt)

        for j in range(0,M):

            dW=dt*np.random.normal(0, 1, N) #removed sqrt(dt)

            for i in range(1,N):
                X[j,0]=X0

                b = np.sqrt( 2*alpha*theta*(1/4)*(X[j,i-1] + 1/2 )*(1-X[j,i-1] + 1/2 )  )

                X[j,i]= (X[j,i-1] +  b*dW[i-1] )/(1+  theta*dt )

            self.data = X

        return()

    def optimize(self, param_initial=np.random.uniform(size=2), bnds = ((None, None), (1e-5, None))):

        mu_initial,sig_initial=param_initial
        #L-BFGS-B
        min_param=minimize(self.likelihood,param_initial,bounds=bnds, method='L-BFGS-B', options={ 'ftol': 1e-25, 'gtol': 1e-25});min_param
        #minimizer_kwargs = {"method": "SLSQP", "bounds":bnds} #, "options":{ 'ftol': 1e-10, 'gtol': 1e-10} }
        #min_param = basinhopping(self.likelihood, x0=param_initial, minimizer_kwargs=minimizer_kwargs)
        self.optimal.mu= min_param.x[0]
        self.optimal.sigma= min_param.x[1]

        return(min_param.x)

    def get_error(self,real):
        err_mu = np.abs( self.optimal.mu - real.mu )/real.mu
        err_sigma = np.abs( self.optimal.sigma - real.sigma )/real.sigma
        return(err_mu,err_sigma)
###########################################################################

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
    max_tries=200
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



def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--', label='slope='+str(slope))
