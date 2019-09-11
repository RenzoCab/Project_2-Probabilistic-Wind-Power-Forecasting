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
from multiprocessing import Process
import time

def f(name, sec):
    info('function f')
    time.sleep(sec)
    print ('hello', name)
    L=sec

def info(title):
    print (title)
    print ('module name:', __name__)
    if hasattr(os, 'getppid'):  # only available on Unix
        print('parent process:', os.getppid())
    print('process id:', os.getpid())

if __name__ == '__main__':
    p1 = Process(target=f, args=('bob',0))
    p2 = Process(target=f, args=('tom',5))
    p3 = Process(target=f, args=('goerge',10))
    p1.start()
    p2.start()
    p3.start()
    p1.join()


def moment_compute(dN, X_prev, X_next, p_prev, theta_current, alpha_current ):
    m_1=X_prev*np.exp(- dN*theta_current)
    m_2=  (X_prev**2 + 2*dN*( X_prev*(alpha_current*theta_current*p_prev*(1-p_prev)*(1-2*p_prev )) + \
                 alpha_current*theta_current*p_prev**2*(1-p_prev)**2) ) /(  1+ dN*2*(theta_current+alpha_current*theta_current*p_prev*(1-p_prev) ))
    a=m_1;
    b=m_2- m_1**2;

    if b==0: b=1e-16
    if b<0:
        b=1e-16

    beta_param_alpha= - ( (1+a)*(a**2 +b -1)    )/(2*b)
    beta_param_beta= ( (a-1)*(a**2 + b -1)  )  /(2*b)

    L_n_current=(beta_param_alpha-1 )*np.log(  (X_next+1)/2 ) +\
     (beta_param_beta-1)*np.log(1-( X_next +1)/2 )-\
     scipy.special.betaln(beta_param_alpha, beta_param_beta)

    return(L_n_current)



moment_compute(0.1,0.4,0.6, 10,0.1)

from multiprocessing import Process, Queue

class Multiprocessor():

    def __init__(self):
        self.processes = []
        self.queue = Queue()

    @staticmethod
    def _wrapper(func, queue, args, kwargs):
        ret = func(*args, **kwargs)
        queue.put(ret)

    def run(self, func, *args, **kwargs):
        args2 = [func, self.queue, args, kwargs]
        p = Process(target=self._wrapper, args=args2)
        self.processes.append(p)
        p.start()

    def wait(self):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        for p in self.processes:
            p.join()
        return rets

# tester
if __name__ == "__main__":
    mp = Multiprocessor()
    num_proc = 10
    for i in range(num_proc): # queue up multiple tasks running `sum`
        mp.run(sum, [2, i])
    ret = mp.wait() # get all results
    print(ret)
