#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:20:37 2022

@author: sampasmann
"""

import numpy as np
#import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import sympy as sym


def shift_samples(LB, RB, x_values, plot=False):
    """
    Parameters
    ----------
    LB : left bound
    RB : right bound
    samps : samples to shift (use np.random.uniform() to generate samps)

    Returns
    -------
    x : moment matched (shifted) samples
    
    This algorithm takes the bounds from a definite integral and calculates
    the corresponding logistic function, inverse of the logistic function, and
    derivative of the inverse (the PDF).

    """
    
    
    #E = (RB-LB)/2 + LB
    E = (RB+LB)/2
    if (LB<0):
        LB *= -1
    
    # logit function
    g = lambda x: np.log((LB+x)/(RB-x))
    # inverse logit function
    ginv = lambda y: (RB*np.exp(y)-LB)/(np.exp(y)+1) 
    # pdf of logit function (derivative of inverse)
    pdf = lambda y: np.exp(y)*(LB+RB)/(np.exp(y)+1)**2 
    # inverse pdf (logit function)
    pinv = lambda x: g(x)  
    
    # Expected value of shifted samples in transformed space
    samps = pinv(x_values)
    f = lambda d: (np.mean(ginv(samps+d))-E) 
    a = -10*(RB+LB)
    b = 10*(RB+LB)
    ds = bisection(f,a,b,tol=1e-8)
    yhat = samps + ds
    x = ginv(yhat)
    
    x_values_avg = x_values.sum()/x_values.size
    x_avg = x.sum()/x.size
    print("Prior Shift Average:", x_values_avg)
    print("Post Shift Average:", x_avg)
    
    if (plot):
    
        plt.figure(dpi=300, figsize=(18,6))
        plt.subplot(131)
        xx = np.linspace(-RB,RB)
        plt.plot(xx,ginv(xx))
        plt.title('Inverse Logit Function')

        
        plt.subplot(132)
        xx = np.linspace(-7.5,7.5)
        plt.plot(xx, pdf(xx))
        plt.title('Derivative of Inverse (PDF)')

        
        plt.subplot(133)
        #plt.hist(g(x_values),bins=20)
        plt.hist(x,bins=10)
        plt.title('Shifted Samples')

    
    return x

def shift_samples_sympy(LB, RB, x_values, plot=False):
    
    E = (RB+LB)/2
    if (LB<0):
        LB *= -1
    
    x, y = sym.symbols("x, y")
    # logistic function
    g = sym.Eq(y, sym.log((LB+x)/(RB-x)))
    # inverse of logistic function
    ginv = sym.solve(g, x)[0].simplify()
    g = g.rhs # weird thing with sympy where it stored the lhs and rhs, need
    # g to just be the rhs 
    # derivative of inverse = Probability Density Function
    pdf = sym.diff(ginv).simplify()
    
    
    # lambdify functions
    g = sym.lambdify(x, g, "numpy")
    ginv = sym.lambdify(y, ginv, "numpy")
    pdf = sym.lambdify(y, pdf, "numpy")
    
    
    samps = g(x_values)
    f = lambda d: (np.mean(ginv(samps+d))-E) 
    ds = bisection(f,-2*(RB+LB),2*(RB+LB),tol=1e-8)
    yhat = samps + ds
    shifted_samps = ginv(yhat)
    
    if (plot):
        
        plt.figure(dpi=300, figsize=(18,6))
        plt.subplot(131)
        xx = np.linspace(LB,RB)
        plt.plot(xx,g(xx))
        plt.title('Logistic Function')
        plt.grid()
        
    
        plt.subplot(132)
        xx = np.linspace(-7.5,7.5)
        plt.plot(xx, pdf(xx))
        plt.title('Derivative of Inverse (PDF)')
        plt.grid()
        
        plt.subplot(133)
        plt.hist(shifted_samps,bins=20)
        plt.title('Shifted Samples')
    
    return shifted_samps


def get_answer(f,LB,RB):
   
    x = sym.symbols("x")
    F = sym.integrate(f,(x,LB,RB))
    A = F.evalf(subs={x: RB}) - F.evalf(subs={x: LB})
    
    return F



def bisection(f,a,b,tol=1e-6):
    
    assert (f(a)*f(b)<0)

    fa = f(a)
    fb = f(b)
    while (b-a > tol):
        mid = (a+b)*0.5
        fm = f(mid)
        if (fa*fm < 0):
            b = mid
            fb = fm
        else:
            a = mid
            fa = fm
    return 0.5*(a+b)



def MC_plotting(F,f,LB,RB,Ns,A,sols,sols_shift):
    """
    Plots:
        - function and area under the curve given bounds
        - standard deviation of estimate
        - mean absolute error of estimate
    """
    fig, ax = plt.subplots(1,3,figsize=(18,6), dpi=500)
    
    ########################################
    # Integral Plot
    ########################################
    a, b = LB, RB  # integral limits
    x = np.linspace(a+0.2*a, b+0.2*b)
    y = f(x)
    
    ax[0].plot(x, y, 'r', linewidth=1)
    
    # Make the shaded region
    ix = np.linspace(a, b)
    iy = f(ix)
    verts = [(a, 0), *zip(ix, iy), (b, 0)]
    poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
    ax[0].add_patch(poly)
    latex = sym.latex(F)
    ax[0].text(0.5 * (a + b), 0.1, r"$\int_{{{0}}}^{{{1}}} {2} dx$".format(round(a,2),round(b,2),latex),
            horizontalalignment='center', fontsize=26)
    
    #ax[0].spines.right.set_visible(False)
    #ax[0].spines.top.set_visible(False)
    ax[0].xaxis.set_ticks_position('bottom')
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")
    
    ax[0].set_xticks((a, b))
    ax[0].set_xticklabels(('{0}'.format(round(a,2)), '{0}'.format(round(b,2))))
    ax[0].set_yticks([])
    ax[0].set_title('Integral',fontdict={'fontsize': 16, 'fontweight': 'medium'})

    ########################################
    ########################################
    
    vrs = np.std(sols,axis=1)
    vrs_shift = np.std(sols_shift,axis=1)
    ax[1].loglog(Ns,vrs,"o-")
    ax[1].loglog(Ns,vrs_shift,"*-")
    #ax[1].spines.right.set_visible(False)
    #ax[1].spines.top.set_visible(False)
    ax[1].set_xlabel("Number of Samples")
    ax[1].set_ylabel("Standard Deviation of Estimate")
    ax[1].set_title('Standard Deviation',fontdict={'fontsize': 16, 'fontweight': 'medium'})
    #plt.grid()
    
    merr = np.mean(np.abs(sols-A),axis=1)
    merr_shift = np.mean(np.abs(sols_shift-A),axis=1)
    ax[2].loglog(Ns,merr,"o-",label='Standard')
    ax[2].loglog(Ns,merr_shift,"*-",label='Shifted')
    #ax[2].spines.right.set_visible(False)
    #ax[2].spines.top.set_visible(False)
    #plt.plot(Ns, np.array(Ns)**-0.5)
    ax[2].set_xlabel("Number of Samples")
    ax[2].set_ylabel("Mean Absolute Error of Estimate")
    ax[2].set_title('Error',fontdict={'fontsize': 16, 'fontweight': 'medium'})
    #plt.grid()
    ax[2].legend()
    
    return



def monte_carlo_test(F, A, LB, RB, Ns=[1000], times=30, sympy=False):
    """
    Parameters
    ----------
    F : SymPy lambda function to be integrated with Monte Carlo.
        F may only have one variable, "x" and may only use sympy functions.
    A : Numerical answer to definite integral
    LB : left-bound
    RB : right-bound
    Ns : Array of number of particles to use The default is [1000].
         If length of array > 1 plots are produced.
    times: number of times each particle count is run.

    Returns
    -------
    None.

    """
    x = sym.symbols("x")
    f = sym.lambdify(x, F, "numpy")
    
    Ns = np.array(Ns)
    d = 0
    sols = np.zeros([len(Ns),times])
    sols_shift = np.zeros([len(Ns),times])
    
    for i in range(len(Ns)):
        for t in range(times):
            # Normal Solve
            samps = np.random.uniform(LB, RB, size=int(Ns[i]))
            sols[i,t] = (RB-LB)*np.sum(f(samps))/Ns[i]
            # Shifted Solve
            if (sympy):
                samps = shift_samples_sympy(LB, RB, samps)
            else:
                samps = shift_samples(LB, RB, samps)
            sols_shift[i,t] = (RB-LB)*np.sum(f(samps))/Ns[i]
            
    if (Ns.size > 1):
        MC_plotting(F,f,LB,RB,Ns,A,sols,sols_shift)
    else:
        print('Standard Error: ', np.mean(np.abs(sols-A)))
        print('Shifted Error: ', np.mean(np.abs(sols_shift-A)))
    
    return 
    
    
    
    
    
    