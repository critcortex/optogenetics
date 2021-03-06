import numpy as np
from scipy.special import erf
from scipy.optimize import curve_fit
from scipy.stats import norm
import pylab
import math

lambdas = np.array([390.,410.,430.,445.,470.,500.,520.,540.,560.,585.,630.])
smoothx = np.linspace(lambdas[0],lambdas[-1], num=len(lambdas))
smooth4x = np.linspace(lambdas[0],lambdas[-1], num=4*len(lambdas))

i_chr2 = np.array([700.,850.,920.,880.,880.,800.,600.,450.,200.,0.,0.])
i_nphr = np.array([175.,275.,300.,300.,400.,600.,800.,950.,1050.,1100.,450.])

i_archt = np.array([0.17,0.17,0.35,0.55,0.65,0.95,1.,0.8,0.35,0.05])
lambdas_archt = np.array([400,435,460,485,520,540,575,600,630,660])
smooth_archtx = np.linspace(lambdas_archt[0],lambdas_archt[-1], num=len(lambdas_archt))

pylab.close('all')

def skew_func(x,a=0.,b=1.,c=1.,d=0.):#*p):
    #lognormal
    #mu = a 
    #sigma= b
    #return 1./(x*sigma*math.sqrt(2*np.pi))*np.exp((-(np.log(x)-mu)**2)/(2*sigma**2))
    
    #skewed
    xscale = (x-b)/c
    phi_1 = 1./(math.sqrt(2*np.pi))*np.exp(-(xscale**2)/2)
    phi_2 = 0.5*(1+erf((a*xscale)/math.sqrt(2)))
    return 2./c*phi_1*phi_2+d
    



def rayleigh_func(x,A,sigma):
    return A*x/sigma*sigma*(np.exp(-x**2/(2*sigma*sigma)))

def gauss_func(x, *p):
    A = p[0]
    mu = p[1]
    sigma = p[2]
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def try_fit():

    p0_chr2 = [0.5,450,1./800,1]
    p0_hr = [-0.5,575,1./1000,1]
    
    
    pylab.figure()
    xs = np.arange(0.,5.1,0.01)
    p0_chr2 = [-0.5,2.,1./2.,1]
    p0_hr = [-1,1.,1.,1]
    th_xs1 = skew_func(xs,p0_chr2[0],p0_chr2[1],p0_chr2[2],p0_chr2[3])
    th_xs2 = skew_func(xs,p0_hr[0],p0_hr[1],p0_hr[2],p0_hr[3])
    ls = '-'
    pylab.plot(xs,th_xs1,c='b',ls=ls)
    pylab.plot(xs,th_xs2,c='y',ls=ls)
    
    p0_ch2_rayleigh = [1.,2.]
    p0_nphr_rayleigh = [2.,0.5]
    th_xs1 = rayleigh_func(xs,p0_ch2_rayleigh[0],p0_ch2_rayleigh[1])
    th_xs2 = rayleigh_func(xs,p0_nphr_rayleigh[0],p0_nphr_rayleigh[1])
    ls = '--'
    pylab.plot(xs,th_xs1,c='b',ls=ls)
    pylab.plot(xs,th_xs2,c='y',ls=ls)
    
    p0_ch2_gauss = [1.,2.,0.2]
    p0_nphr_gauss = [2.,0.5,2.]
    th_xs1 = gauss_func(xs,*p0_ch2_gauss)
    th_xs2 = gauss_func(xs,*p0_nphr_gauss)
    ls = ':'
    pylab.plot(xs,th_xs1,c='b',ls=ls)
    pylab.plot(xs,th_xs2,c='y',ls=ls)
    
    
    """
    th_cr = skew_func(smoothx,p0_chr2[0],p0_chr2[1],p0_chr2[2],p0_chr2[3])
    th_hr = skew_func(smoothx,p0_hr[0],p0_hr[1],p0_hr[2],p0_hr[3])
    """
    """
    pylab.figure()
    pylab.scatter(smoothx,i_chr2,c='b',s=40,marker='o')
    pylab.scatter(smoothx,i_nphr,c='y',s=40,marker='s')
    pylab.plot(smoothx,th_cr,ls='-',c='b',lw=2)
    pylab.plot(smoothx,th_hr,ls='--',c='y',lw=2)
    """
    """
    popt_ch2, pcov_ch2 = curve_fit(func, smoothx, i_chr2,p0=p0_chr2)
    popt_hr, pcov_hr = curve_fit(func, smoothx, i_nphr,p0=p0_hr)
    print popt_ch2
    print popt_hr
    
    pylab.scatter(smoothx,i_chr2,c='b',s=40,marker='o')
    pylab.scatter(smoothx,i_nphr,c='y',s=40,marker='s')
    th_ichr = func(smoothx,popt_ch2[0],popt_ch2[1],popt_ch2[2],popt_ch2[3])
    pylab.plot(smoothx,th_ichr,ls='-',c='b',lw=4)
    th_ihr = func(smoothx,popt_hr[0],popt_hr[1],popt_hr[2],popt_hr[3])
    pylab.plot(smoothx,th_ihr,ls='--',c='y',lw=4)
    """
    pylab.show()


def family_skew():
    alpharange = np.arange(-15,-5,0.5) #5.01,0.5)
    cs = range(1,101,5)
    
    xs = np.arange(0.,2.,0.01)
    pylab.close('all')
    pylab.figure()
    for alpha in alpharange:
        ys = skew_func(xs, alpha)
        pylab.plot(xs,ys)
    pylab.savefig('sample_org.png')
    
    pylab.figure()
    for alpha in alpharange:
        ys = skew_func(xs, alpha,b=2.)
        pylab.plot(xs,ys)
    pylab.savefig('sample_b2.png')

    
    alpharange = [-5, -10, -15]
    cs = range(100,401,50)
    for c in cs:
        pylab.close('all')
        pylab.figure()
        for alpha in alpharange:
            ys = skew_func(smooth4x, b=600,a=alpha,c=c,d=100)
            pylab.plot(smooth4x,ys)
        pylab.savefig('sample_c%g_alpha%g.png'%(c,alpha))    
    
family_skew()
