import numpy as np
import pylab

fileraw = 'doc/sarah_raw_data.txt'
filespikes = 'doc/sarah_spikes.txt'

"""
datraw = np.loadtxt(fileraw)
pylab.figure()
timeindex = len(datraw) -1
offset = 100
lenvector = len(datraw[timeindex,:])
for i in range(timeindex):
    baseline = np.ones(lenvector)*offset*i
    pylab.plot( datraw[timeindex,:],datraw[i,:]+baseline,lw=3)

pylab.savefig('raw_mua.png')

"""

lw=4
"""
lambdas = np.array([390.,410.,430.,445.,470.,500.,520.,540.,560.,585.,630.])
i_chr2 = np.array([700.,850.,920.,880.,880.,800.,600.,450.,200.,0.,0.])
i_nphr = np.array([175.,275.,300.,300.,400.,600.,800.,950.,1050.,1100.,450.])

i_archt = np.array([170,170,350,550,650,950,1000,800,350,50])
lambdas_archt = np.array([400,435,460,485,520,540,575,600,630,660])

pylab.figure()
pylab.plot(lambdas,i_chr2,lw=lw,c='#003399')
pylab.plot(lambdas,i_nphr,lw=lw,c='#FF9900')
pylab.plot(lambdas_archt,i_archt,lw=lw,c='#CCFF33')
pylab.savefig('activation_spectra.png')

"""
iallfile = 'experiments/130711_illumination_uneq_irradiance/dat/130711_illumination_uneq_irradiance_irr5.0_factor1.00_NpHR_whole_ChR_whole.dat'
ichrfile = 'experiments/130711_illumination_uneq_irradiance/dat/130711_illumination_uneq_irradiance_irr5.0_factor1.00_NpHR_whole_ChR_whole_iChR.dat'
ihrfile = 'experiments/130711_illumination_uneq_irradiance/dat/130711_illumination_uneq_irradiance_irr5.0_factor1.00_NpHR_whole_ChR_whole_iNpHR.dat'

datall = np.loadtxt(iallfile)
datchr = np.loadtxt(ichrfile)
dathr = np.loadtxt(ihrfile)

pylab.figure()
#soma
time = datall[0,:]
isoma = datall[4,:]
ichr = datchr[1:,1]
ihr = dathr[1:,1]

pylab.figure()
pylab.plot(time,ichr,lw=lw,c='#003399')
pylab.plot(time,ihr,lw=lw,c='#FF9900')
pylab.savefig('test_icurr.png')
