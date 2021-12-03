
# CAMB script, mostly copied from the readthedocs notebook.

import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

#Assume installed from github using "git clone --recursive https://github.com/cmbant/CAMB.git"
#This file is then in the docs folders
camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))
sys.path.insert(0,camb_path)
import camb
from camb import model, initialpower
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))

# Params for Small Box Runs
h = 0.675
Ob = 0.049 #0.0487
Om = 0.315 #0.31

# CDM
# Setup CAMB
pars = camb.CAMBparams()
pars.set_cosmology(H0=100*h, ombh2=Ob*h*h, omch2=(Om-Ob)*h*h)
pars.InitPower.set_params(As=2.0989e-9, ns=0.965, r=0)

# Get power spectra
pars.set_matter_power(redshifts=[0,199,200,201,997,1000,1003], kmax=10000)
pars.NonLinear = model.NonLinear_none
results = camb.get_transfer_functions(pars)
tf = results.get_matter_transfer_data().transfer_data
size = tf.shape[1]

# Files for CICASS:
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_cdm_transfer_z000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
            tf[0,ii,-1],tf[1,ii,-1],tf[2,ii,-1],tf[3,ii,-1],
            tf[4,ii,-1],tf[5,ii,-1],tf[6,ii,-1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_cdm_transfer_z1003.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
            tf[0,ii,0],tf[1,ii,0],tf[2,ii,0],tf[3,ii,0],
            tf[4,ii,0],tf[5,ii,0],tf[6,ii,0]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_cdm_transfer_z1000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
            tf[0,ii,1],tf[1,ii,1],tf[2,ii,1],tf[3,ii,1],
            tf[4,ii,1],tf[5,ii,1],tf[6,ii,1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_cdm_transfer_z997.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
            tf[0,ii,2],tf[1,ii,2],tf[2,ii,2],tf[3,ii,2],
            tf[4,ii,2],tf[5,ii,2],tf[6,ii,2]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_cdm_transfer_z201.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
            tf[0,ii,3],tf[1,ii,3],tf[2,ii,3],tf[3,ii,3],
            tf[4,ii,3],tf[5,ii,3],tf[6,ii,3]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_cdm_transfer_z200.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
            tf[0,ii,4],tf[1,ii,4],tf[2,ii,4],tf[3,ii,4],
            tf[4,ii,4],tf[5,ii,4],tf[6,ii,4]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_cdm_transfer_z199.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
            tf[0,ii,5],tf[1,ii,5],tf[2,ii,5],tf[3,ii,5],
            tf[4,ii,5],tf[5,ii,5],tf[6,ii,5]))
f.close()

# Now set up for WDM
tf_cdm = np.copy(tf)
m_wdm = 10.0 # keV
nu = 1.12
a_10kev = 0.049 * (1/m_wdm)**1.11 * ((Om-Ob)/0.25)**0.11 * (h/0.7)**1.22
wdm_corr = (1+(a_10kev*tf_cdm[0,:,:])**(2*nu))**(-5/nu) # sqrt of WDM correction to Pk; see Viel et al. 2005

tf[1:-1,:,:] = tf_cdm[1:-1,:,:]*wdm_corr

f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm10_transfer_z000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,-1],tf[1,ii,-1],tf[2,ii,-1],tf[3,ii,-1],
                                                          tf[4,ii,-1],tf[5,ii,-1],tf[6,ii,-1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm10_transfer_z1003.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,0],tf[1,ii,0],tf[2,ii,0],tf[3,ii,0],
                                                          tf[4,ii,0],tf[5,ii,0],tf[6,ii,0]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm10_transfer_z1000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,1],tf[1,ii,1],tf[2,ii,1],tf[3,ii,1],
                                                          tf[4,ii,1],tf[5,ii,1],tf[6,ii,1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm10_transfer_z997.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,2],tf[1,ii,2],tf[2,ii,2],tf[3,ii,2],
                                                          tf[4,ii,2],tf[5,ii,2],tf[6,ii,2]))
f.close()



m_wdm = 5.0 # keV
a_5kev = 0.049 * (1/m_wdm)**1.11 * ((Om-Ob)/0.25)**0.11 * (h/0.7)**1.22
wdm_corr = (1+(a_5kev*tf_cdm[0,:,:])**(2*nu))**(-5/nu) # sqrt of WDM correction to Pk; see Viel et al. 2005

tf[1:-1,:,:] = tf_cdm[1:-1,:,:]*wdm_corr

f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm5_transfer_z000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,-1],tf[1,ii,-1],tf[2,ii,-1],tf[3,ii,-1],
                                                          tf[4,ii,-1],tf[5,ii,-1],tf[6,ii,-1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm5_transfer_z1003.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,0],tf[1,ii,0],tf[2,ii,0],tf[3,ii,0],
                                                          tf[4,ii,0],tf[5,ii,0],tf[6,ii,0]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm5_transfer_z1000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,1],tf[1,ii,1],tf[2,ii,1],tf[3,ii,1],
                                                          tf[4,ii,1],tf[5,ii,1],tf[6,ii,1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm5_transfer_z997.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,2],tf[1,ii,2],tf[2,ii,2],tf[3,ii,2],
                                                          tf[4,ii,2],tf[5,ii,2],tf[6,ii,2]))
f.close()


m_wdm = 2.5 # keV
a_2p5kev = 0.049 * (1/m_wdm)**1.11 * ((Om-Ob)/0.25)**0.11 * (h/0.7)**1.22
wdm_corr = (1+(a_2p5kev*tf_cdm[0,:,:])**(2*nu))**(-5/nu) # sqrt of WDM correction to Pk; see Viel et al. 2005

tf[1:-1,:,:] = tf_cdm[1:-1,:,:]*wdm_corr

f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm30_transfer_z000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,-1],tf[1,ii,-1],tf[2,ii,-1],tf[3,ii,-1],
                                                          tf[4,ii,-1],tf[5,ii,-1],tf[6,ii,-1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm30_transfer_z1003.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,0],tf[1,ii,0],tf[2,ii,0],tf[3,ii,0],
                                                          tf[4,ii,0],tf[5,ii,0],tf[6,ii,0]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm30_transfer_z1000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,1],tf[1,ii,1],tf[2,ii,1],tf[3,ii,1],
                                                          tf[4,ii,1],tf[5,ii,1],tf[6,ii,1]))
f.close()
f = open('/Users/fdavies/CICASS/vbc_transfer/TFs/camb_wdm30_transfer_z997.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,2],tf[1,ii,2],tf[2,ii,2],tf[3,ii,2],
                                                          tf[4,ii,2],tf[5,ii,2],tf[6,ii,2]))
f.close()




# Random other TF output:
pars.set_matter_power(redshifts=[1000], kmax=10000)
pars.NonLinear = model.NonLinear_none
results = camb.get_transfer_functions(pars)
tf = results.get_matter_transfer_data().transfer_data
size = tf.shape[1]
f = open('./camb_cdm_transfer_z1000.dat','w')
f.write('{:d}\n'.format(size))
for ii in range(size):
    f.write('{:E} {:E} {:E} {:E} {:E} {:E} {:E}\n'.format(
                                                          tf[0,ii,0],tf[1,ii,0],tf[2,ii,0],tf[3,ii,0],
                                                          tf[4,ii,0],tf[5,ii,0],tf[6,ii,0]))
f.close()
