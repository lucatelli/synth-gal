import os
import numpy as np
import multiprocessing

from os import listdir
from os.path import isfile,join
from scipy import special
from scipy.special import gamma, gammainc
import astropy.io.fits as pf
import scipy
from scipy import signal

"""
Define some functions
"""
def Fn(n,Rn):
	"""
	Factor that converts luminosity to effectve brightness.
	"""
	return ((bn(n))**(2.0*n))/(2.0*np.pi*Rn**(2.0)*n*(np.exp(bn(n)))*(special.gamma(2.0*n)))

def bn(n):
	"""
	See Ciotti 1999.
	"""
	return 2.*n - 1./3.+(4./405.*n) + (46./25515.*n**2.0)

def w_values(file_name,p_value,path):
	'''
		Write on a file the input parameters given, for each galaxie.
	'''
	f = open(path+file_name+p_value[0]+'_values.dat', 'w')
	# with open(path+file_name+p_value[0]+'_values.glmk', 'w') as f:
	f.write('#gal_name,nb,Ib,Rb,Id,Rd,LT,BT,L_BT,L_DT,N, \n')
	for i in p_value:
		f.write(str(i))
		f.write(',')
	f.write('\n')
	f.close()

"""
Create table of values for bulge+disk profiles
"""
mock_params = ["X0","Y0", "FUNCTION","PA","ell","n","I_e","r_e","X0","Y0", "FUNCTION","PA","ell","n","I_e","r_e"]

number_of_images = 100

#BULGE+DISK RELATIONS
LT = np.random.uniform(1.0e2,1.0e4,size=(number_of_images,))
N  = np.random.randint(105,205,size=(number_of_images,))*2
# N[:] = 412
X0  = N/2.+0.5
Y0  = N/2.+0.5
PAB  = np.random.uniform(0,360,size=(number_of_images,))
ellB = np.random.uniform(0,0.4,size=(number_of_images,))

#BULGE_values
BT = np.random.uniform(0.00,0.9999,size=(number_of_images,))
# BT = np.zeros(number_of_images)
# BT[:] = 0.99999
L_BT=LT*BT
nb = np.random.uniform(0.9,6.0,size=(number_of_images,))
# nb = np.linspace(1.0,6.0,number_of_images)
# nb[:]=3.0
Rb = np.random.uniform(5.0,30.0,size=(number_of_images,))
Ib = Fn(nb,Rb)*L_BT 

#DISK_values
DT  =1.0-BT
L_DT=LT*DT
Rd = np.random.uniform(20.0,60.0,size=(number_of_images,))
nd = np.random.uniform(0.6,1.2,size=(number_of_images,))
# nd = 1.0
Id = Fn(nd,Rd)*L_DT

PAD = np.random.uniform(0,360,size=(number_of_images,))
ellD = np.random.uniform(0,0.5,size=(number_of_images,))

psf = "psf_new_less.fits"
path = "/run/media/sagauga/data/data/mock_imfit/"
psf_data = pf.getdata(psf)

for i in range(1,number_of_images):
	"""
	Now, each i is a mock image.
	"""
	print("Creating object ",i)
	mock_values = [X0[i],Y0[i],"Sersic",PAB[i],ellB[i],nb[i],Ib[i],Rb[i],X0[i],Y0[i],"Sersic",PAD[i],ellD[i],nd[i],Id[i],Rd[i]]
	#Now, to to each value of the model.
	f = open("make_image.dat", "w+")
	for k in range(len(mock_params)):
		f.write(str(mock_params[k])+" "+str(mock_values[k])+"\n")
	f.close()
	os.system("makeimage make_image.dat --nrows "+str(N[i])+" --ncols "+str(N[i])+" --psf "+psf+" -o "+path+"/gal_"+str(i)+".fits --max-threads 1")
	noise_level = 0.05
	noise = noise_level*np.random.rand(N[i],N[i])
	noise2 = 0.005*np.random.rand(N[i],N[i])
	noise_bkg = scipy.signal.fftconvolve(noise.copy(),psf_data,'same')
	g = pf.getdata(path+"/gal_"+str(i)+".fits")
	gn = g.copy() + noise_bkg+noise2
	pf.writeto(path+"/gal_"+str(i)+".fits",gn,overwrite=True)
	# os.system("makeimage make_image.dat --nrows "+str(N[i])+" --ncols "+str(N[i])+" -o "+path+"/gal_"+str(i)+".fits --max-threads 1")
	w_values("",("gal_"+str(i),nb[i],Ib[i],Rb[i],nd[i],Id[i],Rd[i],LT[i],BT[i],L_BT[i],L_DT[i],N[i]),path)

