# this is used to generate a phase map from the EOBD E field using data from COMSOL
# and will also attempt to do get the mode coupling factors.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pykat.commands import *
import pykat as pk
from scipy.interpolate import griddata

# import data and put them into arrays
data = pd.read_csv('EOBD_1_v.txt', delimiter='\s+', skiprows=9, header=None)
x = data[0]
y = data[1]
Ey = data[2]
Ex = data[3]

# center field about the origin
x = np.subtract(x, 0.0025)
y = np.subtract(y, 0.002)

# We dont necessarily care for the field that corresponds to the silicon
# so I remove it.

n_x = []
for i in range(int(len(x))):
    if (x[i] >= -.002 and x[i] <= 0.00201) == True:
        n_x.append(i)
        
x_c = []
y_c = []
Ey_c = []
Ex_c = []

for j in n_x:
    x_c.append(x[j])
    y_c.append(y[j])
    Ey_c.append(Ey[j])
    Ex_c.append(Ex[j])

x_c = np.array(x_c)
y_c = np.array(y_c)
Ey_c = np.array(Ey_c)
Ex_c = np.array(Ex_c)

# same for y
n_y = []
for i in range(int(len(x_c))):
    if (y_c[i] >= -.0020 and y_c[i] <= 0.0020) == True:
        n_y.append(i)
        
x_cc = []
y_cc = []
Ey_cc = []
Ex_cc = []

for j in n_y:
    x_cc.append(x_c[j])
    y_cc.append(y_c[j])
    Ey_cc.append(Ey_c[j])
    Ex_cc.append(Ex_c[j])

x_cc = np.array(x_cc)
y_cc = np.array(y_cc)
Ey_cc = np.array(Ey_cc)
Ex_cc = np.array(Ex_cc)


# convert E-field components to phase 
n0 = [1.774 for i in x_cc] # index of refraction for RTP in the y direction (y-cut)
r33 = 38.5e-12 # [pm/V] electro optic coefficient
n_E = [j+0.5*(j**3)*r33*i for i,j in zip(Ey_cc, n0)] # change of n in E_field
dOPL = [i*40e-3 - 7.096e-2 for i in n_E] # optical path length change
phase = [k*2*np.pi/1064e-9 for k in dOPL]

# now we will interpolate the data to make numerical integration easier

X = np.linspace(-0.002, 0.002, 101)
Y = X
XX, YY = np.meshgrid(X, Y)
phase_m = griddata((x_cc, y_cc), phase, (XX, YY), method='cubic', fill_value=0)

# will save grid to csv for importing to Matlab
# if using simtool you will need to make a separate grid in matlab
#df = pd.DataFrame(data=phase_m.astype(float))
#df.to_csv('Phase_map.csv', sep=' ', header=False, float_format='%.2f', index=False)


def HG_mode_content(w0, z, beam_data, x, y, n, m):
    
    import pykat.optics.gaussian_beams as gb
    q = gb.gauss_param(w0=w0, z=z)
    
    for i in range(0, n+1):
        for j in range(0, m+1):
            HG_mode = gb.HG_mode(q, n=i, m=j)
            HG_field = HG_mode.Unm(X,Y)
            k_nm = np.vdot(beam_data, HG_field)*np.diff(X)[0]*np.diff(Y)[0]
            print('%i%i: %.20F   %.20F' % (i, j, np.real(k_nm), np.imag(k_nm)))

 def HG_mode_content_1(w0, z, beam_data, x, y, n, m):
    
    import pykat.optics.gaussian_beams as gb
    q = gb.gauss_param(w0=w0, z=z)
    
    for i in range(0, n+1):
        for j in range(0, m+1):
            HG_mode = gb.HG_mode(q, n=i, m=j)
            HG_field = HG_mode.Unm(X,Y)
            mode_power = np.sum(np.multiply(np.conj(beam_data), HG_field))*np.diff(X)[0]*np.diff(Y)[0]
            print('%i%i: %.20F   %.20F' % (i, j, np.real(k_nm), np.imag(k_nm)))

HG_mode_content(250e-6, 0, EOBDXHG, X, Y, 3, 3)