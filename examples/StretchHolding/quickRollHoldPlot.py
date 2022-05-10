import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from math import *
a = 0.34
dk = 100
steps = 10
pos = 6
stiffnesses = []

fig, axs = plt.subplots(2,2, figsize=(20,20))

for k in range(0, steps+1, 1):
	i = dk*k
	PlumedRoll=np.loadtxt('DDD_TiltHold'+str(i)+'_Tilt.dat')
	axs[0,0].plot(PlumedRoll[::10,0], PlumedRoll[::10,pos]*180/pi, label=str(i))
	stiffnesses.append(1/np.var(PlumedRoll[:,pos]))
	hists, edges = np.histogram(PlumedRoll[:,pos], bins = 20)
	centers = (edges[1:]+edges[:-1])/2
	axs[0,1].plot(centers, -1*np.log(hists))
	axs[1,0].plot(np.average(PlumedRoll[:,1:],axis =0)*180/pi)
	print(np.average(PlumedRoll[:,pos]))

	
	
	
axs[0,0].legend()
axs[0,0].set_xlabel('time [ps]')
axs[0,0].set_ylabel('tilt [deg]')


axs[1,1].plot(dk*np.array(range(0, steps+1, 1)), stiffnesses)
axs[1,1].set_xlabel('kappa $kj/mol/(rad)^{-2}$')
axs[1,1].set_ylabel('gemeten tilt stiffness $k_{b}T$')
fig.savefig('quickHoldPlot.png')
