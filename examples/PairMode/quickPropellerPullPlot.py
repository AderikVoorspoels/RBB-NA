import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from math import *
a = 0.34
dx = 0.01
steps = 10
pos = 6
stiffnesses = []

fig, axs = plt.subplots(2,2, figsize=(20,20))

for k in range(0, steps+1, 1):
	i = dx*k
	PlumedPropeller=np.loadtxt('DDD_PropellerAt'+str(i)+'_Propeller.dat')
	axs[0,0].plot(PlumedPropeller[::10,0], PlumedPropeller[::10,pos]*180/pi, label=str(i))
	stiffnesses.append(1/np.var(PlumedPropeller[:,pos]))
	hists, edges = np.histogram(PlumedPropeller[:,pos], bins = 20)
	centers = (edges[1:]+edges[:-1])/2
	axs[0,1].plot(centers, -1*np.log(hists))
	axs[1,0].plot(np.average(PlumedPropeller[:,1:],axis =0)*180/pi)
	print(np.average(PlumedPropeller[:,pos]))

	
	
	
axs[0,0].legend()
axs[0,0].set_xlabel('time [ps]')
axs[0,0].set_ylabel('propeller [deg]')


axs[1,1].plot(dx*np.array(range(0, steps+1, 1))*180/pi, stiffnesses)
axs[1,1].set_xlabel('at [deg]')
axs[1,1].set_ylabel('measured propeller stiffness [$k_{b}T$]')
fig.savefig('quickPullPlot.png')
