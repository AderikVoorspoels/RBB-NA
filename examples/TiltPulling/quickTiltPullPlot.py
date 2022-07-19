import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from math import *
a = 0.34
dx = 0.04
steps = 28
pos = 6
stiffnesses = []

fig, axs = plt.subplots(2,2, figsize=(20,20))
T=0

for k in range(0, steps+1, 1):
	i = dx*k
	PlumedRoll=np.loadtxt('DDD_TiltAt'+str(i)+'_Tilt.dat')
	axs[0,0].plot(PlumedRoll[::10,0]+T, PlumedRoll[::10,pos]*180/pi, label=str(i*180/pi))
	T+=PlumedRoll[-1,0]
	stiffnesses.append(1/np.var(PlumedRoll[:,pos]))
	hists, edges = np.histogram(PlumedRoll[:,pos], bins = 20)
	centers = (edges[1:]+edges[:-1])/2
	axs[0,1].plot(centers*180/pi, -1*np.log(hists))
	axs[1,0].plot(np.average(PlumedRoll[:,1:]*180/pi,axis =0))
	print(np.average(PlumedRoll[:,pos]))

	
	
	
axs[0,0].legend()
axs[0,0].set_xlabel('time [ps]')
axs[0,0].set_ylabel('tilt [deg]')


axs[1,1].plot(dx*np.array(range(0, steps+1, 1))*180/pi, stiffnesses)
axs[1,1].set_xlabel('at [deg]')
axs[1,1].set_ylabel('measured tilt stiffness [$k_{b}T$]')
axs[1,1].set_ylim([0,750])
fig.savefig('quickMovePlot.png')
