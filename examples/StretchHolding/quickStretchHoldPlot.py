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
	PlumedStretch=np.loadtxt('DDD_StretchHold'+str(i)+'_Stretch.dat')
	axs[0,0].plot(PlumedStretch[::10,0], PlumedStretch[::10,pos], label=str(i))
	stiffnesses.append(1/np.var(PlumedStretch[:,pos]))
	hists, edges = np.histogram(PlumedStretch[:,pos], bins = 20)
	centers = (edges[1:]+edges[:-1])/2
	axs[0,1].plot(centers, -1*np.log(hists))
	axs[1,0].plot(np.average(PlumedStretch[:,1:],axis =0))
	print(np.average(PlumedStretch[:,pos]))

	
	
	
axs[0,0].legend()
axs[0,0].set_xlabel('time [ps]')
axs[0,0].set_ylabel('Stretch [nm]')


axs[1,1].plot(dk*np.array(range(0, steps+1, 1)), stiffnesses)
axs[1,1].set_xlabel('kappa [$kj/mol/(nm)^{-2}$]')
axs[1,1].set_ylabel('measured Stretch stiffness [$k_{b}T/(nm)^{-2}$]')
fig.savefig('quickHoldPlot.png')
