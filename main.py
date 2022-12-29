import sys
import os
from datetime import date
import shutil
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from operator import add
import time
import imageio
import matplotlib.ticker
#import classes and functions
from cells import *
from functions import *
import globals
from params import *
import seaborn as sns
import pandas as pd 
import matplotlib.gridspec as gridspec


####################### SCATTER DATA ###########################################
def scatter_receptor_il2(mylist):
	'''Prepare the data for the scatterplot: receptor as a function of date of birth'''
	ncc=len(mylist)
	conv_ir=np.zeros((ncc,2))  
	for i in range(ncc):
		T=mylist[i]
		conv_ir[i,0]=T.timelist[0]
		conv_ir[i,1]=T.r
	return(conv_ir)
####################### UPDATE FUNCTION ########################################
def advance(c_Tc,c_Tr,mu_Tc,mu_Tr,lb_tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,a_Tc,a_Tr,dt,timenow,receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr):
	'''Update local and population variables'''
	print('t=',timenow)
	
	receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr=death_division_starvation(mu_Tc,mu_Tr,lb_tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,dt,timenow,receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr)
	immigration(a_Tc,a_Tr,dt,timenow)
	R=IL2consumption_IL2production_IL2Rupregulation(c_Tc,c_Tr,p_Tc,u_Tc,u_Tr,dt,timenow)
	return(R,receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr)
################### ANIMATE and PLOT ###########################################

def animateGraph(timenow):
	'''Update animated graph'''
	convTc_ir=scatter_receptor_il2([T for T in globals.profile_Tc])
	convTr_ir=scatter_receptor_il2([T for T in globals.profile_Tr])
	Trcells=len(globals.profile_Tr)
	Tccells=len(globals.profile_Tc)
	alivecells=globals.profile_Tc+globals.profile_Tr
	if len(alivecells)>0:
		indexold=np.argmin(np.array([T.timelist[0] for T in alivecells]))
		oldestcell=alivecells[indexold]
		oldage=oldestcell.timelist[0]
	else:
		oldage=0
	
	ax_xhist.cla()
	ax_yhist.cla()
	ax_xhist.set_ylim((0,ymax))
	ax_yhist.set_xlim((0,ymax))
	ax_yhist.set_xlabel('Number of cells')
	ax_xhist.set_ylabel('Number of cells')
	axrt=ax_main.set(xlabel='Date of birth $t_T^{in}$',ylabel='IL-2R expression level  $r_T$',xlim=(-10+oldage, timenow+100),ylim=(xmin1, xmax1),yscale='log')
	scattTc_convall.set_offsets(convTc_ir)
	
		

	ax_xhist.hist(convTc_ir[:,0],bins=50,histtype='step',label='Conventional T cells',color='darkturquoise')
	ax_yhist.hist(convTc_ir_init[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),alpha=0.1,color='blue',orientation='horizontal')
	ax_yhist.hist(convTc_ir[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),histtype='step',label='Conventional T cells',color='darkturquoise',orientation='horizontal')

	if Nr>0 or a_Tr>0:
		scattTr_convall.set_offsets(convTr_ir)
		ax_xhist.hist(convTr_ir[:,0],bins=50,histtype='step',label='Regulatory T cells',color='crimson')
		ax_yhist.hist(convTr_ir_init[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),alpha=0.1,color='red',orientation='horizontal')
		ax_yhist.hist(convTr_ir[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),histtype='step',label='Regulatory T cells',color='crimson',orientation='horizontal')
		ax_yhist.text(0,10**(17),r'$t=$'+str(timenow)+'\n'+r'$N_c=$'+str(Tccells)+r' $N_r=$'+str(Trcells), fontsize=25)

	else:
		ax_yhist.text(0,10**(17),r'$t=$'+str(timenow)+'\n'+r'$N_c=$'+str(Tccells), fontsize=25)

	
	
	
	
def animateAndSave(tmax,dt):
	#create a folder called Figure with the creation date.#counter for the day? 
	today=date.today()
	mydate = today.strftime("_%d%m%Y")
	i = 0
	while os.path.exists( '%s_Figure'%i+mydate+'/'):
		i += 1
	script_dir = os.path.dirname(__file__)
	results_dir = os.path.join(script_dir, '%s'%i+'_Figure'+mydate+'/')
	if not os.path.isdir(results_dir):
		os.makedirs(results_dir)
		#copy paste params.py in the new folder so save the parameter values we used
	original = os.path.join(script_dir, 'params.py')
	target = os.path.join(results_dir,'params.py')
	shutil.copyfile(original, target)
	
	#Initialise lists to track what we want
	timenow=0
	timelist=[0]
	sumRtc=np.sum(np.array([T.r for T in globals.profile_Tc]))
	sumRtr=np.sum(np.array([T.r for T in globals.profile_Tr]))
	Rlist=[c_Tc*sumRtc+c_Tr*sumRtr]
	Tregnumberlist=[len(globals.profile_Tr)]
	Tconvnumberlist=[len(globals.profile_Tc)]
	cellnumber=list( map(add, Tregnumberlist, Tconvnumberlist) )
	receptors_Tc=[]
	timelives_Tc=[]
	receptors_Tr=[]
	timelives_Tr=[]
	if animation==True: # Save picture of the inital conditions
		#images=[]
		file_name = "animation"+'%s'%timenow+'.png'
		#images.append(results_dir+file_name)	
		plt.savefig(results_dir + file_name)
	
	#Simulation
	while timenow<tmax:
		#update
		timenow+=dt
		R,receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr=advance(c_Tc,c_Tr,mu_Tc,mu_Tr,lb_Tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,a_Tc,a_Tr,dt,timenow,receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr)
		timelist.append(timenow)
		Rlist.append(R)
		
		#save figure at time t (timenow)
		if animation==True:
			animateGraph(timenow)
			file_name = "animation"+'%s'%timenow+'.png'
			plt.savefig(results_dir + file_name)
			
		Tregnumberlist.append(len(globals.profile_Tr))
		Tconvnumberlist.append(len(globals.profile_Tc))
		cellnumber=list( map(add, Tregnumberlist, Tconvnumberlist) )
		#update time
		if len(globals.profile_Tc+globals.profile_Tr)==0 and a_Tc==0 and a_Tr==0: #If no cells anymore then stop the simulation
			break
		
	#plot N, R and trajectories
	file_name='0_NR.png'
	R0=Rlist[0]
	Nth,Rth= theoretical_NR(c_Tc,c_Tr,mu_Tc,mu_Tr,lb_Tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,a_Tc,a_Tr,dt,timelist,Nc,Nr,R0,meanr0_Tc)
	fig,axs=plt.subplots(3,1,figsize=(15,10))
	params_plt = {'axes.labelsize': 24,
          'axes.titlesize': 30,
          'legend.fontsize': 24,
         'xtick.labelsize': 24,
         'ytick.labelsize': 24,
         'xtick.major.width': 2,
         'xtick.major.size': 10,
         'xtick.minor.width': 1,
         'xtick.minor.size': 5,
         'ytick.major.width': 2,
         'ytick.major.size': 10,
         'ytick.minor.width': 1,
         'ytick.minor.size': 5}
	plt.rcParams.update(params_plt)
	
	for i in range(2):
		axs[i].tick_params(width=2, length=7)
		for lab in axs[i].get_xticklabels():
			lab.set_fontsize(24)
		for lab in axs[i].get_yticklabels():
			lab.set_fontsize(24)
	# N
	#axs[0].set(ylabel='Number of cells',xlim=(0, tmax),ylim=(10**(-1),max( max(Tregnumberlist),max(Tconvnumberlist))+10000), yscale='log')
	axs[0].set(xlim=(0, tmax),ylim=(10**(-1),max(cellnumber)+10**3))#, yscale='log')
	axs[0].set_ylabel('Number of cells',fontsize=24)
	#axs[0].set_yscale('log')
	
	axs[0].plot(timelist,Tconvnumberlist,color='darkturquoise', lw=1,label='Conventional T cells')
	if Nth[-1]>0:
		axs[0].plot(timelist, Nth,linestyle='--',color='lightgray',label='Theory ')
		#axs[0].set_title(r'$N$')
	if Tregnumberlist[0]>0:
		axs[0].plot(timelist,Tregnumberlist,color='crimson', lw=1,label='Regulatory T cells')[0]
		axs[0].plot(timelist,cellnumber,color='black', lw=1,label='All T cells')
	axs[0].legend(loc='lower right')
	axs[0].set_xticks([])
	

	
	axs[1].set(ylabel=r'$R$', xlim=(0, tmax),ylim=(10**(1), max(Rlist)+10000),yscale='log')
	if Rth[-1]>0:
		axs[1].plot(timelist,Rth,linestyle='--',color='lime',label=r'Theoretical $\bar{R}$') 
	axs[1].plot(timelist,Rlist,color='black',lw=1,label=r'$R$')
	axs[1].legend(loc='best')
	axs[1].legend(loc='lower right')
	axs[1].set_xticks([])
	
	axs[2].set(xlim=(0, tmax),ylim=(10**(-1), 10**11), yscale='log')
	#axs[2].axhline(y=up_Tc*p_Tc/(d_Tc)**2,linestyle='--',color='pink',label=r'$\frac{up}{\mu^2}$') # do I want up/mu^2 or r0+up/mu^2?
	iTc=0
	iTr=0
	for T in globals.profile_Tc:
		if T.track==True:
			if iTc==0:
				axs[2].plot(T.timelist,T.rlist, color='blue',alpha=0.3,label='Tconvs')
			else:
				axs[2].plot(T.timelist,T.rlist, color='blue',alpha=0.3)
			iTc+=1
	if len(receptors_Tc)>0:
		for i in range(len(receptors_Tc)):
			if iTc==0:
				axs[2].plot(timelives_Tc[i],receptors_Tc[i],color='blue',alpha=0.3,label='Tconvs')
			else:
				axs[2].plot(timelives_Tc[i],receptors_Tc[i],color='blue',alpha=0.3)
			iTc+=1
		
	for T in globals.profile_Tr:
		if T.track==True:
			if iTr==0:
				axs[2].plot(T.timelist,T.rlist, color='red',alpha=0.3,label='Tregs')
			else:
				axs[2].plot(T.timelist,T.rlist, color='red',alpha=0.3)
			iTr+=1
	if len(receptors_Tr)>0:
		for i in range(len(receptors_Tr)):
			if iTr==0:
				axs[2].plot(timelives_Tr[i],receptors_Tr[i],color='red',alpha=0.3,label='Tregs')
			else:
				axs[2].plot(timelives_Tr[i],receptors_Tr[i],color='red',alpha=0.3)
			iTr+=1
	axs[2].legend(loc='upper right')
	axs[2].set_ylabel(r'$r_T$',fontsize=24)
	axs[2].set_xlabel(r'$t$',fontsize=24)
	
	
	print('saving figure...')
	plt.savefig(results_dir + file_name)
	print('Figure saved.')
	plt.close()
	return(results_dir)
	
################# RUN ##########################


#initialise variables	
globals.initialize(Nr,Nc)
print('Test initialization',len(globals.profile_Tc)==Nc,len(globals.profile_Tr)==Nr)



# Prepare axes for animation.
fig = plt.figure(figsize=(15,15))

params_plt = {'axes.labelsize': 24,
          'axes.titlesize': 30,
          'legend.fontsize': 24,
         'xtick.labelsize': 24,
         'ytick.labelsize': 24,
         'xtick.major.width': 2,
         'xtick.major.size': 10,
         'xtick.minor.width': 1,
         'xtick.minor.size': 5,
         'ytick.major.width': 2,
         'ytick.major.size': 10,
         'ytick.minor.width': 1,
         'ytick.minor.size': 5}
plt.rcParams.update(params_plt)

convTc_ir=scatter_receptor_il2([T for T in globals.profile_Tc])
convTr_ir=scatter_receptor_il2([T for T in globals.profile_Tr])
convTc_ir_init=convTc_ir
convTr_ir_init=convTr_ir

gs = gridspec.GridSpec(3, 3)
ax_main = plt.subplot(gs[1:3, :2])
ax_xhist = plt.subplot(gs[0, :2],sharex=ax_main)
ax_yhist = plt.subplot(gs[1:3, 2],sharey=ax_main)
xmin1=10**(-2)
xmax1=10**(11)
ymax=2000
axrt=ax_main.set(xlabel="Date of birth $t_T^{in}$", ylabel="IL-2R expression level $r_T$",xlim=(-10, 100),ylim=(xmin1,xmax1 ),yscale='log')
scattTc_convall=ax_main.scatter(convTc_ir[:,0],convTc_ir[:,1], color='darkturquoise',alpha=0.3, marker='o',s=5,label='Conventional T cells')

#ax_xhist.hist(conv_ir_init[:,0],bins=100,alpha=0.3,color='gray')
ax_xhist.hist(convTc_ir[:,0],bins=50,histtype='step',label='Conventional T cells',color='darkturquoise')

ax_xhist.set_ylabel('Number of cells')
ax_xhist.set_ylim((0,ymax))

ax_yhist.hist(convTc_ir_init[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),alpha=0.1,color='blue',orientation='horizontal')
ax_yhist.hist(convTc_ir[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),orientation='horizontal',histtype='step',label='Conventional T cells',color='darkturquoise')

ax_yhist.set_xlabel('Number of cells')
ax_yhist.set_xlim((0,ymax))

ax_yhist.yaxis.set_visible(False)
ax_xhist.xaxis.set_visible(False)
#ax_main.axhline(y=meanr0,label=r'$\bar{r_0}$',color='gray',linestyle='--')
#ax_yhist.axhline(y=meanr0,label=r'$\bar{r_0}$',color='gray',linestyle='--')
#ax_xhist.axvline(x=0,label=r'$t$',color='darkgray',linestyle='-',alpha=0.3)
if Nr>0 or a_Tr>0:
	ax_xhist.hist(convTr_ir[:,0],bins=50,histtype='step',label='Regulatory T cells',color='Crimson')
	ax_yhist.hist(convTr_ir[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),orientation='horizontal',histtype='step',label='Regulatory T cells',color='darkturquoise')
	ax_yhist.hist(convTr_ir_init[:,1],bins=10**np.linspace(np.log10(xmin1), np.log10(xmax1),50),alpha=0.1,color='red',orientation='horizontal')
	scattTr_convall=ax_main.scatter(convTr_ir[:,0],convTr_ir[:,1], color='crimson',alpha=0.3, marker='o',s=5,label='Regulatory T cells')
	ax_yhist.text(0,10**(17),r'$t=0$'+'\n'+r'$N_c=$'+str(Nc)+r' $N_r=$'+str(Nr), fontsize=25)
else:
	ax_yhist.text(0,10**(17),r'$t=0$'+'\n'+r'$N_c=$'+str(Nc), fontsize=25)

ax_main.legend(bbox_to_anchor=[1.65,1.3],markerscale=6)


results_dir=animateAndSave(tmax,dt)
print('All figures are created in folder '+results_dir)

