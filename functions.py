import sys
import numpy as np
from cells import *
import globals
from params import *

def dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc):
	''' Return death rate of Tconv T
	'''
	rr=mu_Tc
	if T.i<thrstarv_Tc:
		rr+=starv_Tc
	return(rr)
def divr_Tr(T,lb_Tr,thrdiv_Tr):
	''' Return the division rate of Treg T
	'''
	rr=0
	if T.i>thrdiv_Tr:
		rr=lb_Tr
	return(rr)
	
def death_division_starvation(mu_Tc,mu_Tr,lb_tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,dt,timenow,receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr):
	''' Choose cells that will divide or die during dt and update the profiles
	'''
	if mu_Tc+starv_Tc+lb_Tc==0:
		print('No death or division for Tconvs')
	else:
		# Determine cells that will do something
		tc=len(globals.profile_Tc)
		ran_Tc = np.random.random(tc)
		action_Tc=[T for i,T in enumerate(globals.profile_Tc) if ran_Tc[i]<(dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc)+lb_Tc)*dt]
		# Choose between death or division
		ran2_Tc = np.random.random(len(action_Tc))
		mothers_Tc=[T for i,T in enumerate(action_Tc) if np.searchsorted(np.cumsum([lb_Tc/(dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc)+lb_Tc),dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc)/(dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc)+lb_Tc)]),ran2_Tc[i])==0]
		dyingcells_Tc=[T for i,T in enumerate(action_Tc) if np.searchsorted(np.cumsum([lb_Tc/(dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc)+lb_Tc),dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc)/(dr_Tc(T,starv_Tc,mu_Tc,thrstarv_Tc)+lb_Tc)]),ran2_Tc[i])==1]
		# Update tracked lists
		for T in dyingcells_Tc:
			if T.track==True:
				receptors_Tc.append(T.rlist)
				timelives_Tc.append(T.timelist)
				
		# Update T conv profile
		globals.profile_Tc=[T for T in globals.profile_Tc if  T not in dyingcells_Tc]
		daughters_Tc=[Tconv(timenow) for i_ in range(len(mothers_Tc))]
		for i in range(len(mothers_Tc)): #daughter cells receive half the receptors from their mother
			T=mothers_Tc[i]
			Td=daughters_Tc[i]	
			Td.r=T.r/2
			T.r=T.r/2
			T.reset_tracking(timenow)
       			
		globals.profile_Tc.extend(daughters_Tc)
	if mu_Tr+lb_Tr==0:
		print('No death or division for Tregs')
	else:
		# Determine cells that will do something
		tr=len(globals.profile_Tr)
		ran_Tr = np.random.random(tr)
		action_Tr=[T for i,T in enumerate(globals.profile_Tr) if ran_Tr[i]<(divr_Tr(T,lb_Tr,thrdiv_Tr)+mu_Tr)*dt]
		# Choose between division or death
		ran2_Tr = np.random.random(len(action_Tr))
		mothers_Tr=[T for i,T in enumerate(action_Tr) if np.searchsorted(np.cumsum([divr_Tr(T,lb_Tr,thrdiv_Tr)/(divr_Tr(T,lb_Tr,thrdiv_Tr)+mu_Tr),mu_Tr/(mu_Tr+divr_Tr(T,lb_Tr,thrdiv_Tr))]),ran2_Tr[i])==0]
		dyingcells_Tr=[T for i,T in enumerate(action_Tr) if np.searchsorted(np.cumsum([divr_Tr(T,lb_Tr,thrdiv_Tr)/(divr_Tr(T,lb_Tr,thrdiv_Tr)+mu_Tr),mu_Tr/(mu_Tr+divr_Tr(T,lb_Tr,thrdiv_Tr))]),ran2_Tr[i])==1]
		# Update tracked lists
		for T in dyingcells_Tr:
			if T.track==True:
				receptors_Tr.append(T.rlist)
				timelives_Tr.append(T.timelist)
		# Update Tregs profile
		globals.profile_Tr=[T for T in globals.profile_Tr if T not in dyingcells_Tr]
		daughters_Tr=[Treg(timenow) for i_ in range(len(mothers_Tr))]
		for i in range(len(mothers_Tr)): #daughter cells receive half the receptors from their mother
			T=mothers_Tr[i]
			Td=daughters_Tr[i]	
			Td.r=T.r/2
			T.r=T.r/2
			T.reset_tracking(timenow)
		globals.profile_Tr.extend(daughters_Tr)
		return(receptors_Tc,timelives_Tc,receptors_Tr,timelives_Tr)
############################ BIRTH #############################################
def immigration(a_Tc,a_Tr,dt,timenow):
	'''A certain number of cells are created during dt at rate a_Tc for Tconvs, a_Tr for Tregs'''
	new_Tc = np.random.poisson(lam = a_Tc*dt)
	newTc=[Tconv(timenow) for _ in range(new_Tc)]
	globals.profile_Tc.extend(newTc)
	new_Tr = np.random.poisson(lam = a_Tr*dt)
	newTr=[Treg(timenow) for _ in range(new_Tr)]
	globals.profile_Tr.extend(newTr)
	

############### IL-2 CONSUMPTION  AND IL-2R upregulation ####################### 
def IL2consumption_IL2production_IL2Rupregulation(c_Tc,c_Tr,p_Tc,u_Tc,u_Tr,dt,timenow): #here, all T cells are consuming il-2 (no probability) 
	'''At each step, Tconvs produce some IL-2. All this IL-2 will be consumed by Tregs and Tconvs proportionally to their number of IL-2R
	During dt, Tregs and Tconvs upregulate their number of IL-2R depending on the amount of IL-2 they consumed since birth (Tregs and Tconvs) or activation (Tconvs only)'''
	production=p_Tc*len(globals.profile_Tc)*dt 
	sumRtc=np.sum(np.array([T.r for T in globals.profile_Tc]))
	sumRtr=np.sum(np.array([T.r for T in globals.profile_Tr]))
	R=c_Tc*sumRtc+c_Tr*sumRtr
	if R>0:
		[T.consume(c_Tr*T.r*production/R) for T in globals.profile_Tr ] # all Tregs consume IL-2
		[T.consume(c_Tc*T.r*production/R) for T in globals.profile_Tc] #only activated t conv consume IL-2
	else:
		print('normalisation constant equal to 0: no cells so no consumption')

	[T.il2R_upregulation(u_Tc*dt*T.i) for T in globals.profile_Tc] #only activated Tconv upregulate IL-2R expression
	[T.il2R_upregulation(u_Tr*dt*T.i) for T in globals.profile_Tr] # all Tregs upregulate I-2R expression
	
	[T.update_list(timenow) for T in globals.profile_Tc]
	[T.update_list(timenow) for T in globals.profile_Tr]
	return(R)


################### Theory #################################################
def theoretical_NR(c_Tc,c_Tr,mu_Tc,mu_Tr,lb_tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,a_Tc,a_Tr,dt,timelist,Nc,Nr,R0,meanr0_Tc):
	'''Return the theoretical global functions R and N evaluated at each time of the simulation as an Numpy array '''
	timelist=np.array(timelist)
	#Conventional T cells only
	if Nr==0 and a_Tr==0:
		#fixed population
		if mu_Tc==0 and a_Tc==0 and lb_Tc==0 and starv_Tc==0:
			print('Deterministic model')
			N=np.array([Nc]*len(timelist))
			R=u_Tc*c_Tc*p_Tc*Nc*timelist**2/2+R0
		#pure death process
		elif a_Tc==0 and lb_Tc==0 and starv_Tc==0 and mu_Tc>0:
			print('Hybrid model with death only')
			N=Nc*np.exp(-mu_Tc*timelist)
			R=(u_Tc*c_Tc*p_Tc*Nc*timelist**2/2+R0)*np.exp(-mu_Tc*timelist)
		#death and immigration
		elif a_Tc>0 and mu_Tc>0 and lb_Tc==0 and starv_Tc==0:
			print('Hybrid model with death and immigration')
			N=(Nc-a_Tc/mu_Tc)*np.exp(-mu_Tc*timelist)+a_Tc/mu_Tc
			Rinf=Rinf=(a_Tc*c_Tc/mu_Tc)*(p_Tc*u_Tc/(mu_Tc**2)+meanr0_Tc)
			R=(c_Tc*u_Tc*p_Tc*(Nc-a_Tc/mu_Tc)/2*(timelist)**2-c_Tc*u_Tc*p_Tc*a_Tc/(mu_Tc**2)*np.array(timelist)-Rinf+R0)*np.exp(-mu_Tc*timelist)+Rinf
		#death, immigration, division
		elif a_Tc>0 and mu_Tc>0 and lb_Tc>0 and starv_Tc==0:
			print('Hybrid model with death, immigration and division')
			Ninf=a_Tc/(mu_Tc-lb_Tc)
			N=(Nc-Ninf)*np.exp((lb_Tc-mu_Tc)*timelist)+Ninf
			Rinf=(a_Tc*c_Tc/mu_Tc)*(p_Tc*u_Tc/(mu_Tc**2-lb_Tc**2)+meanr0_Tc)
			R=c_Tc*u_Tc*p_Tc*(((np.exp(lb_Tc*timelist)-1)/(2*lb_Tc**2))*(Nc-Ninf)+((np.exp(-lb_Tc*timelist)-1)/(2*lb_Tc**2))*(Nc-(a_Tc)/(mu_Tc+lb_Tc)))*np.exp(-mu_Tc*timelist)+(R0-Rinf)*np.exp(-mu_Tc*timelist)+Rinf
		#others
		else:
			N=np.array([-1]*len(timelist))
			R=np.array([-1]*len(timelist))
			print('We did not computed the global variables')
	else:
		R=np.array([-1]*len(timelist))
		N=np.array([-1]*len(timelist))
		print('We did not computed the global values')
	return(N,R)
