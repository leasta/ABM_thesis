import sys
import numpy as np
from cells import *
import globals
from params import *

def death_division_starvation(mu_Tc,mu_Tr,lb_tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,dt,timenow):
	if mu_Tc>0 or lb_Tc>0 or starv_Tc>0:
		tc=len(globals.profile_Tc)
		# select the cells that will do something during dt
		ran_Tc = np.random.random(tc)
		prob_Tc=lb_Tc+mu_Tc+starv_Tc
		action_Tc=[T for i,T in enumerate(globals.profile_Tc) if ran_Tc[i]<prob_Tc*dt]
		
		if mu_Tc+starv_Tc==0 or lb_Tc==0:
			ran1_Tc=np.random.random(len(action_Tc))
			prob1_Tc=max(mu_Tc+starv_Tc,lb_Tc)
			do_Tc=[T for i,T in enumerate(action_Tc) if ran1_Tc[i]<prob1_Tc*dt] #wait starv only if il2<threshold...
			if mu_Tc+Starv_Tc==0:
				#cells divide
			else:
				#cells die
				dyingcells=[T for i,T in enumerate(action_Tr) if ran1_Tr[i]<mu_Tc]
			
		else:
			#Discrimitate cells that will die (naturally or by starvation) or divide
	if mu_Tr>0 or lb_Tr>0:
		tr=len(globals.profile_Tr)
		# select the cells that will do something during dt
		ran_Tr = np.random.random(tr)
		prob_Tr=lb_Tr+mu_Tr
		action_Tr=[T for i,T in enumerate(globals.profile_Tr) if ran_Tr[i]<prob_Tr*dt]
		
		#Discrimitate cells that will die or divide
		

def death_division_starvation2(mu_Tc,mu_Tr,lb_tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,dt,timenow):
	'''A certain number of cells dies, starves or divide during dt at rate mu_Tc for convention T cells, mu_tr for Tregs'''
	#There is a need for spliting all the cases... but I have 64 different cases.. help.
	tc=len(globals.profile_Tc)
	tr=len(globals.profile_Tr)
	# select the cells that will do something during dt
	ran_Tc = np.random.random(tc)
	ran_Tr = np.random.random(tr)
	prob_Tc=lb_Tc+mu_Tc+starv_Tc
	prob_Tr=lb_Tr+mu_Tr
	action_Tc=[T for i,T in enumerate(globals.profile_Tc) if ran_Tc[i]<prob_Tc*dt]
	action_Tr=[T for i,T in enumerate(globals.profile_Tr) if ran_Tr[i]<prob_Tr*dt]
	
	#select cells that die and make them die
	ran1_Tc=np.random.random(len(action_Tc))
	ratio_death_Tc=mu_Tc/prob_Tc
	dyingcells_Tc=[T for i,T in enumerate(action_Tc) if ran1_Tc[i]<ratio_death_Tc]
	others_Tc=[T for i,T in enumerate(action_Tc) if ran1_Tc[i]>=ratio_death_Tc]
	globals.profile_Tc=[T for T in globals.profile_Tc if  T not in dyingcells_Tc] #cells that died are erased from the list
	
	ran1_Tr=np.random.random(len(action_Tr))
	ratio_death_Tr=mu_Tr/prob_Tr
	dyingcells_Tr=[T for i,T in enumerate(action_Tr) if ran1_Tr[i]<ratio_death_Tr]
	mothers_Tr=[T for i,T in enumerate(action_Tr) if ran1_Tr[i]>=ratio_death_Tr and (T.i>thrdiv_Tr)]
	globals.profile_Tr=[T for T in globals.profile_Tr if  T not in dyingcells_Tr] # cells that died are erased from the list
	
	# Remaining cells divide or starve (Tconv only)
	ran2_Tc=np.random.random(len(others_Tc))
	ratio_div_Tc=lb_Tc/(lb_Tc+starv_Tc) #if starv =0 then all cells divide- if div=0 then...all cells starv ok. if div+starv=0 then error FIX THIS
	mothers_Tc=[T for i,T in enumerate(others_Tc) if ran2_Tc[i]<ratio_div_Tc]
	starv_Tc=[T for i,T in enumerate(others_Tc) if (ran2_Tc[i]>=ratio_div_Tc) and (T.i<thrstarv_Tc)]
	globals.profile_Tc=[T for T in globals.profile_Tc if  T not in starv_Tc] #cells that starved are erased from the list
	
	daughters_Tc=[Tconv(timenow) for i_ in range(len(mothers_Tc))]
	daughters_Tr=[Treg(timenow) for i_ in range(len(mothers_Tr))]
	for i in range(len(mothers_Tc)): #daughter cells receive half the receptors from their mother
		T=mothers_Tc[i]
		Td=daughters_Tc[i]	
		Td.r=T.r/2
		T.r=T.r/2
	for i in range(len(mothers_Tr)):  #daughter cells receive half the receptors from their mother
		T=mothers_Tr[i]
		Td=daughters_Tr[i]	
		Td.r=T.r/2
		T.r=T.r/2
	
############################ BIRTH #############################################
def immigration(a_Tc,a_Tr,dt,timenow):
	'''A certain number of cells are created during dt at rate im_Tc for conventional T cells, im_Tr for Tregs'''
	new_Tc = np.random.poisson(lam = a_Tc*dt)
	newTc=[Tconv(timenow) for _ in range(new_Tc)]
	globals.profile_Tc.extend(newTc)
	new_Tr = np.random.poisson(lam = a_Tr*dt)
	newTr=[Treg(timenow) for _ in range(new_Tr)]
	globals.profile_Tr.extend(newTr)
	

############### IL-2 CONSUMPTION  AND IL-2R upregulation ####################### 
def IL2consumption_IL2production_IL2Rupregulation(c_Tc,c_Tr,p_Tc,u_Tc,u_Tr,dt): #here, all T cells are consuming il-2 (no probability) 
	'''At each step,  activated T cells produce some IL-2. All this IL-2 will be consumed by Tregs and activated Tconv proportionally to their number of IL-2R
	During dt, Tregs and activated T conv upregulates their number of IL-2R depending on the amount of IL-2 they consumed since birth (Tregs and Tconvs) or activation (Tconvs only)'''
	production=p_Tc*len(globals.profile_Tc)*dt 
	sumRtc=np.sum(np.array([T.r for T in globals.profile_Tc]))
	sumRtr=np.sum(np.array([T.r for T in globals.profile_Tr]))
	norm=c_Tc*sumRtc+c_Tr*sumRtr
	if norm>0:
		[T.consume(c_Tr*T.r*production/norm) for T in globals.profile_Tr ] # all Tregs consume IL-2
		[T.consume(c_Tc*T.r*production/norm) for T in globals.profile_Tc] #only activated t conv consume IL-2
	else:
		print('normalisation constant equal to 0: no cells so no consumption')

	[T.il2R_upregulation(u_Tc*dt*T.i) for T in globals.profile_Tc] #only activated Tconv upregulate IL-2R expression
	[T.il2R_upregulation(u_Tr*dt*T.i) for T in globals.profile_Tr] # all Tregs upregulate I-2R expression
	return(norm)


################### Theory #################################################
def theoretical_NR(c_Tc,c_Tr,mu_Tc,mu_Tr,lb_tc,lb_Tr,thrdiv_Tr,starv_Tc,thrstarv_Tc,a_Tc,a_Tr,dt,timelist,Nc,Nr,R0,meanr0_Tc):
	'''Return the theoretical global functions R and N evaluated at each time of the simulation as an  numpy array '''
	timelist=np.array(timelist)
	#Conventional T cells only
	if Nr==0 and im_Tr==0:
		#fixed population
		if mu_Tc==0 and a_Tc==0 and lb_Tc==0 and starv_Tc==0:
			N=np.array([Nc]*len(timelist))
			R=u_Tc*c_Tc*p_Tc*Nc*timelist/2+R0
		#pure death process
		elif a_Tc==0 and lb_Tc==0 and starv_Tc==0 and mu_Tc>0:
			N=Nc*np.exp(-mu_Tc*timelist)
			R=(u_Tc*c_Tc*p_Tc*Nc*timelist/2+R0)*np.exp(-mu_tc*timelist)
		#death and immigration
		elif a_Tc>0 and my_Tc>0 and lb_Tc==0 and starv_Tc==0:
			N=(Nc-a_Tc/mu_Tc)*np.exp(-mu_Tc*timelist)+im_Tc/mu_Tc
			Rinf=Rinf=(a_Tc*c_Tc/mu_Tc)*(p_Tc*u_Tc/(mu_Tc**2)+meanr0_Tc)
			R=(c_Tc*u_Tc*p_Tc*(Nc-a_Tc/mu_Tc)/2*(timelist)**2-c_Tc*u_Tc*p_Tc*a_Tc/(mu_Tc**2)*np.array(timelist)-Rinf+R0)*np.exp(-mu_Tc*timelist)+Rinf
		#death, immigration, division
		elif a_Tc>0 and mu_Tc>0 and lb_Tc>0 and starv_Tc==0:
			Ninf=a_Tc/(mu_Tc-lb_Tc)
			N=(Nc-Ninf)*np.exp((lb_Tc-mu_Tc)*timelist)+Ninf
			Rinf=(a_Tc*c_Tc/mu_Tc)*(p_Tc*u_Tc/(mu_Tc**2-lb_Tc**2)+meanr0_Tc)
			R=c_Tc*u_Tc*p_Tc*(((np.exp(lb_Tc*timelist)-1)/(2*lb_Tc**2))*(Nc-Ninf)+((np.exp(-lb_Tc*timelist)-1)/(2*lb_Tc**2))*(Nc-(im_Tc)/(mu_Tc+lb_Tc))-Rinf+R0)*np.exp(-mu_Tc*timelist)+Rinf
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
