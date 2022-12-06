######################## PARAMETERS ############################################
import numpy as np
############### Tconvs ###############
a_Tc = 30 			# immigration/activation rate
mu_Tc = 0.01  		# death rate
c_Tc = 1 	# comsumption rate of IL-2
p_Tc=10				# production rate of IL-2
lb_Tc = 0 		# division rate
thrstarv_Tc = 0		# threshold under which a cell can die of IL-2 starvation
starv_Tc=0			# death rate by starvation
u_Tc=10			# upregulation rate
m0_Tc=3				# parameter of the initial receptor distribution: mean
sig0_Tc=1			# parameter of the initial receptor distribution: std

############### Tregs ################
a_Tr = 0			# immigration/activation rate
mu_Tr = 0.01  		# death rate
c_Tr = 1 	# comsumption rate of IL-2
lb_Tr = 0.007 		# division rate
thrdiv_Tr = 200		# threshold of il-2 which make the cell eligible to division
u_Tr=5				# upregulation rate
m0_Tr=3				# parameter of the initial receptor distribution: mean
sig0_Tr=2			# parameter of the initial receptor distribution: std

############### Time ###############
tmax=1000
dt=10
###################### INITIALISATION ##########################################

#nc = int(im_Tc/(d_Tc-div_Tc))#initial number of Tconv
Nc=5000

Nr=0
#nr = 5000 #initial number of Tregs


meanr0_Tc=np.exp(m0_Tc+sig0_Tc**2/2)
meanr0_Tr=np.exp(m0_Tr+sig0_Tr**2/2)
