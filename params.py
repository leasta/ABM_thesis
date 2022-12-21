######################## PARAMETERS ############################################
import numpy as np

############### Tconvs ###############

a_Tc = 30			# immigration/activation rate
mu_Tc = 0.01  		# death rate
c_Tc = 1 	# comsumption rate of IL-2
p_Tc=10				# production rate of IL-2
lb_Tc =0.003 		# division rate (should be smaller than mu_Tc)
thrstarv_Tc =0	# threshold under which a cell can die of IL-2 starvation
starv_Tc=0		# death rate by starvation
u_Tc=10			# upregulation rate
m0_Tc=1			# parameter of the initial receptor distribution: mean
sig0_Tc=1			# parameter of the initial receptor distribution: std
tracked_Tc=0.003     # percentage of T convs tracked

############### Tregs ################
a_Tr = 0			# immigration/activation rate
mu_Tr = 0.01  		# death rate
c_Tr = 1 	# comsumption rate of IL-2
lb_Tr = 0.007 		# division rate (should be smaller than mu_Tr)
thrdiv_Tr = 500		# threshold of il-2 which make the cell eligible to division
u_Tr=5				# upregulation rate
m0_Tr=3				# parameter of the initial receptor distribution: mean
sig0_Tr=2			# parameter of the initial receptor distribution: std
tracked_Tr=0.003     # percentage of Tregs tracked
############### Time ###############
tmax=3000
dt=1
############### ANIMATION ##########

animation=True #Save a scatterplot for each time step. 
#animation=False # Faster simulation

###################### INITIALISATION ##########################################

Nc = int(a_Tc/(mu_Tc-lb_Tc))#initial number of Tconv
#Nc=10000

Nr=0 #initial number of Tregs


meanr0_Tc=np.exp(m0_Tc+sig0_Tc**2/2)
meanr0_Tr=np.exp(m0_Tr+sig0_Tr**2/2)
