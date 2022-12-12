import numpy as np

from params import *
###############CLASSES######################################
class Tconv():
    def __init__(self,t):
        self.r = np.random.lognormal(m0_Tc,sig0_Tc) # IL2R level
        self.i=0 #quantity of IL-2 consumed since its birth
        self.db=t #date of birth
        
        ran=np.random.random()
        if ran<tracked_Tc:
        	self.track=True #trak 10% of cells
        else:
       		self.track=False
       	self.rlist=[self.r]
       	self.timelist=[t]
    
    def consume(self,ci):
    	self.i+=ci
    	
    def il2R_upregulation(self,nbR):
    	self.r+=nbR
    	
    def update_list(self,t):
    	if self.track==True:
    		self.rlist.append(self.r)
    		self.timelist.append(t)
    		
    def reset_tracking(self,t):
    	self.rlist=[self.r]
    	self.timelist=[t]

        
class Treg():
    def __init__(self,t):
        self.r = np.random.lognormal(m0_Tr,sig0_Tr) # IL2R level
        self.i=0 #quantity of IL-2 consumed since its birth
        self.db=t #date of birth
        
        ran=np.random.random()
        if ran<tracked_Tr:
        	self.track=True #trak 10% of cells
        else:
       		self.track=False
       	self.rlist=[self.r]
       	self.timelist=[t]
    	
    def consume(self,ci):
    	self.i+=ci

    
    def il2R_upregulation(self,nbR):
    	self.r+=nbR
    	
    def update_list(self,t):
    	if self.track==True:
    		self.rlist.append(self.r)
    		self.timelist.append(t)
    
    def reset_tracking(self,t):
    	self.rlist=[self.r]
    	self.timelist=[t]

        
