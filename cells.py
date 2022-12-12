import numpy as np

from params import *
###############CLASSES######################################
class Tconv():
    def __init__(self,t):
        self.r = np.random.lognormal(m0_Tc,sig0_Tc) # IL2R level
        self.i=0 #quantity of IL-2 consumed since its birth
        self.db=t #date of birth
      	
    
    def consume(self,ci):
    	self.i+=ci
    	
    def il2R_upregulation(self,nbR):
    	self.r+=nbR

        
class Treg():
    def __init__(self,t):
        self.r = np.random.lognormal(m0_Tr,sig0_Tr) # IL2R level
        self.i=0 #quantity of IL-2 consumed since its birth
        self.db=t #date of birth
    	
    def consume(self,ci):
    	self.i+=ci

    
    def il2R_upregulation(self,nbR):
    	self.r+=nbR
        
