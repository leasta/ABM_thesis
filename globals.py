from cells import *
from params import *
def initialize(Nr,Nc): 
    global profile_Tc
    global profile_Tr
    profile_Tc = [Tconv(0) for i in range(Nc)]
    profile_Tr = [Treg(0) for i in range(Nr)]
