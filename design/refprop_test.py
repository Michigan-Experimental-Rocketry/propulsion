import os, numpy as np
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

def NBP():
    RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    RP.SETPATHdll(os.environ['RPPREFIX'])
    print(RP.RPVersion())
    
    MOLAR_BASE_SI = RP.GETENUMdll(0,"MOLAR BASE SI").iEnum
    p_Pa = 101325
    Q = 0.0
    r = RP.REFPROPdll("Water","PQ","T",MOLAR_BASE_SI,0,0,p_Pa,Q,[1.0])
    print(r.Output[0])

    # r = RP.REFPROPdll("PROPANE","PQ","T",MOLAR_BASE_SI, 0,0,101325, 0, [1.0])
    # print(r.ierr, r.herr, r.Output[0])

if __name__=='__main__':
    # If the RPPREFIX environment variable is not already set by your installer (e.g., on windows), 
    # then uncomment this line and set the absolute path to the location of your install of 
    # REFPROP
    os.environ['RPPREFIX'] = r'/home/jasonyc/REFPROP-cmake/build'
    
    # Print the version of REFPROP in use and the NBP
    NBP()
