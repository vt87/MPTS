import os
import copy
import glob
import math
import numpy as np
import random
import time
import sys
import BOO

class AnlysKMUtilsError(Exception):
    """ Raised if invalid thermodynamic ensemble is passed """

    def __init__(self, message):
        super().__init__(message)

'''
get the fraction of particles with solidlike order from steinhardt bond order parameters
datafile : datafile
trjfile : trajectory file
anlysparams : anlys params
'''
def getOPSTEINHARDT_AG_from_trj(datafile,trjfile,anlysparams):
    output = None    
    vardict = copy.deepcopy(anlysparams)

    #HARDCODED
    print(vardict["fortranPATH"])
    vardict["compilestr"] = "gfortran"

    #analyze 1
    vardict["selatoms"] = "type 1"
    opb = BOO.BOO(datafile,trjfile)
    oprunst = opb.BOOest(vardict,dict())
    clusdict = copy.deepcopy(opb.clusdict)
    output1 = clusdict["f_Totstar_avg"]
    if output1 is None:
        raise AnlysUtilsError("Q6 analysis has failed for %s"%trjfile)
    
    #analyze 2
    vardict["selatoms"] = "type 2"
    opb = BOO.BOO(datafile,trjfile)
    oprunst = opb.BOOest(vardict,dict())
    clusdict = copy.deepcopy(opb.clusdict)
    output2 = clusdict["f_Totstar_avg"]
    if output2 is None:
        raise AnlysUtilsError("Q6 analysis has failed for %s"%trjfile)

    #output
    if output1 is not None and output2 is not None:
        output = (output1 + output2)/2.0
        print(output1,output2,output)

    return output
