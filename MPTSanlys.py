#!/usr/bin/env python
#-----------------------------------------------------------------------------------------------
# Vikram Thapar, vt87@cornell.edu, Prof. Escobedo's lab, Cornell University
# Analysis of the framework, Model Parameter Targeted Search (MPTS)
# File Name: MPTSanlys.py
#-----------------------------------------------------------------------------------------------

import os
import sys
import argparse
import time
import yaml
import itertools
import csv
import pandas as pd
import copy
import shutil
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Global plot setting
plt.rc('font', family='serif')
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


class MPTSAnlysError(Exception):
    """ Raised the error in MPTSanlys.py """

    def __init__(self, message):
        super().__init__(message)

def initialize():
    parser = argparse.ArgumentParser(description='Analyzing MPTS')
    parser.add_argument('-y', '--yaml', type=str, help='A .yaml configuration file')
    return parser

'''
Collecting proposal files
mpts_params : parameters for mpts
'''
def collectprpsl(anlys_params):
    
    #list of dictionaries and a dictionary
    outlst = []
    alldict = dict()
    isSkip = False
    stpoint = 0
    anlys_params["skipstr"] = ""

    #copy args
    rundir = os.path.abspath(mpts_params["out_path"])
    
    #change directory
    os.chdir(rundir)

    #collect prpslfiles
    files = glob.glob("prpsl_*")
    if len(files) == 0:
        print("No proposal files found. Returning....")
        return outlst

    #read the files
    for i in range(stpoint,len(files)):
        tmpdict = dict()
        #filename
        f = "prpsl_%s.txt"%i
        tmpdict["filename"] = f
        tmpdict["id"] = i
        #pandas
        df = pd.read_csv(f,delimiter = " ", dtype = str, na_filter = False)
        tmpdict["out"] = df.to_dict('list')
        if i == stpoint:
            alldict = copy.deepcopy(tmpdict["out"])
        else:
            for key,value in tmpdict["out"].items():
                for j in range(len(value)):
                    alldict[key].append(value[j])
        outlst.append(copy.deepcopy(tmpdict))        

    print("********************************************************************")
    
    return outlst,alldict

'''
Plotting data in proposal files
mpts_params : parameters for the mpts
inlst : List of dictionaries containing data
'''
def pltprpsl(mpts_params,inlst):
    
    #copy args
    rundir = os.path.abspath(mpts_params["out_path"])


    #output dir    
    outdir = "%s/Anlys"%(rundir)
    os.makedirs(outdir,exist_ok = True)

    #change directory
    os.chdir(outdir)

    #length check
    if len(inlst) == 0:
        print("No proposal data found. Returning....")
        return

    #read the files
    X = []  #iteration number
    Y = []  # data
    hyp = []
    for i in range(len(inlst)):
        d = inlst[i]
        Output = d["out"]["Output"]
        index = d["id"] 
        for j in range(len(Output)):
            out = Output[j]
            if out != "NA" and out != "FAIL":
                X.append(index)
                Y.append(out)
    #array
    X = np.asarray(X).astype(int)
    Y = np.asarray(Y).astype(float)

    #save file
    outstr = "Prpsl"
    outfile = "%s.txt"%(outstr)
    outimg = "%s.png"%(outstr)
    pltdict = dict()
    pltdict["iter"] = X
    pltdict["Output"] = Y


    print("Saving the iteration vs output data in %s/%s"%(outdir,outfile))
    print("Plotting the iteration vs output data in %s/%s"%(outdir,outfile))

    #write 
    df = pd.DataFrame.from_dict(pltdict)
    df.to_csv(outfile,index = False, sep = " ")

    #write image
    fig = plt.figure()
    ax1 = fig.add_subplot(111)    
    ax1.set_xlabel("ML Iteration")
    ax1.set_ylabel("Output")
    ax1.plot(X,Y,'ro')
    plt.savefig(outimg)
    plt.close()
        
    print("********************************************************************")


#---------------------------------------------------------------------------------
# MAIN FUNCTION
#---------------------------------------------------------------------------------     
def main(mpts_params):

    workpath = os.getcwd()

    #collect and plot data in prpsfiles
    outlst,alldict = collectprpsl(mpts_params)
    os.chdir(workpath)
    pltprpsl(mpts_params,outlst)
    os.chdir(workpath)


# Execute main function if this file is the starting point
if __name__ == '__main__':
    start_time = time.time()

    #Read arguments
    args = initialize().parse_args()    

    #making 
    if args.yaml:
        with open(args.yaml, 'r') as yml:
            cfg = yaml.load(yml, Loader=yaml.FullLoader)
            mpts_params = cfg['mpts']
    else:
        mpts_params = dict()
    main(mpts_params)

    print("--- %s seconds ---" % (time.time() - start_time))
