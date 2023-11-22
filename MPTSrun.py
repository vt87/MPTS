#!/usr/bin/env python
#-----------------------------------------------------------------------------------------------
# Vikram Thapar, vt87@cornell.edu, Prof. Escobedo's lab, Cornell University
# Execution of the framework, Model Parameter Targeted Search (MPTS)
# File Name: MPTS_main.py
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
import importlib
import pandas as pd
import numpy as np

class MPTSRunError(Exception):
    """ Raised the error in MPTSmain.py """

    def __init__(self, message):
        super().__init__(message)

def initialize():
    parser = argparse.ArgumentParser(description='Running MPTS')
    parser.add_argument('-y', '--yaml', type=str, help='A .yaml configuration file')
    return parser

'''
default values of arguments for this main script in case they are not defined in yaml
mpts_params : Input dictionary for this mpts script
'''
def defaults(mpts_params):
    #print("********************************************************************")
    #print('Setting default values of arguments in mpts_params.')
    if "workpath" not in mpts_params.keys():
        mpts_params["workpath"] = os.getcwd()
    if "sleeptime" not in mpts_params.keys():
        mpts_params["sleeptime"] = 60
    #print("Done....")
    #print("********************************************************************")
    return mpts_params


'''
error check
mpts_params : Input dictionary for this mptsrun script
'''
def errorcheck(mpts_params):
    #print("********************************************************************")
    print('Checking for errors in the initial set up.')

    #inp_path
    if "inp_path" not in mpts_params.keys():
        raise MPTSRunError("Argument named inp_path not found in MPTS input file.")
    else:
        inp_path = mpts_params["inp_path"]
        if not os.path.exists(inp_path):
            raise MPTSRunError("Input path named %s does not exist."%inp_path)
        mpts_params["inp_path"] = os.path.abspath(mpts_params["inp_path"])

    #out_path
    if "out_path" not in mpts_params.keys():
        raise MPTSRunError("Argument named out_path not found in MPTS input file.")
    else:
        out_path = mpts_params["out_path"]
        if not os.path.exists(out_path):
            raise MPTSRunError("Output path named %s does not exist."%out_path)
        mpts_params["out_path"] = os.path.abspath(mpts_params["out_path"])

    #module_path
    if "module_path" not in mpts_params.keys():
        raise MPTSRunError("Argument named module_path not found in MPTS input file.")
    else:
        module_path = mpts_params["module_path"]
        if not os.path.exists(module_path):
            raise MPTSRunError("Module_path named %s does not exist."%module_path)
        mpts_params["module_path"] = os.path.abspath(mpts_params["module_path"])

    #module_inppath
    if "module_inppath" not in mpts_params.keys():
        raise MPTSRunError("Argument named module_inppath not found in MPTS input file.")
    else:
        module_inppath = mpts_params["module_inppath"]
        if not os.path.exists(module_inppath):
            raise MPTSRunError("Module_inppath named %s does not exist."%module_inppath)
        mpts_params["module_inppath"] = os.path.abspath(mpts_params["module_inppath"])

    #niter
    if "niter" not in mpts_params.keys():
        raise MPTSRunError("Argument named niter (Number of iterations) not found in MPTS input file.")

    #mliter
    if "mliter" not in mpts_params.keys():
        raise MPTSRunError("Argument named mliter (Number of machine learning iterations) not found in MPTS input file.")
            
    #nextsize
    if "nextsize" not in mpts_params.keys():
        raise MPTSRunError("Argument named nextsize (Batch size for the next best batch) not found in MPTS input file.")

    #libname
    if "libname" not in mpts_params.keys():
        raise MPTSRunError("Argument named libname (Library name) not found in MPTS input file.")
    else:
        libname = mpts_params["libname"]
        libdatafile = "%s/%s"%(mpts_params["inp_path"],mpts_params["libname"])
        if not os.path.exists(libdatafile):
            raise MPTSRunError("Library file named %s does not exist"%libdatafile)
        mpts_params["libdatafile"] = libdatafile

    #initfile
    if "initfile" not in mpts_params.keys():
        raise MPTSRunError("Argument named initfile (initial batch library file) not found in MPTS input file.")
    else:
        initfile = mpts_params["initfile"]
        libinitdatafile = "%s/%s"%(mpts_params["inp_path"],mpts_params["initfile"])
        if not os.path.exists(libinitdatafile):
            raise MPTSRunError("Initial batch library file named %s does not exist"%libinitdatafile)
        mpts_params["libinitdatafile"] = libinitdatafile
        
        
    #sim modules
    if "module_sim" not in mpts_params.keys():
        raise MPTSRunError("Argument of dictionary type named module_sim not found in MPTS input file.")
    module_sim = mpts_params["module_sim"]
    if "name" not in module_sim.keys():
        raise MPTSRunError("Argument named name in module_sim dictionary not found in MPTS input file.")
    else:
        pyfile = "%s/%s.py"%(mpts_params["module_path"],module_sim["name"])
        if not os.path.exists(pyfile):
            raise MPTSRunError("Simulation module py file named %s does not exist"%pyfile)
        mpts_params["module_sim"]["pyfile"] = pyfile
    if "inpname" not in module_sim.keys():
        raise MPTSRunError("Argument named inpname in module_sim dictionary not found in MPTS input file.")
    else:
        inpfile = "%s/%s"%(mpts_params["module_inppath"],module_sim["inpname"])
        if not os.path.exists(inpfile):
            raise MPTSRunError("Simulation module input file named %s does not exist"%inpfile)
        mpts_params["module_sim"]["inpfile"] = inpfile
    if "vars" not in module_sim.keys():
        mpts_params["module_sim"]["vars"] = dict()


    #anlys modules
    if "module_anlys" not in mpts_params.keys():
        raise MPTSRunError("Argument of dictionary type named module_anlys not found in MPTS input file.")
    module_anlys = mpts_params["module_anlys"]
    if "name" not in module_anlys.keys():
        raise MPTSRunError("Argument named name in module_anlys dictionary not found in MPTS input file.")
    else:
        pyfile = "%s/%s.py"%(mpts_params["module_path"],module_anlys["name"])
        if not os.path.exists(pyfile):
            raise MPTSRunError("Analysis module py file named %s does not exist"%pyfile)
        mpts_params["module_anlys"]["pyfile"] = pyfile
    if "inpname" not in module_anlys.keys():
        raise MPTSRunError("Argument named inpname in module_anlys dictionary not found in MPTS input file.")
    else:
        inpfile = "%s/%s"%(mpts_params["module_inppath"],module_anlys["inpname"])
        if not os.path.exists(inpfile):
            raise MPTSRunError("Analysis module input file named %s does not exist"%inpfile)
        mpts_params["module_anlys"]["inpfile"] = inpfile
    if "vars" not in module_anlys.keys():
        mpts_params["module_anlys"]["vars"] = dict()

    #ml modules
    if "module_ml" not in mpts_params.keys():
        raise MPTSRunError("Argument of dictionary type named module_ml not found in MPTS input file.")
    module_ml = mpts_params["module_ml"]
    if "name" not in module_ml.keys():
        raise MPTSRunError("Argument named name in module_ml dictionary not found in MPTS input file.")
    else:
        pyfile = "%s/%s.py"%(mpts_params["module_path"],module_ml["name"])
        if not os.path.exists(pyfile):
            raise MPTSRunError("ML module py file named %s does not exist"%pyfile)
        mpts_params["module_ml"]["pyfile"] = pyfile
    if "inpname" not in module_ml.keys():
        raise MPTSRunError("Argument named inpname in module_ml dictionary not found in MPTS input file.")
    else:
        inpfile = "%s/%s"%(mpts_params["module_inppath"],module_ml["inpname"])
        if not os.path.exists(inpfile):
            raise MPTSRunError("ML module input file named %s does not exist"%inpfile)
        mpts_params["module_ml"]["inpfile"] = inpfile
    if "vars" not in module_ml.keys():
        mpts_params["module_ml"]["vars"] = dict()


    #errors in library data file
    df = pd.read_csv(mpts_params["libdatafile"],delimiter = " ", dtype = str, na_filter = False, nrows = 1)
    datadict = df.to_dict('list')
    ndescs = 0
    for key in datadict.keys():
        if "desc_" in key:
            ndescs = ndescs + 1
    if ndescs == 0:
        raise MPTSRunError("Number of descriptors in the Library data file is 0. No keys with desc_ found.")
    ntimes = 0
    for key in datadict.keys():
        if "id" in key:
            ntimes += 1
    if ntimes == 0:
        raise MPTSRunError("dataid key not found in the data file of initial candidate batch")
    if ntimes > 1:
        raise MPTSRunError("dataid key found multiple times in the data file of initial candidate batch")

    #errors in library init file
    df = pd.read_csv(mpts_params["libinitdatafile"],delimiter = " ", dtype = str, na_filter = False, nrows = 1)
    datadict = df.to_dict('list')
    ndescs = 0
    for key in datadict.keys():
        if "desc_" in key:
            ndescs = ndescs + 1
    if ndescs == 0:
        raise MPTSRunError("Number of descriptors in the data file of initial candidate batch is 0. No keys with desc_ found.")
    ntimes = 0
    for key in datadict.keys():
        if "id" in key:
            ntimes += 1
    if ntimes == 0:
        raise MPTSRunError("dataid key not found in the data file of initial candidate batch")
    if ntimes > 1:
        raise MPTSRunError("dataid key found multiple times in the data file of initial candidate batch")

    print('No errors found.')
    #print("********************************************************************")
    return mpts_params


'''
load modules
mpts_params : input dict
'''
def load_modules(mpts_params):
    #copy args
    module_path = mpts_params["module_path"]
    sys.path.insert(1,module_path)

    #simulation module
    sim_name = mpts_params["module_sim"]["name"]
    simmodule = importlib.import_module(sim_name)
    mpts_params["module_sim"]["pymodule"] = simmodule
    print("Loading simulation module ... %s"%simmodule)

    #analysis module
    anlys_name = mpts_params["module_anlys"]["name"]
    anlysmodule = importlib.import_module(anlys_name)
    mpts_params["module_anlys"]["pymodule"] = anlysmodule
    print("Loading analysis module ... %s"%anlysmodule)

    #simulation module
    ml_name = mpts_params["module_ml"]["name"]
    mlmodule = importlib.import_module(ml_name)
    mpts_params["module_ml"]["pymodule"] = mlmodule
    print("Loading ML module ... %s"%mlmodule)

    return mpts_params

'''
concatenate two files into one
infile1: inputfile1
infile2: inputfile2
outfile: outfile
'''
def concatenatefiles(infile1,infile2,outfile):
    #read infile1
    df1 = pd.read_csv(infile1,delimiter = " ", dtype = str, na_filter = False)
    indict1 = df1.to_dict('list')
    
    #read infile2
    df2 = pd.read_csv(infile2,delimiter = " ", dtype = str, na_filter = False)
    indict2 = df2.to_dict('list')

    #concatenate
    outdict = {key: indict1[key] + indict2[key] for key in indict1}
    df = pd.DataFrame.from_dict(outdict)
    df.to_csv(outfile,index = False, sep = " ")

'''
Running many iterations of machine learning
mpts_params : input dict
'''
def run_mliter(mpts_params):
    #copy args
    out_path = mpts_params["out_path"]
    mliter = mpts_params["mliter"]
    libinitdatafile = mpts_params["libinitdatafile"]
    libdatafile = mpts_params["libdatafile"]

    #output dictionary 
    odict = dict()
    odict["prpslfile"] = "NULL"

    #change directory
    os.chdir(out_path)

    #an integer run flag
    # 0 in progress
    # 1 means succesful completion, 
    # 2 means training stopped in between
    # -1 something fish happened
    odict["runflag"] = 0

    #counting the number of proposal files
    prpslfiles = glob.glob("prpsl*.txt")
    nprpsl = len(prpslfiles)

    #check points with maximum iterations
    if nprpsl >= mliter + 1:
        odict["runflag"] = 1
        odict["iterid"] = nprpsl - 1        
        odict["prpslfile"] = "prpsl_%s.txt"%odict["iterid"]
        print("Maximum number of ML iterations completed, %d "%(nprpsl-1))
        print("Returning.....")
        print("********************************************************************")
        return odict

    #looping
    if nprpsl == 0:
        prpslfile = "prpsl_0.txt"
        print("ML Iteration 0 : Proposal file not found. Copying the initial batch file as first proposal file")
        shutil.copyfile(libinitdatafile,prpslfile)
        odict["iterid"] = -1
        odict["prpslfile"] = prpslfile
        return odict
    else:
        #generate the training file by concatenating the prpsl file
        iterid = nprpsl - 1
        prpslfile = "prpsl_%s.txt"%iterid
        trainfile = "train_%s.txt"%iterid

        #existence of train file
        if os.path.exists(trainfile):
            print("ML Iteration %s : Something fishy has happened here. Training file already exists before it is being written by the MPTS"%iterid)
            odict["runflag"] = -1
            odict["prpslfile"] = prpslfile
            odict["iterid"] = iterid
            print("********************************************************************")
            return odict

        # check for NA in proposal file
        # read proposal file
        df = pd.read_csv(prpslfile,delimiter = " ", dtype = str, na_filter = False)
        tmpdict = df.to_dict('list')
        tmparr = np.asarray(tmpdict["Output"]).astype("str")
        if len(np.where(tmparr == "NA")[0]) != 0:
            print("ML Iteration %s : Proposal file still contain NA files. Data required"%iterid)
            odict["runflag"] = 2
            odict["iterid"] = iterid
            odict["prpslfile"] = prpslfile
            return odict

        #generate train_iterid.txt by combining prpslfile and previous train file
        if iterid == 0:
            shutil.copyfile(prpslfile,trainfile)
        else:
            trainprevfile = "train_%d.txt"%(iterid-1)
            concatenatefiles(trainprevfile,prpslfile,trainfile)

        #run ML module
        print("ML Iteration %s : Running ML module."%iterid)
        libfile = libdatafile
        curr_train_file = trainfile
        next_batch_size = mpts_params["nextsize"]
        next_batch_file = "prpsl_%d.txt"%(iterid+1)
        inp_file = mpts_params["module_ml"]["inpfile"]
        vardict = mpts_params["module_ml"]["vars"]
        if "initseed" in vardict.keys():
            vardict["rseed"] = vardict["initseed"] + iterid
        mpts_params["module_ml"]["pymodule"].main(libfile,curr_train_file,next_batch_size,next_batch_file,inp_file,vardict)
        
        #odict
        odict["iterid"] = iterid + 1
        return odict

'''
Running main function of MPTS
mpts_params : input dict
'''
def run_MPTSmain(mpts_params):

    #copy args
    out_path = mpts_params["out_path"]
    workpath = mpts_params["workpath"]
    
    #odict
    odict = run_mliter(mpts_params)
    os.chdir(workpath)
    runflag = odict["runflag"]
    iterid = odict["iterid"]
    prpslfile = odict["prpslfile"]

    #runflag status
    if runflag == 1:
        print("All ML iterations finished.")
    if runflag == -1:
        print("MPTS main is abrupty finished at ML iteration %s"%iterid)
    if runflag == 2:
        print("ML Iteration %s : Running simulaton module."%iterid)
    
        #submitting simulations
        batch_file = "%s/%s"%(out_path,prpslfile)
        inp_file = mpts_params["module_sim"]["inpfile"]
        vardict = mpts_params["module_sim"]["vars"]
        simflag = mpts_params["module_sim"]["pymodule"].main(batch_file,inp_file,vardict)
        if simflag == 1:
            print("ML Iteration %s : Simulations finished. "%iterid)
        else:
            print("ML Iteration %s : Simulations still in progress."%iterid)

        #analyzing them
        print("ML Iteration %s : Running analysis module."%iterid)
        batch_file = "%s/%s"%(out_path,prpslfile)
        inp_file = mpts_params["module_anlys"]["inpfile"]
        vardict = mpts_params["module_anlys"]["vars"]
        mpts_params["module_anlys"]["pymodule"].main(batch_file,inp_file,vardict)
        print("ML Iteration %s : Analysis finished."%iterid)
        
    mpts_params["runflag"] = runflag

    return mpts_params


'''
Writing of dictionary in yaml
mpts_params : input dict
filename : Name of file
'''
def writelog(mpts_params,filename):
    print("********************************************************************")
    print("Writing mpts_params in %s"%filename)
    fulldict = {'mpts' : mpts_params}
    with open(filename,'w') as outfile:
        yaml.dump(fulldict,outfile)
    print("********************************************************************")




#---------------------------------------------------------------------------------
# MAIN FUNCTION
#---------------------------------------------------------------------------------     
def main(mpts_params):

    print("********************************************************************")
    print("********************************************************************")
    print("RUNNING MODEL PARAMETER TARGETED SEARCH (MPTS)")

    #workpath
    workpath = os.getcwd()

    #set defaults
    mpts_params = defaults(mpts_params)

    #error check
    mpts_params = errorcheck(mpts_params)

    #load modules
    mpts_params = load_modules(mpts_params)

    #running MPTS iteration
    niter = mpts_params["niter"]
    sleeptime = mpts_params["sleeptime"]
    i = 0
    isRun = True
    while isRun:

        print("MPTS ITERATION %d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"%(i+1))

        #run mlstages
        mpts_params = run_MPTSmain(mpts_params)
        print("SLEEPING FOR %s seconds"%sleeptime)
        time.sleep(sleeptime)
        os.chdir(workpath)

        #update i
        i = i + 1

        #stop when reaching maximum MPTSiterations
        if i == niter:
            isRun = False
            print("Maximum iterations of running MPTS reached...")

        #stop when ml is finished or abruptly ended
        runflag = mpts_params["runflag"]
        if runflag == 1 or runflag == -1:
            isRun = False


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
