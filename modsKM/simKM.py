import yaml
import pandas as pd
import os
import numpy as np
import time
import argparse
import sys
import shutil
sys.path.append(os.path.dirname(__file__))
import simutils
import simbruteKM

def initialize():
    parser = argparse.ArgumentParser(description='Options to do analysis')

    parser.add_argument('--inp_file',  type=str, help='Input file to the simulation module',default = "inp.txt")
    parser.add_argument('--batch_file',  type=str, help='Library file of batch of candidates',default = "inp.txt")

    return parser

class SimKMError(Exception):

    def __init__(self, message):
        super().__init__(message)

'''
batch_file : Library file of a batch of canddiates
inp_file : Input file to the simulation module
vardict : Dictionary to add/overwrite the variables in inp_file
'''
def main(batch_file,inp_file,vardict):

    #workpath and flag
    workpath = os.getcwd()
    flag = -1

    #read inp_file
    with open(inp_file, 'r') as yml:
        cfg = yaml.load(yml, Loader=yaml.FullLoader)
        simparams = cfg['sim']

    #add/overwrite the variables in simparams using vardict
    for key in vardict.keys():
        simparams[key] = vardict[key]
        
    #read batch_file
    df = pd.read_csv(batch_file,delimiter = " ", dtype = str, na_filter = False)
    data = df.to_dict('list')

    #flag 1
    simlogdir = os.path.abspath(simparams["simlogdir"])
    sim_path = os.path.abspath(simparams["sim_path"])
    bfile = "%s/batch.txt"%simlogdir
    if os.path.exists(bfile):
        flag = 1
        #replace batch_file
        shutil.copyfile(bfile,batch_file)
        return flag

    #first write siminp.yaml in simlogdir
    if os.path.exists(simlogdir):
        return flag

    #job submit
    os.makedirs(simlogdir,exist_ok = True)
    simparams["simlogdir"] = simlogdir
    simparams["sim_path"] = sim_path
    inpl_file = "%s/siminp.yaml"%(simlogdir)
    fulldict = {'sim' : simparams}
    with open(inpl_file,'w') as outfile:
        yaml.dump(fulldict,outfile)
    submit_style = simparams["submit_style"]
    simtmplname = simparams["simtmplname"]
    if submit_style == "local" :
        if simtmplname == "simbruteKM":
            simbruteKM.main(batch_file,inpl_file)
    if submit_style == "cluster" :
        jobfile = "%s/job_local.sh"%simlogdir
        clus = simparams["clus"]
        py_name = "%s/%s.py --inp_file %s --batch_file %s --run_style recursive"%(os.path.dirname(__file__),simtmplname,inpl_file,os.path.abspath(batch_file))
        simutils.local_py(jobfile,1,clus,py_name)
        os.chdir(simlogdir)
        cmd = "sbatch job_local.sh"
        os.system(cmd)
        os.chdir(workpath)


    return flag
                                        

# Execute main function if this file is the starting point
if __name__ == '__main__':
    start_time = time.time()

    #Read arguments
    args = initialize().parse_args()    
    params = vars(args)
    vardict = dict()
    flag = main(params["batch_file"],params["inp_file"],dict())
    print(flag)

    print("--- %s seconds ---" % (time.time() - start_time))
