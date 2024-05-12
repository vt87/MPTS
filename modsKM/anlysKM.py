import yaml
import pandas as pd
import os
import numpy as np
import time
import argparse
import sys
import shutil
import importlib
sys.path.append(os.path.dirname(__file__))
import anlysKMutils

def initialize():
    parser = argparse.ArgumentParser(description='Options to do analysis')

    parser.add_argument('--inp_file',  type=str, help='Input file to the analysis  module',default = "inp.txt")
    parser.add_argument('--batch_file',  type=str, help='Library file of batch of candidates',default = "inp.txt")

    return parser

class AnlysKMError(Exception):

    def __init__(self, message):
        super().__init__(message)

'''
batch_file : Library file of a batch of canddiates
inp_file : Input file to the analysis module
vardict : Dictionary to add/overwrite the variables in inp_file
'''
def main(batch_file,inp_file,vardict):
    #workpath and flag
    workpath = os.getcwd()
    flag = -1

    #read inp_file
    with open(inp_file, 'r') as yml:
        cfg = yaml.load(yml, Loader=yaml.FullLoader)
        anlysparams = cfg['anlys']

    #add/overwrite the variables in anlysarams using vardict
    for key in vardict.keys():
        anlysparams[key] = vardict[key]

    #read batch_file
    df = pd.read_csv(batch_file,delimiter = " ", dtype = str, na_filter = False)
    data = df.to_dict('list')

    #flag 1
    OPstyle = anlysparams["OPstyle"]
    anlyslogdir = "%s/%s"%(os.path.abspath(anlysparams["anlyslogdir"]),OPstyle)
    bfile = "%s/batch.txt"%anlyslogdir
    if os.path.exists(bfile):
        flag = 1
        #replace batch_file
        shutil.copyfile(bfile,batch_file)
        return flag

    #set fortranpath
    if OPstyle == "OPSTEINHARDT_AG":
        oppath = "%s/aux/fortranBOO"%os.path.dirname(os.path.abspath(__file__))
        tmppath = "%s/fortranBOO"%(anlyslogdir)
        if not os.path.exists(tmppath):
            shutil.copytree(oppath,tmppath)
        anlysparams["fortranPATH"] = tmppath
    else:
        raise AnlysKMError("Invalid OP style named %s"%OPstyle)

    #fill Output column
    os.makedirs(anlyslogdir,exist_ok = True)
    cand_pathlst = data['cand_path']
    cand_statuslst = data["cand_status"]
    Outputlst = []
    for i in range(len(cand_pathlst)):
        cand_status = int(cand_statuslst[i])
        if cand_status == -1:
            Output = "FAIL"
        else:
            datafile = "%s/in.data"%cand_pathlst[i]
            trjfile = "%s/trajfinal.xtc"%cand_pathlst[i]
            if OPstyle == "OPSTEINHARDT_AG":
                anlysparams["outclusfile"] = "%s/opst_log.txt"%cand_pathlst[i]
                Output = anlysKMutils.getOPSTEINHARDT_AG_from_trj(datafile,trjfile,anlysparams)
            else:
                raise AnlysKMError("Invalid OP style named %s"%OPstyle)
        Outputlst.append(Output)
    data["Output"] = Outputlst
    df = pd.DataFrame.from_dict(data)
    df.to_csv(bfile,index = False, sep = " ")

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
