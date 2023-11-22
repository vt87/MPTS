import yaml
import pandas as pd
import os
import numpy as np

'''
batch_file : Library file of a batch of canddiates
inp_file : Input file to the analysis module
vardict : Dictionary to add/overwrite the variables in inp_file
'''
def main(batch_file,inp_file,vardict):
    #read inp_file
    with open(inp_file, 'r') as yml:
        cfg = yaml.load(yml, Loader=yaml.FullLoader)
        simparams = cfg['anlys']

    #add/overwrite the variables in anlysarams using vardict
    for key in vardict.keys():
        anlysparams[key] = vardict[key]

    #read batch_file
    df = pd.read_csv(batch_file,delimiter = " ", dtype = str, na_filter = False)
    data = df.to_dict('list')

    #fill Output column
    cand_pathlst = data['cand_path']
    cand_statuslst = data["cand_status"]
    Outputlst = []
    for i in range(len(cand_pathlst)):
        cand_status = int(cand_statuslst[i])
        if cand_status == -1:
            Output = "FAIL"
        else:
            opfile = "%s/op.txt"%cand_pathlst[i]
            f = open(opfile,'r')
            lst = f.readlines()
            f.close()
            Output = float(lst[0].split(":")[1])
        Outputlst.append(Output)
    data["Output"] = Outputlst
    df = pd.DataFrame.from_dict(data)
    df.to_csv(batch_file,index = False, sep = " ")
