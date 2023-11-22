import yaml
import pandas as pd
import os
import numpy as np

'''
batch_file : Library file of a batch of canddiates
inp_file : Input file to the simulation module
vardict : Dictionary to add/overwrite the variables in inp_file
'''
def main(batch_file,inp_file,vardict):
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

    #run the model
    cand_namelst = data['cand_name']
    xlst = data['desc_x']
    sim_path = os.path.abspath(simparams["sim_path"])
    cand_pathlst = []
    cand_statuslst = []
    for i in range(len(cand_namelst)):
        cand_path  = "%s/%s"%(sim_path,cand_namelst[i])
        #if data path already exists for a given candidate, dont run the model
        if os.path.exists(cand_path):
            cand_status = 1
        else:
            #make directory
            os.makedirs(cand_path,exist_ok = True)
            #simulate the model
            xval = float(xlst[i])
            yval = np.exp(-(xval - 2) ** 2) + np.exp(-(xval - 6) ** 2 / 10) + 1/ (xval ** 2 + 1)
            #store the yval
            f = open("%s/op.txt"%cand_path,'w')
            f.write("y : %s\n"%yval)
            f.close()
            cand_status = 1
        cand_statuslst.append(cand_status)
        cand_pathlst.append(cand_path)

    #check the status
    flag = 1   
    for i in range(len(cand_statuslst)):
        cand_status = cand_statuslst[i]
        if cand_status == 0:
            flag = -1

    #if simulations finished, rewrite the batch_file
    if flag == 1:
        data["cand_path"] = cand_pathlst
        data["cand_status"] = cand_statuslst
        df = pd.DataFrame.from_dict(data)
        df.to_csv(batch_file,index = False, sep = " ")
        
    #return flag
    return flag
        
