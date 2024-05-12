import yaml
import pandas as pd
import os
import numpy as np
import time
import argparse
import sys
sys.path.append(os.path.dirname(__file__))
import mlutils
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import Matern
import warnings
warnings.filterwarnings("error")

def initialize():
    parser = argparse.ArgumentParser(description='Options to do analysis')

    parser.add_argument('--libfile',  type=str, help='Library File',default = "inp.txt")
    parser.add_argument('--curr_train_file',  type=str, help='Library file containing data gathered from all previous iterations',default = "inp.txt")
    parser.add_argument('--next_batch_size',  type=int, help='Next batch size',default = 20)
    parser.add_argument('--next_batch_file',  type=str, help='Next batch file',default = "inp.txt")
    parser.add_argument('--inp_file',  type=str, help='Input file to the ML  module',default = "inp.txt")
    return parser

class MLKMError(Exception):
    def __init__(self, message):
        super().__init__(message)

'''
libfile : Library file of all candidates
curr_train_file : Library file containing data gathered from all previous iterations
next_batch_size : Batch size for the next batch of candidates
next_batch_file : Library file for the next batch of candidates
inp_file : Input file to the ML module
vardict : Dictionary to add/overwrite the variables in inp_file
'''
def main(libfile,curr_train_file,next_batch_size,next_batch_file,inp_file,vardict):
    #read inp_file
    with open(inp_file, 'r') as yml:
        cfg = yaml.load(yml, Loader=yaml.FullLoader)
        mlparams = cfg['ml']

    #add/overwrite the variables in mlparams using vardict
    for key in vardict.keys():
        mlparams[key] = vardict[key]

    #distance parameter
    if "hyp" in mlparams.keys():
        hyp_params = mlparams["hyp"]
        isHyp = True
        hyp_name = hyp_params["name"]
        if hyp_name == "dist":
            hyp_key = hyp_params["key"]
            hyp_s = hyp_params["scale"]
        else:
            raise MLKMError("Invalid hyperparameter name %s"%hyp_name)
    else:
        isHyp = False
        hyp_key = None

    #read the libfile
    fulldict,indict = mlutils.readdatafile(libfile, hyp_key = hyp_key)
    X,y,hypd,idall = mlutils.parseindict(indict, hyp_key = hyp_key)

    #scalar
    if "isScaled" not in mlparams.keys():
        mlparams["isScaled"] = False
    if mlparams["isScaled"]:
        mlparams["scaler"] = mlutils.scale_features(X)
    else:
        mlparams["scaler"] = None

    #read the curr_train_file
    fulldicttrain,indicttrain = mlutils.readdatafile(curr_train_file, hyp_key = hyp_key)
    Xtrain,ytrain,hypdtrain,idtrain = mlutils.parseindict(indicttrain, hyp_key = hyp_key)

    #add the assertion
    orderflag = mlutils.assert_order_in_keys(indict,indicttrain)
    if not orderflag :
        raise MLV1Error("Order of descriptors in training and full data set is not correct")
    
    #predict data points
    idpredict = np.where(~np.in1d(idall,idtrain))[0]
    Xpredict = X[idpredict]

    #remove the failed points
    failid = np.where(ytrain == "FAIL")[0]
    successid = np.where(ytrain != "FAIL")[0]
    if len(failid) == len(ytrain):
        raise Exception("ML module : Cannot build the model. Every simulation is failed.")
    if len(failid) > 0:
        print("ML module : Out of %s simulations, %s simulations failed. removing failed points"%(len(ytrain),len(failid)))
        Xtrain = Xtrain[successid]
        ytrain = ytrain[successid]
        if hypd is not None:
            hypdtrain = hypdtrain[successid]    

    #train the model
    model = mlutils.gpr_model(Xtrain,ytrain,Xpredict,mlparams,hypdtrain = hypdtrain)

    #active learning optimization
    opt_technique = mlparams["opt_technique"]
    if opt_technique == "EXPT":
        idprpsl = mlutils.exploitation(model,Xpredict,next_batch_size,mlparams)
        idnext = idpredict[idprpsl]
    if opt_technique == "EXPR":
        idprpsl = mlutils.exploration(model,Xpredict,next_batch_size,mlparams)
        idnext = idpredict[idprpsl]
    if opt_technique == "BEE":
        if hypdtrain is not None:
            if mlparams["hyp"]["scale"] == "mult" :
                ymax = np.max((ytrain.astype(float)*hypdtrain.astype(float)))
        else:
            ymax = np.max((ytrain.astype(float)))
        idprpsl = mlutils.bee(model,Xpredict,next_batch_size,ymax,mlparams)
        idnext = idpredict[idprpsl]

    #write the new proposal file
    Xnext = X[idnext]
    ynext = y[idnext]
    if hypd is not None:
        hypdnext = hypd[idnext]
    else:
        hypdnext = None
    mlutils.writeXy_with_keys(Xnext,ynext,indict,next_batch_file,idtrain = idnext, fulldict = fulldict, hyp_key = hyp_key, hypd =  hypdnext)

# Execute main function if this file is the starting point
if __name__ == '__main__':
    start_time = time.time()

    #Read arguments
    args = initialize().parse_args()    
    params = vars(args)
    vardict = dict()
    flag = main(params["libfile"],params["curr_train_file"],params["next_batch_size"],params["next_batch_file"],params["inp_file"],dict())

    print("--- %s seconds ---" % (time.time() - start_time))
