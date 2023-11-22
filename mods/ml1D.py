import yaml
import pandas as pd
import os
import numpy as np
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import Matern

'''
Reading the data file containing descriptors and output
datafile : name of data file
returns the fulldict and indict containing descriptors and output
'''
def readdatafile(datafile):
    df = pd.read_csv(datafile,delimiter = " ", dtype = str, na_filter = False)
    fulldict = df.to_dict('list')
    indict = dict()
    for key,value in fulldict.items():
        if "desc_" in key:
            indict[key] = value
        if  key == "Output":
            indict[key] = value
        if  key == "id":
            indict[key] = value
    return fulldict,indict


'''
parse the dictionary and return X,y,idlst as list
indict : Input dictionary
'''
def parseindict(indict):
    Xlst = []
    ylst = []
    hypdlst = []
    idlst = None
    for key,value in indict.items():
        if "desc_" in key:
            Xlst.append(value)
        if key == "Output":
            ylst.append(value)
        if key == "id":
            idlst = value
    X = np.asarray(Xlst).T
    y = np.asarray(ylst).T
    return X,y,idlst


'''
gaussian process regression model
Xtrain : training dataset X
ytrain : training dataset,y
inpparams : input parameters for the model
'''
def gpr_model(Xtrain,ytrain,Xpredict,inpparams):
    #float
    Xtrain = Xtrain.astype(float)
    ytrain = ytrain.astype(float)
    Xpredict = Xpredict.astype(float)

    #inputparams
    rseed = inpparams["rseed"]
    length_scale = float(inpparams["length_scale"])
    normalize_y = inpparams["normalize_y"]
    length_scale_bounds = (float(inpparams["length_scale_bounds_low"]),float(inpparams["length_scale_bounds_up"]))
    nu = inpparams["nu"]
    n_restarts_optimizer = inpparams["n_restarts_optimizer"]
    alpha = float(inpparams["alpha"])

    #gpr_model
    rng = np.random.RandomState(rseed)
    kernel = 1 * Matern(nu = nu, length_scale=length_scale, length_scale_bounds=length_scale_bounds)
    model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=n_restarts_optimizer, random_state = rng, normalize_y = True, alpha = alpha)
    model.fit(Xtrain, ytrain)

    return model


'''
Exploitation strategy under GPR
model : model
Xpredict : input dataset for prediction
next_batch_size : size of next batch
'''
def exploitation(model,Xpredict,next_batch_size):
    #copy args
    Xpredict = Xpredict.astype(float)
    
    #model predict
    mean, std = model.predict(Xpredict,return_std = True)
    
    #take the maximum mean value
    mean = mean[:,0]
    idx = np.argpartition(mean, -next_batch_size)
    idx = idx[-next_batch_size:]
    return idx


'''
Exploration strategy under GPR
model : model
Xpredict : input dataset for prediction
next_batch_size : size of next batch
'''
def exploration(model,Xpredict,next_batch_size):
    #copy args
    Xpredict = Xpredict.astype(float)
    
    #model predict
    mean, std = model.predict(Xpredict,return_std = True)
    
    #take the maximum standard deviation value
    std = std[:]
    idx = np.argpartition(std, -next_batch_size)
    idx = idx[-next_batch_size:]
    return idx

'''
Balanced Exploitation/Exploration strategy under GPR
model : model
Xpredict : input dataset for prediction
next_batch_size : size of next batch
ymax : Maximum y value obtained till now
'''
def bee(model,Xpredict,next_batch_size,ymax):
    #copy args
    Xpredict = Xpredict.astype(float)

    #model predict
    mean, std = model.predict(Xpredict,return_std = True)

    #take the maximum of expected improvement
    mean  = mean[:,0]
    pix = np.zeros((len(mean),))
    eix = np.zeros((len(mean),))
    for k in range(len(mean)):
        s = std[k]
        m = mean[k]
        if s == 0.0:
            eix[k] = 0.0
        else:
            z = (m - ymax)/s
            eix[k] = (m - ymax)*norm.cdf(z) + s*norm.pdf(z)
    idx = np.argpartition(eix, -next_batch_size)
    idx = idx[-next_batch_size:]
    return idx


'''
Store X,y with keys associated with indict in a file
X : numpy array of descriptors
y : numpy array of output
indict : dictionary with key information
outfile : name of file to be written
'''
def writeXy_with_keys(X,y,indict,outfile,idtrain = None, fulldict = None):
    outdict = dict()
    desc_index = 0
    for key in indict.keys():
        if "desc_" in key:
            outdict[key] = X[:,desc_index]
            desc_index += 1
        if key == "Output":
            outdict[key] = y[:,0]

    #also add training ids
    if idtrain is not None:
        outdict["id"] = idtrain

    #also add the attributes other than desc and output
    if fulldict is not None:
        for key in fulldict.keys():
            if key not in indict.keys():
                lst = []
                value = fulldict[key]
                for i in range(len(idtrain)):
                    lst.append(value[idtrain[i]])
                outdict[key] = lst
                
    #write 
    df = pd.DataFrame.from_dict(outdict)
    df.to_csv(outfile,index = False, sep = " ")


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

    #read the libfile
    fulldict,indict = readdatafile(libfile)
    X,y,idall = parseindict(indict)

    #read the curr_train_file
    fulldicttrain,indicttrain = readdatafile(curr_train_file)
    Xtrain,ytrain,idtrain = parseindict(indicttrain)
    
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
    
    #train the model
    model = gpr_model(Xtrain,ytrain,Xpredict,mlparams)

    #active learning optimization
    opt_technique = mlparams["opt_technique"]
    if opt_technique == "EXPT":
        idprpsl = exploitation(model,Xpredict,next_batch_size)
        idnext = idpredict[idprpsl]
    if opt_technique == "EXPR":
        idprpsl = exploration(model,Xpredict,next_batch_size)
        idnext = idpredict[idprpsl]
    if opt_technique == "BEE":
        ymax = np.max((ytrain.astype(float)))
        idprpsl = bee(model,Xpredict,next_batch_size,ymax)
        idnext = idpredict[idprpsl]


    #write the new proposal file
    Xnext = X[idnext]
    ynext = y[idnext]
    writeXy_with_keys(Xnext,ynext,indict,next_batch_file,idtrain = idnext, fulldict = fulldict)
    



