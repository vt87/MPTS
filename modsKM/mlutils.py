import yaml
import pandas as pd
import os
import numpy as np
import time
import sys
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import Matern
from sklearn.preprocessing import MinMaxScaler
import warnings
warnings.filterwarnings("error")

class MLUtilsError(Exception):
    """ Raised if invalid thermodynamic ensemble is passed """

    def __init__(self, message):
        super().__init__(message)


'''
assertion
check whether desc_ keys are in same order
'''
def assert_order_in_keys(dict1,dict2):
    lst1 = []
    for key in dict1.keys():
        if "desc_" in key:
            lst1.append(key)
    lst2 = []
    for key in dict2.keys():
        if "desc_" in key:
            lst2.append(key)
    return lst1 == lst2


'''
minmaxscale
'''
def scale_features(X):
    X = X.astype(float)
    scaler = MinMaxScaler()
    scaler.fit(X)
    #print(scaler.data_max_)
    return scaler

'''
Reading the data file containing descriptors and output
datafile : name of data file
hyp_key : hyperparameter key
returns the fulldict and indict containing descriptors and output
'''
def readdatafile(datafile, hyp_key = None):
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
        if hyp_key is not None:
            if key == hyp_key :
                indict[key] = value
    return fulldict,indict


'''
parse the dictionary and return X,y,idlst as list
indict : Input dictionary
'''
def parseindict(indict, hyp_key = None):
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
        if hyp_key is not None:
            if key == hyp_key :
                hypdlst.append(value)
    X = np.asarray(Xlst).T
    y = np.asarray(ylst).T
    #distance list
    if len(hypdlst) == 0:
        hypd = None
    else:
        hypd = np.asarray(hypdlst).T
    return X,y,hypd,idlst


'''
gaussian process regression model
Xtrain : training dataset X
ytrain : training dataset,y
inpparams : input parameters for the model
'''
def gpr_model(Xtrain,ytrain,Xpredict,inpparams,hypdtrain = None):
    #float
    Xtrain = Xtrain.astype(float)
    ytrain = ytrain.astype(float)
    Xpredict = Xpredict.astype(float)

    #distance hyperparameter
    if hypdtrain is not None:
        if inpparams["hyp"]["scale"] == "mult" :
            hypdtrain = hypdtrain.astype(float)
            ytrain *= hypdtrain
            
    #scale
    #print(Xtrain,ytrain)
    if inpparams["isScaled"]:
        Xtrain = inpparams["scaler"].transform(Xtrain)
    #print(Xtrain)

    #kernel
    kernelstyle  = inpparams["kernel"]

    if kernelstyle == "Matern":
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
        #for hyperparameter in kernel.hyperparameters: print(hyperparameter)
        #params = kernel.get_params()
        #for key in sorted(params): print("%s : %s" % (key, params[key]))
        #sys.exit(-1)
        #length_scale_vec = kernel.get_params()['k2__length_scale']
        #print('length_scale_vec=',length_scale_vec)

        model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=n_restarts_optimizer, random_state = rng, normalize_y = normalize_y, alpha = alpha)
        model.fit(Xtrain, ytrain)

    elif kernelstyle == "RBF":
        #inputparams
        rseed = inpparams["rseed"]
        length_scale = float(inpparams["length_scale"])
        normalize_y = inpparams["normalize_y"]
        length_scale_bounds = (float(inpparams["length_scale_bounds_low"]),float(inpparams["length_scale_bounds_up"]))
        n_restarts_optimizer = inpparams["n_restarts_optimizer"]
        alpha = float(inpparams["alpha"])

        #gpr_model
        rng = np.random.RandomState(rseed)
        kernel = 1 * RBF(length_scale=length_scale, length_scale_bounds=length_scale_bounds)

        model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=n_restarts_optimizer, random_state = rng, normalize_y = normalize_y, alpha = alpha)
        model.fit(Xtrain, ytrain)

    else:
        raise MLUtilsError("Invalid Kernel for GPR")

    print("GPR  score :",model.score(Xtrain,ytrain))

    return model


'''
Exploitation strategy under GPR
model : model
Xpredict : input dataset for prediction
next_batch_size : size of next batch
mlparams : pass additional parameters
'''
def exploitation(model,Xpredict,next_batch_size,mlparams):
    #copy args
    Xpredict = Xpredict.astype(float)
    
    #scale
    if mlparams["isScaled"]:
        Xpredict = mlparams["scaler"].transform(Xpredict)

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
mlparams : pass additional parameters
'''
def exploration(model,Xpredict,next_batch_size,mlparams):
    #copy args
    Xpredict = Xpredict.astype(float)

    #scale
    if mlparams["isScaled"]:
        Xpredict = mlparams["scaler"].transform(Xpredict)
    
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
mlparams : pass additional parameters
'''
def bee(model,Xpredict,next_batch_size,ymax,mlparams):
    #copy args
    Xpredict = Xpredict.astype(float)

    #scale
    if mlparams["isScaled"]:
        Xpredict = mlparams["scaler"].transform(Xpredict)

    #model predict
    mean, std = model.predict(Xpredict,return_std = True)
    #print(mean)

    #take the maximum of expected improvement
    mean  = mean[:,0]
    pix = np.zeros((len(mean),))
    eix = np.zeros((len(mean),))
    #std[0] = 0.0
    # for k in range(len(mean)):
    #     s = std[k]
    #     m = mean[k]
    #     if s == 0.0:
    #         eix[k] = 0.0
    #     else:
    #         z = (m - ymax)/s
    #         eix[k] = (m - ymax)*norm.cdf(z) + s*norm.pdf(z)
    #numpy where. Using the errstate function to suppress divide by 0 warnings
    with np.errstate(divide='ignore'):
        eix = np.where(std == 0.0,0.0,((mean-ymax)*norm.cdf((mean-ymax)/std)) + (std*norm.pdf((mean-ymax)/std)))
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
def writeXy_with_keys(X,y,indict,outfile,idtrain = None, fulldict = None, hyp_key = None, hypd = None):
    outdict = dict()
    desc_index = 0
    for key in indict.keys():
        if "desc_" in key:
            outdict[key] = X[:,desc_index]
            desc_index += 1
        if key == "Output":
            outdict[key] = y[:,0]
        if hyp_key is not None:
            if key == hyp_key:
                outdict[key] = hypd[:,0]

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
