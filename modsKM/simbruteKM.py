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

def initialize():
    parser = argparse.ArgumentParser(description='Options to do analysis')

    parser.add_argument('--inp_file',  type=str, help='Input file to the simulation module',default = "inp.txt")
    parser.add_argument('--batch_file',  type=str, help='Library file of batch of candidates',default = "inp.txt")
    parser.add_argument('--run_style',  type=str, help='run style',default = "single")

    return parser

class SimBruteKMError(Exception):

    def __init__(self, message):
        super().__init__(message)

'''
batch_file : Library file of a batch of canddiates
inp_file : Input file to the simulation module
'''
def main(batch_file,inp_file):
    #workpath
    workpath = os.getcwd()

    #read inp_file
    with open(inp_file, 'r') as yml:
        cfg = yaml.load(yml, Loader=yaml.FullLoader)
        simparams = cfg['sim']

    #read batch_file
    df = pd.read_csv(batch_file,delimiter = " ", dtype = str, na_filter = False)
    data = df.to_dict('list')

    #lmp_scripts
    lmp_scripts = simparams["lmp_scripts"]
    lmp_path = simparams["lmp_path"]
    for i in range(len(lmp_scripts)):
        l = "%s/%s"%(lmp_path,lmp_scripts[i])
        if not os.path.exists(l):
            raise SimBruteKMError("File named %s does not exist"%l)

    #lmp_exec
    lmp_exec = simparams["lmp_exec"]
    if not os.path.exists(lmp_exec):
        raise SimBruteKMError("Lammps executable File named %s does not exist"%lmp_exec)


    #look for the lattice descriptor
    cand_namelst = data['cand_name']
    sim_path = os.path.abspath(simparams["sim_path"])
    cand_pathlst = ["NULL"]*len(cand_namelst)
    cand_statuslst = [0]*len(cand_namelst)
    for i in range(len(cand_namelst)):
        cand_path  = "%s/%s"%(sim_path,cand_namelst[i])
        #equi_path is same as cand_path
        equi_path = cand_path
        if "sim_specialstr" in simparams.keys():
            equi_path = "%s_%s"%(equi_path,simparams["sim_specialstr"])
        cand_path = equi_path
        cand_status = cand_statuslst[i]
        if cand_status == 0:
            
            #status string
            ststr = ""

            #individual candidate descriptor dict
            descdict = dict()
            for key in data.keys():
                if "desc_" in key:
                    key2 = key.split("desc_")[1]
                    descdict[key2] = data[key][i]
            #replace values 
            candparams = simutils.replacevars(simparams,descdict)

            #run equilibration calculation
            pop = int(candparams["pop"])
            if not os.path.exists(equi_path):
                
                #simulation parameters file for every candidate
                os.makedirs(equi_path,exist_ok = True)
                cand_paramfile = "%s/params"%(equi_path)
                if "sim_specialstr" in simparams.keys():
                    cand_paramfile = "%s_%s"%(cand_paramfile,simparams["sim_specialstr"])
                cand_paramfile = "%s.yaml"%cand_paramfile
                if not os.path.exists(cand_paramfile):
                    fulldict = {'sim' : candparams}
                    with open(cand_paramfile,'w') as outfile:
                        yaml.dump(fulldict,outfile)

                #read args
                frac = float(candparams["frac"])
                ndens = float(candparams["initndens"])
                tol = candparams["tol1"]
                pause = candparams["pause"]
                clus = candparams["clus"]
                nprocs = candparams["job1_procs"]
                tolmargin = float(candparams["tolmargin1"])
                writetrj = candparams["writetrj1"]
                inpdir = "%s/inputdata"%equi_path
                ofile = "%s/in"%inpdir
                os.makedirs(inpdir,exist_ok = True)

                #make config
                simutils.makelmpconfig_packmol(pop,frac,ndens,tol,tolmargin,ofile,writetrj)

                #lmpscript,run,job
                lmp1 = lmp_scripts[0]
                shutil.copyfile("%s/%s"%(lmp_path,lmp1),"%s/%s"%(inpdir,lmp1))
                lmp1_vars = candparams["lmp1_vars"]
                if "SEED" in lmp1_vars.keys():
                    SEED = lmp1_vars["SEED"]
                    if SEED == "RANDOM":
                        lmp1_vars["SEED"] = simutils.get_randint()
                shfile = "%s/run.sh"%(inpdir)
                simutils.writerunsh(lmp1,lmp1_vars,lmp_exec,nprocs,shfile)
                jobfile = "%s/job.sh"%equi_path
                simutils.job_classic(jobfile,nprocs,clus)
                
                #submit job
                os.chdir(equi_path)
                time.sleep(pause)
                cmd = "sbatch job.sh"
                os.system(cmd)
                os.chdir(workpath)
            else:
                jobstatus_equi, outdata = simutils.trackjob(equi_path)
                if jobstatus_equi == -1:
                    cand_statuslst[i] = -1
                    ststr += "CAND_STATUS :  %s\n"%cand_statuslst[i]
                    ststr += "CAND_PATH : %s\n"%cand_pathlst[i]
                    ststr +=  "FAIL : WHILE DOING EQUILIBRATION\n"
                if jobstatus_equi == 1:
                    cand_statuslst[i] = 1
                    cand_pathlst[i] = os.path.dirname(outdata)                        
                    ststr += "CAND_STATUS :  %s\n"%cand_statuslst[i]
                    ststr += "CAND_PATH : %s\n"%cand_pathlst[i]
                    ststr +=  "SUCCESS\n"
        
        #check the candidate status again
        cand_status_after = cand_statuslst[i]
        if cand_status_after  !=  0 and cand_status == 0:
            f = open("%s/status.txt"%cand_path,'w')
            f.write(ststr)
            f.close()

    #check the status
    flag = 1   
    for i in range(len(cand_statuslst)):
        cand_status = cand_statuslst[i]
        if cand_status == 0:
            flag = -1
                        
    #if simulations finished, rewrite the batch_file, inp_file in simlogdir
    if flag == 1:
        simlogdir = os.path.abspath(simparams["simlogdir"])
        os.makedirs(simlogdir,exist_ok = True)
        bfile = "%s/batch.txt"%simlogdir
        data["cand_path"] = cand_pathlst
        data["cand_status"] = cand_statuslst
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

    #run
    run_style = params["run_style"]
    if run_style == "single":
        flag = main(params["batch_file"],params["inp_file"])
    if run_style == "recursive":
        flag = -1
        while flag == -1:
            flag = main(params["batch_file"],params["inp_file"])        
            time.sleep(30)

    print("--- %s seconds ---" % (time.time() - start_time))
