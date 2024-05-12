import os
import copy
import glob
import math
import numpy as np
import random
import time

class SimUtilsError(Exception):
    """ Raised if invalid thermodynamic ensemble is passed """

    def __init__(self, message):
        super().__init__(message)


'''
Replacing variables in indict and making outdict
indict : input dictionary
repdict : replacement dictionary
depth : search depth
'''
def replacevars(indict,repdict, depth = 1):
    outdict = copy.deepcopy(indict)
    
    #key search dict
    keyS = dict()
    for key in repdict.keys():
        keyS[key] = False

    #search at depth level 0
    for key in repdict.keys():
        count = 0
        if key in outdict.keys():
            count += 1
            outdict[key] = repdict[key]
        if count > 1:
            raise SimUtilsError("%s has been found multiple times"%key)
        if count == 1:
            keyS[key] = True
            
    #search at depth level 1
    for key1 in repdict.keys():
        count = 0
        for key2 in outdict.keys():
            if isinstance(outdict[key2], dict):
                if key1 in outdict[key2].keys():
                    count += 1
                    keyS[key1] = True
                    outdict[key2][key1] = repdict[key1]

    #check if search is succesful
    isSuccess = True
    for key,value in keyS.items():
        if value == False:
            print("%s not replaced successfully"%key)
            isSuccess = False
    if not isSuccess:
        raise SimUtilsError("All keys are not replaced.")

    return outdict


'''
Tracking the job progress based on a directory
inpdir : Input directory
datadir : data directory inside the inpdir
'''
def trackjob(inpdir,datadir = "inputdata",outname = "out.data", raiseerror = True):
    print("Tracking job progress in %s"%inpdir)    
    #initialize
    jobstatus = 0
    outdata = None
    
    isExist,slurmdir = get_slurmdir(inpdir)
    jobid = os.path.basename(slurmdir)
    if isExist:
        #job check on err and out files transfer (extra inexpensive checks (redundant at times)!!)
        errfile = "%s/slurm-%s.err"%(inpdir,os.path.basename(slurmdir))
        if not os.path.exists(errfile):
            jobstatus = 0
            print("Job finished but in the process of doing the data transfer")
            return jobstatus,outdata
        outfile = "%s/slurm-%s.out"%(inpdir,os.path.basename(slurmdir))
        if not os.path.exists(outfile):
            jobstatus = 0
            print("Job finished but in the process of doing the data transfer")
            return jobstatus,outdata
        #canceled job check errfile
        for line in reversed(open(errfile).readlines()):
            if "JOB" in line and "CANCELLED" in line and jobid in line:
                jobstatus = -1
                return jobstatus,outdata
        #finished job check outfile
        isFinish = False
        for line in reversed(open(outfile).readlines()):
            if "remove data from /tmp space" in line:
                isFinish = True
        if  isFinish :
            print("Job finished")
            outdata = "%s/%s/%s"%(slurmdir,datadir,outname)
            if os.path.exists(outdata):
                print("Success : Output data file successfully created")
                time.sleep(5)
                jobstatus = 1
                return jobstatus,outdata
            else:
                print("Failure : Output file is not created")
                jobstatus = -1
                return jobstatus,outdata                
        else:
            jobstatus = 0
            return jobstatus,outdata
            print("Job in progress")
    else:
        jobstatus = 0
        return jobstatus,outdata
        print("No slurm files exist in %s"%inpdir)
    return jobstatus,outdata


'''
Get the directory associated with slurm job id:
Directory is searched inside directory named inpdir
If more than one slurm job id directories are found, it is an ERROR
Function returns a boolean flag isExist and slurmdir     
'''
def get_slurmdir(inpdir):
    files = glob.glob("%s/slurm-*.err"%inpdir)
    slurmdir = ""
    isExist = False
    if(len(files) > 1):
        print("ERROR : More than one slurm directory exists")
        sys.exit(-1)
    if(len(files) == 1):
        isExist = True
        jobid = os.path.basename(files[0]).split("slurm-")[1].split(".")[0]
        slurmdir = "%s/%s"%(inpdir,jobid)
    return isExist,slurmdir


'''
make the packmol config in lammps
pop : population
frac : fraction (binary)
ndens : number density
tol : tolerance in packmol
tolmargin : tolerance margin used when tol is set to auto
ofile : name of output file
'''
def makelmpconfig_packmol(pop,frac,ndens,tol,tolmargin,ofile, writetrj = False):
    popA = int(round(pop*frac))
    popB = pop - popA
    box = math.ceil(math.pow((pop*1.0)/ndens,1.0/3.0))
    if tol == "auto":
        vol = box*box*box
        tol =2*tolmargin*math.pow((3.0*vol)/(math.pi*4.0*pop),1.0/3.0)
    else:
        tol = float(tol)

    #packmol
    pack_packmol(tol,box,popA,popB,ofile)

    #get coordinates
    f = open("%s.xyz"%ofile,'r')
    lst = f.readlines()
    f.close()
    ptyp = np.zeros((pop,)).astype(int)
    pc = np.zeros((pop,3))
    for i in range(2,len(lst)):
        data = lst[i].split()
        ptyp[i-2] = int(data[0])
        pc[i-2,0] = float(data[1])
        pc[i-2,1] = float(data[2])
        pc[i-2,2] = float(data[3])
    pc = pc - box/2.0

    #nbeads
    nbeads = len(pc)

    #data  str 
    datastr = ""
    datastr += "Lammps\n\n"
    datastr += "%s atoms\n\n"%nbeads
    datastr += "2 atom types\n\n"
    blow = -box/2.0
    bhigh = box/2.0
    datastr += "%s %s xlo xhi\n"%(blow,bhigh)
    datastr += "%s %s ylo yhi\n"%(blow,bhigh)
    datastr += "%s %s zlo zhi\n"%(blow,bhigh)
    datastr += "\n"
    datastr += "Masses\n\n"
    datastr += "1 1.0\n"
    datastr += "2 1.0\n"
    datastr += "\n"
    datastr += "Atoms\n\n"
    for i in range(len(pc)):
        datastr += "%d 1 %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.data"%ofile,'w')
    f.write(datastr)
    f.close()

    if not writetrj:
        return

    #trj str
    trjstr = ""
    trjstr += "ITEM: TIMESTEP\n"
    trjstr += "0\n"
    trjstr += "ITEM: NUMBER OF ATOMS\n"
    trjstr += "%s\n"%nbeads
    trjstr += "ITEM: BOX BOUNDS pp pp pp\n"
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "ITEM: ATOMS id type x y z ix iy iz\n"
    for i in range(len(pc)):
        trjstr += "%d %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.lammpstrj"%ofile,'w')
    f.write(trjstr)
    f.close()

    return



'''
make the packmol config in lammps for ternary systems
pop : population
fracA : fraction of A (binary)
fracB : fraction of B
ndens : number density
tol : tolerance in packmol
tolmargin : tolerance margin used when tol is set to auto
ofile : name of output file
'''
def makelmpconfig_packmol_ternary(pop,fracA,fracB,ndens,tol,tolmargin,ofile, writetrj = False):
    popA = int(round(pop*fracA))
    popB = int(round(pop*fracB))
    popC = pop - popA - popB
    box = math.ceil(math.pow((pop*1.0)/ndens,1.0/3.0))
    if tol == "auto":
        vol = box*box*box
        tol =2*tolmargin*math.pow((3.0*vol)/(math.pi*4.0*pop),1.0/3.0)
    else:
        tol = float(tol)

    #packmol
    pack_packmol_ternary(tol,box,popA,popB,popC,ofile)

    #get coordinates
    f = open("%s.xyz"%ofile,'r')
    lst = f.readlines()
    f.close()
    ptyp = np.zeros((pop,)).astype(int)
    pc = np.zeros((pop,3))
    for i in range(2,len(lst)):
        data = lst[i].split()
        ptyp[i-2] = int(data[0])
        pc[i-2,0] = float(data[1])
        pc[i-2,1] = float(data[2])
        pc[i-2,2] = float(data[3])
    pc = pc - box/2.0

    #nbeads
    nbeads = len(pc)

    #data  str 
    datastr = ""
    datastr += "Lammps\n\n"
    datastr += "%s atoms\n\n"%nbeads
    datastr += "3 atom types\n\n"
    blow = -box/2.0
    bhigh = box/2.0
    datastr += "%s %s xlo xhi\n"%(blow,bhigh)
    datastr += "%s %s ylo yhi\n"%(blow,bhigh)
    datastr += "%s %s zlo zhi\n"%(blow,bhigh)
    datastr += "\n"
    datastr += "Masses\n\n"
    datastr += "1 1.0\n"
    datastr += "2 1.0\n"
    datastr += "3 1.0\n"
    datastr += "\n"
    datastr += "Atoms\n\n"
    for i in range(len(pc)):
        datastr += "%d 1 %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.data"%ofile,'w')
    f.write(datastr)
    f.close()

    if not writetrj:
        return

    #trj str
    trjstr = ""
    trjstr += "ITEM: TIMESTEP\n"
    trjstr += "0\n"
    trjstr += "ITEM: NUMBER OF ATOMS\n"
    trjstr += "%s\n"%nbeads
    trjstr += "ITEM: BOX BOUNDS pp pp pp\n"
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "ITEM: ATOMS id type x y z ix iy iz\n"
    for i in range(len(pc)):
        trjstr += "%d %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.lammpstrj"%ofile,'w')
    f.write(trjstr)
    f.close()

    return


'''
make the packmol config in lammps
box : box
images : number of images
pop : population
frac : fraction (binary)
ndens : number density
tol : tolerance in packmol
tolmargin : tolerance margin used when tol is set to auto
ofile : name of output file
'''
def makelmpconfig_box_packmol(box,images,pop,frac,ndens,tol,tolmargin,ofile, writetrj = False):
    #pops
    boxI = box*images
    pop = int(round(boxI*boxI*boxI*ndens))
    popA = int(round(pop*frac))
    popB = pop - popA

    #pack
    if tol == "auto":
        vol = boxI*boxI*boxI
        tol =2*tolmargin*math.pow((3.0*vol)/(math.pi*4.0*pop),1.0/3.0)
    else:
        tol = float(tol)
    #packmol
    pack_packmol(tol,boxI,popA,popB,ofile)

    #get coordinates
    f = open("%s.xyz"%ofile,'r')
    lst = f.readlines()
    f.close()
    ptyp = np.zeros((pop,)).astype(int)
    pc = np.zeros((pop,3))
    for i in range(2,len(lst)):
        data = lst[i].split()
        ptyp[i-2] = int(data[0])
        pc[i-2,0] = float(data[1])
        pc[i-2,1] = float(data[2])
        pc[i-2,2] = float(data[3])
    pc = pc - boxI/2.0

    #nbeads
    nbeads = len(pc)

    #data  str 
    datastr = ""
    datastr += "Lammps\n\n"
    datastr += "%s atoms\n\n"%nbeads
    datastr += "2 atom types\n\n"
    blow = -boxI/2.0
    bhigh = boxI/2.0
    datastr += "%s %s xlo xhi\n"%(blow,bhigh)
    datastr += "%s %s ylo yhi\n"%(blow,bhigh)
    datastr += "%s %s zlo zhi\n"%(blow,bhigh)
    datastr += "\n"
    datastr += "Masses\n\n"
    datastr += "1 1.0\n"
    datastr += "2 1.0\n"
    datastr += "\n"
    datastr += "Atoms\n\n"
    for i in range(len(pc)):
        datastr += "%d 1 %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.data"%ofile,'w')
    f.write(datastr)
    f.close()

    if not writetrj:
        return

    #trj str
    trjstr = ""
    trjstr += "ITEM: TIMESTEP\n"
    trjstr += "0\n"
    trjstr += "ITEM: NUMBER OF ATOMS\n"
    trjstr += "%s\n"%nbeads
    trjstr += "ITEM: BOX BOUNDS pp pp pp\n"
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "ITEM: ATOMS id type x y z ix iy iz\n"
    for i in range(len(pc)):
        trjstr += "%d %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.lammpstrj"%ofile,'w')
    f.write(trjstr)
    f.close()

    return

'''
make the packmol config in lammps for ternary systems
box : box
images : number of images
pop : population
fracA : fraction of A (binary)
fracB : fraction of B (binary)
ndens : number density
tol : tolerance in packmol
tolmargin : tolerance margin used when tol is set to auto
ofile : name of output file
'''
def makelmpconfig_box_packmol_ternary(box,images,pop,fracA,fracB,ndens,tol,tolmargin,ofile, writetrj = False):
    #pops
    boxI = box*images
    pop = int(round(boxI*boxI*boxI*ndens))
    popA = int(round(pop*fracA))
    popB = int(round(pop*fracB))
    popC = pop - popA - popB

    #pack
    if tol == "auto":
        vol = boxI*boxI*boxI
        tol =2*tolmargin*math.pow((3.0*vol)/(math.pi*4.0*pop),1.0/3.0)
    else:
        tol = float(tol)
    #packmol
    pack_packmol_ternary(tol,boxI,popA,popB,popC,ofile)

    #get coordinates
    f = open("%s.xyz"%ofile,'r')
    lst = f.readlines()
    f.close()
    ptyp = np.zeros((pop,)).astype(int)
    pc = np.zeros((pop,3))
    for i in range(2,len(lst)):
        data = lst[i].split()
        ptyp[i-2] = int(data[0])
        pc[i-2,0] = float(data[1])
        pc[i-2,1] = float(data[2])
        pc[i-2,2] = float(data[3])
    pc = pc - boxI/2.0

    #nbeads
    nbeads = len(pc)

    #data  str 
    datastr = ""
    datastr += "Lammps\n\n"
    datastr += "%s atoms\n\n"%nbeads
    datastr += "3 atom types\n\n"
    blow = -boxI/2.0
    bhigh = boxI/2.0
    datastr += "%s %s xlo xhi\n"%(blow,bhigh)
    datastr += "%s %s ylo yhi\n"%(blow,bhigh)
    datastr += "%s %s zlo zhi\n"%(blow,bhigh)
    datastr += "\n"
    datastr += "Masses\n\n"
    datastr += "1 1.0\n"
    datastr += "2 1.0\n"
    datastr += "3 1.0\n"
    datastr += "\n"
    datastr += "Atoms\n\n"
    for i in range(len(pc)):
        datastr += "%d 1 %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.data"%ofile,'w')
    f.write(datastr)
    f.close()

    if not writetrj:
        return

    #trj str
    trjstr = ""
    trjstr += "ITEM: TIMESTEP\n"
    trjstr += "0\n"
    trjstr += "ITEM: NUMBER OF ATOMS\n"
    trjstr += "%s\n"%nbeads
    trjstr += "ITEM: BOX BOUNDS pp pp pp\n"
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "%s %s\n"%(blow,bhigh)
    trjstr += "ITEM: ATOMS id type x y z ix iy iz\n"
    for i in range(len(pc)):
        trjstr += "%d %s %s %s %s 0 0 0\n"%(i+1,ptyp[i],pc[i,0],pc[i,1],pc[i,2])
    count = len(pc)
    f = open("%s.lammpstrj"%ofile,'w')
    f.write(trjstr)
    f.close()

    return

'''
pack with packmol
'''
def pack_packmol(tol,boxL,pop1,pop2,ofile):
    #make tmp1.xyz
    f = open("tmp1.xyz","w")
    f.write("1\n\n")
    f.write("1\t0.1\t0.1\t0.1\n")
    f.write("end\n")
    f.close()
    #make tmp2.xyz
    f = open("tmp2.xyz","w")
    f.write("1\n\n")
    f.write("2\t0.2\t0.2\t0.2\n")
    f.write("end\n")
    f.close()

    filestr = ""
    filestr += "seed 1\n"
    filestr += "tolerance %s\n"%tol
    filestr += "filetype xyz\n"
    filestr += "output %s.xyz\n\n"%ofile

    filestr += "structure tmp1.xyz\n"
    filestr += "  number  %s\n"%pop1
    filestr += "  inside box %s %s %s %s %s %s\n"%(tol,tol,tol,boxL,boxL,boxL)
    filestr += "end structure\n\n"
    
    filestr += "structure tmp2.xyz\n"
    filestr += "  number  %s\n"%pop2
    filestr += "  inside box %s %s %s %s %s %s\n"%(tol,tol,tol,boxL,boxL,boxL)
    filestr += "end structure\n\n"
    
    f = open("packmol_input.inp",'w')
    f.write(filestr)
    f.close()

    #execute
    cmd = "packmol < packmol_input.inp"
    os.system(cmd)

    #remove tmp files
    os.remove("tmp1.xyz")
    os.remove("tmp2.xyz")



'''
pack with packmol for ternary systems
'''
def pack_packmol_ternary(tol,boxL,pop1,pop2,pop3,ofile):
    #make tmp1.xyz
    f = open("tmp1.xyz","w")
    f.write("1\n\n")
    f.write("1\t0.1\t0.1\t0.1\n")
    f.write("end\n")
    f.close()
    #make tmp2.xyz
    f = open("tmp2.xyz","w")
    f.write("1\n\n")
    f.write("2\t0.2\t0.2\t0.2\n")
    f.write("end\n")
    f.close()
    #make tmp3.xyz
    f = open("tmp3.xyz","w")
    f.write("1\n\n")
    f.write("3\t0.3\t0.3\t0.3\n")
    f.write("end\n")
    f.close()

    filestr = ""
    filestr += "seed 1\n"
    filestr += "tolerance %s\n"%tol
    filestr += "filetype xyz\n"
    filestr += "output %s.xyz\n\n"%ofile

    filestr += "structure tmp1.xyz\n"
    filestr += "  number  %s\n"%pop1
    filestr += "  inside box %s %s %s %s %s %s\n"%(tol,tol,tol,boxL,boxL,boxL)
    filestr += "end structure\n\n"
    
    filestr += "structure tmp2.xyz\n"
    filestr += "  number  %s\n"%pop2
    filestr += "  inside box %s %s %s %s %s %s\n"%(tol,tol,tol,boxL,boxL,boxL)
    filestr += "end structure\n\n"

    filestr += "structure tmp3.xyz\n"
    filestr += "  number  %s\n"%pop3
    filestr += "  inside box %s %s %s %s %s %s\n"%(tol,tol,tol,boxL,boxL,boxL)
    filestr += "end structure\n\n"
    
    f = open("packmol_input.inp",'w')
    f.write(filestr)
    f.close()

    #execute
    cmd = "packmol < packmol_input.inp"
    os.system(cmd)

    #remove tmp files
    os.remove("tmp1.xyz")
    os.remove("tmp2.xyz")
    os.remove("tmp3.xyz")

    
'''
write run.sh for lammps scripts
lmp_script : lammps script
lmp_vars :lammps variables
lmp_exec : executable
shfile : bash file to be written
nprocs : number of processors
'''
def writerunsh(lmp_script,lmp_vars,lmp_exec,nprocs,shfile):
    shstr = ""
    shstr += "#!/bin/bash\n\n"
    varstr = ""
    for key,value in lmp_vars.items():
        varstr += " -var %s %s"%(key,value)
    shstr += "mpiexec -np %s %s %s -log log.lammps -in %s > output.out\n"%(nprocs,lmp_exec,varstr,lmp_script)    
    f = open(shfile,"w")
    f.write(shstr)
    f.close()

'''
Generates an integer random number :
Return a random number in the range 1 to 32-bit integer max value
'''
def get_randint():
    maxInt = 2147483647
    return random.randint(1,maxInt)


'''
Classic writing of a job file
jobfile : name of the job file
nprocs : number of processors
clus : cluster
'''
def job_classic(jobfile,nprocs,clus, job_name = "IDP", bash_name = "run.sh" ):

    #make the bash script
    jobstring = "#!/bin/bash\n"
    jobstring += "#SBATCH  -J %s\n"%job_name
    jobstring += "#SBATCH -p %s\n"%clus
    jobstring += "#SBATCH -N 1\n"
    jobstring += "#SBATCH --exclude=c0048\n"
    jobstring += "#SBATCH --ntasks=%s\n"%nprocs
    jobstring += "#SBATCH --time=168:00:00\n"
    jobstring += "#SBATCH -o slurm-%j.out\n"
    jobstring += "#SBATCH -e slurm-%j.err\n\n"

    #echo
    jobstring += 'echo "starting $SLURM_JOBID at `date` on `hostname`"\n'
    jobstring += 'echo "my home dir is $HOME"\n\n'
    
    #copying
    jobstring += '## copying my data to a local tmp space on the compute node to reduce I/O\n'
    jobstring += 'MYTMP=/tmp/$USER/$SLURM_JOB_ID\n'
    jobstring += '/usr/bin/mkdir -p $MYTMP || exit $?\n'
    jobstring += 'echo "Copying my data over..."\n'
    jobstring += 'cp -rp $SLURM_SUBMIT_DIR/inputdata/ $MYTMP || exit $?\n\n'


    #running executable
    jobstring += "## run your job executables here...\n"
    jobstring += 'echo "Starting run jobs"\n'
    jobstring += 'cd $MYTMP/inputdata\n'
    jobstring += 'bash %s\n\n'%bash_name

    jobstring += 'echo "ended at `date` on `hostname`"\n'
    jobstring += 'echo "copy your data back to your $HOME"\n'
    jobstring += 'cp -rp $MYTMP $SLURM_SUBMIT_DIR || exit $?\n\n'

    #removing tmp data
    jobstring += '## remove your data from the compute node /tmp space\n'
    jobstring += 'echo "finish copy back"\n'
    jobstring += 'rm -rf $MYTMP\n'
    jobstring += 'echo "remove data from /tmp space"\n\n'

    jobstring += 'exit 0\n'
    
    f = open(jobfile,'w')
    f.writelines(jobstring)
    f.close()

'''
Writing of a job file in a local_py style
jobfile : name of the job file
nprocs : number of processors
clus : cluster
job_name : name of the job
py_name : name of the python code
'''
def local_py(jobfile,nprocs,clus,py_name,job_name = "SIM",):
    #make the bash script
    jobstring = "#!/bin/bash\n"
    jobstring += "#SBATCH  -J %s\n"%job_name
    jobstring += "#SBATCH -p %s\n"%clus
    jobstring += "#SBATCH -N 1\n"
    jobstring += "#SBATCH --ntasks=%s\n"%nprocs
    jobstring += "#SBATCH --time=168:00:00\n"
    jobstring += "#SBATCH -o slurm-%j.out\n"
    jobstring += "#SBATCH -e slurm-%j.err\n\n"

    #echo
    jobstring += 'echo "starting $SLURM_JOBID at `date` on `hostname`"\n'
    jobstring += 'echo "my home dir is $HOME"\n\n'

    #running
    jobstring += 'python %s\n\n'%py_name

    jobstring += 'echo "ended at `date` on `hostname`"\n'
    jobstring += 'exit 0\n'

    f = open(jobfile,'w')
    f.writelines(jobstring)
    f.close()


