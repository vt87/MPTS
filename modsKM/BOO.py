#-----------------------------------------------------------------------------------------------
# Snippets toolkit.
# Vikram Thapar, vt87@cornell.edu, Cornell University
# Helper functions for parsing the files using python
# File Name: BOO.py
# Estimation of Bond Order parameters in conjuction with MDAnalysis tools
#--------------------------------------------------------------------------------------

import MDAnalysis as mda
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import copy

#Global plot setting
plt.rc('font', family='serif')
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


#modify grmeters
def modifygrmeters(fortranPATH,grvars = dict()):
    f = open("%s/grmetersFIX.h"%fortranPATH,'r')
    lst = f.readlines()
    f.close()

    #modify NA
    subs  = "PARAMETER(NA="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(NA=%s)\n"%grvars["NA"]
    else:
        raise Exception("NA parameter not found in fortran parameters file")

    #modify grdr
    subs  = "PARAMETER(grdr="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(grdr=%sd0)\n"%grvars["grdr"]
    else:
        raise Exception("grdr parameter not found in fortran parameters file")

    #modify grbins
    subs  = "PARAMETER(grbins="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(grbins=%s)\n"%grvars["grbins"]
    else:
        raise Exception("grbins parameter not found in fortran parameters file")

    #write
    f = open("%s/grmeters.h"%fortranPATH,'w')
    f.writelines(lst)
    f.close()


#modify BOmeters
def modifyBOmeters(fortranPATH,BOvars = dict()):
    f = open("%s/BOmetersFIX.h"%fortranPATH,'r')
    lst = f.readlines()
    f.close()

    #modify NA
    subs  = "PARAMETER(NA="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(NA=%s)\n"%BOvars["NA"]
    else:
        raise Exception("NA parameter not found in fortran parameters file")

    #modify rcutnb
    subs  = "PARAMETER(rcutnb="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(rcutnb=%sd0)\n"%BOvars["rcutnb"]
    else:
        raise Exception("rcutnb parameter not found in fortran parameters file")

    #modify rcut
    subs  = "PARAMETER(rcut="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(rcut=%sd0)\n"%BOvars["rcut"]
    else:
        raise Exception("rcut parameter not found in fortran parameters file")

    #modify rcutcl
    subs  = "PARAMETER(rcutcl="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(rcutcl=%sd0)\n"%BOvars["rcutcl"]
    else:
        raise Exception("rcutcl parameter not found in fortran parameters file")

    #modify Mq 
    subs  = "PARAMETER(Mq="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(Mq=%s)\n"%BOvars["Mq"]
    else:
        raise Exception("Mq parameter not found in fortran parameters file")

    #modify epsstar
    subs  = "PARAMETER(epsstar="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(epsstar=%s)\n"%BOvars["epsstar"]
    else:
        raise Exception("epsstar parameter not found in fortran parameters file")

    #modify epsMAX
    subs  = "PARAMETER(epsMAX="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(epsMAX=%s)\n"%BOvars["epsMAX"]
    else:
        raise Exception("epsMAX parameter not found in fortran parameters file")

    #modify qmavgFLAG 
    subs  = "PARAMETER(qmavgFLAG="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(qmavgFLAG=%s)\n"%BOvars["qmavgFLAG"]
    else:
        raise Exception("qmavgFLAG parameter not found in fortran parameters file")

    #modify qLOW
    subs  = "PARAMETER(qLOW="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(qLOW=%sd0)\n"%BOvars["qLOW"]
    else:
        raise Exception("qLOW parameter not found in fortran parameters file")

    #modify qHIGH
    subs  = "PARAMETER(qHIGH="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(qHIGH=%sd0)\n"%BOvars["qHIGH"]
    else:
        raise Exception("qHIGH parameter not found in fortran parameters file")

    #modify dcstar
    subs  = "PARAMETER(dcstar="
    res = [i for i,s in enumerate(lst) if subs in s]
    if len(res) == 1:
        lst[res[0]] = "      PARAMETER(dcstar=%sd0)\n"%BOvars["dcstar"]
    else:
        raise Exception("dcstar parameter not found in fortran parameters file")

    #write
    f = open("%s/BOmeters.h"%fortranPATH,'w')
    f.writelines(lst)
    f.close()


#read Cluster file
def readclus(filename,clusdict):
    f = open(filename,'r')
    lst = f.readlines()
    f.close()
    
    #take deep copy of clusdict
    clusdict = copy.deepcopy(clusdict)

    #read
    string = ["N_Totstar","f_Totstar","N_GRstar","f_GRstar","N_Totavg","f_Totavg","N_GRavg","f_GRavg"]
    #if "TotOrdStar" in clusdict.keys():
    for k in range(len(string)):
        res = [i for i,s in enumerate(lst) if string[k] in s]
        if len(res) != 1:
            raise Exception("%s parameter not found or found multiple times in %s file"%(string,filename))
        val = lst[res[0]].split(":")[1].rstrip()
        if "N_" in string[k]:
            val = int(val)
        else:
            val = float(val)
        if string[k] in clusdict.keys():
            clusdict[string[k]].append(val)
        else:
            clusdict[string[k]] = [val]
            
    #return
    return string,clusdict
            

#Class definition
class BOO:
    def __init__(self,topfile,trrfile):
        #Class variables
        self.topfile = topfile  #Topology file
        self.trrfile = trrfile  #Trajectory File

        #Error check
        self.errorcheck()

        #clusdict
        self.clusdict = dict()
        
    #Error check
    def errorcheck(self):
        
        #Files error check
        if not os.path.exists(self.topfile):
            raise Exception("Topology file does not exist")
        if not os.path.exists(self.trrfile):
            raise Exception("Trajectory file does not exist")


    '''
    Radial distribution function
    '''
    def GRest(self,vardict,plotdict = dict()):

        #copy vars in vardict
        Nprod = vardict["Nprod"]
        outfile = vardict["outfile"]
        if "selatoms" in vardict.keys():
            selatoms = vardict["selatoms"]

        #Nstep
        if "Nstep" in vardict.keys():
            Nstep = vardict["Nstep"]
        else:
            Nstep = 1

        #interval size 
        if "dr" in vardict.keys():
            dr = vardict["dr"]
        else:
            dr = 0.05

        #MDAnalysis universe variable
        u = mda.Universe(self.topfile,self.trrfile)
        trajlen = len(u.trajectory)

        #Nprod reset
        if(trajlen <= Nprod):
            print("Resetting Nprod from %d to %d"%(Nprod,trajlen)) 
            Nprod = trajlen
            N1 = 0
            N2 = Nprod
        else:
            N1 = trajlen - Nprod
            N2 = trajlen
        #replace N1 and N2 if they are present in vardict.keys()
        if "N1" in vardict.keys():
            N1 = vardict["N1"]
            print("Rewriting N1 to  %d"%(N1)) 
        if "N2" in vardict.keys():
            N2 = vardict["N2"]
            if N2 >= trajlen:
                N2 = trajlen
            print("Rewriting N2 to  %d"%(N2)) 
        Nframes = int((N2 - N1-1)/float(Nstep)) + 1
        
        #trajectory looping for box length
        boxlst = []
        for ts in u.trajectory[N1:N2]:
            time = u.trajectory.time
            boxdim = np.copy(u.dimensions)[0:3]
            boxlst.append(boxdim)
        boxarr = np.asarray(boxlst)
        avgbox = np.mean(boxarr,axis = 0)
        stdbox = np.std(boxarr,axis = 0)

        #maxr 
        if "r_max" in vardict.keys():
            r_max = vardict["r_max"]
        else:
            r_max = avgbox[0]/2.0

        #arange
        nbins = math.floor(r_max/dr)
        #fortranPATH
        if "fortranPATH" in vardict.keys():
            fortranPATH = vardict["fortranPATH"]
            isfortran = True
        else:
            isfortran = False

        #selection of atoms and size
        if "selatoms" in vardict.keys():
            atoms = u.atoms.select_atoms(selatoms)
            selstr = selatoms
        else:
            selstr = "all"
            atoms =  u.atoms
        size = len(atoms)    

        #isfortran
        if isfortran:
            workpath = os.getcwd()
            #make grvars
            grvars = dict()
            grvars["grdr"] = dr
            grvars["grbins"] = nbins
            grvars["NA"] = size 

            #modify grmeters.h
            modifygrmeters(fortranPATH,grvars)

            #change directory, compile 
            os.chdir(fortranPATH)
            cmd = "gfortran gr3D.f -o gr3D_execute"
            os.system(cmd)
            os.chdir(workpath)
        else:
            raise Exception("OP evaluation can only be done under fortran framework")


        #trajectory looping
        frame = -1
        for ts in u.trajectory[N1:N2:Nstep]:
            time = u.trajectory.time
            frame += 1

            print("Gr estimation for frame number %s"%frame)
            
            #positions
            boxdim = np.expand_dims(np.copy(u.dimensions)[0:3],axis = 0)
            r = np.copy(atoms.positions)

            #print("Correcting small precision errors to run in fortran")
            r += boxdim/2
            r = r - boxdim*np.floor(r/boxdim)
            r -= boxdim/2
            #print(np.min(r,axis = 0)/boxdim)
            #print(np.max(r,axis = 0)/boxdim)

            #concatenate array
            inpc = np.vstack((boxdim,r))

            #get gr through fortran
            os.chdir(fortranPATH)
            #write positions 
            np.savetxt('coord.inp',inpc,fmt = "%s")
            #run gr code in fortran
            cmd = "./gr3D_execute"
            os.system(cmd)
            #read gr.txt
            arr = np.loadtxt("gr.txt")
            os.chdir(workpath)

            #update avgarr
            if frame == 0:
                avgarr = np.copy(arr)
            else:
                avgarr[:,1] += arr[:,1]

        #average
        avgarr[:,1] /= Nframes

        #save the results
        f = open(outfile,'w')
        f.write("# Gr estimation : Courtesy of MDAnalysis, FortranOP and OPboo\n")
        f.write("# Inputs\n")
        f.write("# dr : %s\n"%dr)
        f.write("# nbins : %s\n"%nbins)
        f.write("# selstr (selatoms in input dict) : %s\n"%(selstr))        
        f.write("# Outputs\n")
        f.write("# Frames Read :  %s to %s\n"%(N1,N2))        
        f.write("# Average Lattice Parameter : %s %s %s\n"%(avgbox[0],avgbox[1],avgbox[2]))
        f.write("# Standard deviation in Lattice Parameter : %s %s %s\n"%(stdbox[0],stdbox[1],stdbox[2]))
        f.write("# r gr\n")
        for i in range(nbins):
            f.write("%s\t%s\n"%(avgarr[i,0],avgarr[i,1]))
        f.close()


        #plotting
        if len(plotdict) == 0:
            return

        imname = plotdict["imname"]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("Radial distribution function")    
        ax1.set_xlabel("r")
        ax1.set_ylabel("g(r)")
        ax1.plot(avgarr[:,0],avgarr[:,1])
        plt.savefig(imname)
        plt.close()


    '''
    Bond orientational order parameters
    '''
    def BOOest(self,vardict,plotdict = dict()):

        #copy vars in vardict
        Nprod = vardict["Nprod"]
        if "selatoms" in vardict.keys():
            selatoms = vardict["selatoms"]

        #avg file check
        if "outavgfile" in vardict.keys():
            isOutqavg = True
            outavgfile = vardict["outavgfile"]
        else:
            isOutqavg = False

        #clus file check
        if "outclusfile" in vardict.keys():
            isOutclus = True
            outclusfile = vardict["outclusfile"]
        else:
            isOutclus = False

        #star file check
        if "outstarfile" in vardict.keys():
            isOutqstar = True
            outstarfile = vardict["outstarfile"]
        else:
            isOutqstar = False

        #ncon file check
        if "outnconfile" in vardict.keys():
            isOutncon = True
            outnconfile = vardict["outnconfile"]
        else:
            isOutncon = False

        #Nstep
        if "Nstep" in vardict.keys():
            Nstep = vardict["Nstep"]
        else:
            Nstep = 1

        #MDAnalysis universe variable
        u = mda.Universe(self.topfile,self.trrfile)
        trajlen = len(u.trajectory)

        #Nprod reset
        if(trajlen <= Nprod):
            print("Resetting Nprod from %d to %d"%(Nprod,trajlen)) 
            Nprod = trajlen
            N1 = 0
            N2 = Nprod
        else:
            N1 = trajlen - Nprod
            N2 = trajlen
        #replace N1 and N2 if they are present in vardict.keys()
        if "N1" in vardict.keys():
            N1 = vardict["N1"]
            print("Rewriting N1 to  %d"%(N1)) 
        if "N2" in vardict.keys():
            N2 = vardict["N2"]
            if N2 >= trajlen:
                N2 = trajlen
            print("Rewriting N2 to  %d"%(N2)) 
        Nframes = int((N2 - N1-1)/float(Nstep)) + 1

        #fortranPATH
        if "fortranPATH" in vardict.keys():
            fortranPATH = vardict["fortranPATH"]
            isfortran = True
        else:
            isfortran = False

        #selection of atoms and size
        if "selatoms" in vardict.keys():
            atoms = u.atoms.select_atoms(selatoms)
            selstr = selatoms
        else:
            selstr = "all"
            atoms =  u.atoms
        size = len(atoms)    

        
        #variables for BOmeters
        rcutnb = vardict["rcutnb"]
        Mq = vardict["Mq"]
        if "qLOW" in vardict.keys():
            qLOW = vardict["qLOW"]
        else:
            qLOW = 0.05
        if "qHIGH" in vardict.keys():
            qHIGH = vardict["qHIGH"]
        else:
            qHIGH = 0.1
        if "dcstar" in vardict.keys():
            dcstar = vardict["dcstar"]
        else:
            dcstar = 0.7
        if "epsstar" in vardict.keys():
            epsstar = vardict["epsstar"]
        else:
            epsstar = 7
        if "rcut" in vardict.keys():
            rcut = vardict["rcut"]
        else:
            rcut = rcutnb + 1.0
        if "rcutcl" in vardict.keys():
            rcutcl = vardict["rcutcl"]
        else:
            rcutcl = rcutnb
        if "qmavgFLAG" in vardict.keys():
            qmavgFLAG = vardict["qmavgFLAG"]
        else:
            qmavgFLAG = 0
        if "epsMAX" in vardict.keys():
            epsMAX = vardict["epsMAX"]
        else:
            epsMAX = -1

        #isfortran
        if isfortran:
            #compiler flag
            if "compilestr" in vardict.keys():
                compilestr = vardict["compilestr"]
            else:
                compliestr = "gfortran"

            workpath = os.getcwd()
            #make BOvars
            BOvars = dict()
            BOvars["NA"] = size
            BOvars["rcutnb"] = rcutnb
            BOvars["rcut"] = rcut
            BOvars["rcutcl"] = rcutcl
            BOvars["Mq"] = Mq
            BOvars["qmavgFLAG"] = qmavgFLAG
            BOvars["qLOW"] = qLOW
            BOvars["qHIGH"] = qHIGH
            BOvars["dcstar"] = dcstar
            BOvars["epsstar"] = epsstar
            BOvars["epsMAX"] = epsMAX

            #modify BOmeters.h
            modifyBOmeters(fortranPATH,BOvars)

            #change directory, compile 
            os.chdir(fortranPATH)
            if compilestr == "gfortran":
                cmd = "gfortran BOest.f -o BOest_execute"
            else:
                cmd = "ifort -O3 BOtest.f -o BOest_execute"
            os.system(cmd)
            os.chdir(workpath)
        else:
            raise Exception("OP evaluation can only be done under fortran framework")


        #trajectory looping
        frame = -1
        frame_id = N1
        frameids = []
        for ts in u.trajectory[N1:N2:Nstep]:
            time = u.trajectory.time
            frame += 1

            print("BOO estimation for frame number %s"%frame_id)

            #update frame_id
            frameids.append(frame_id)
            frame_id += Nstep
            
            #positions
            boxdim = np.expand_dims(np.copy(u.dimensions)[0:3],axis = 0)
            r = np.copy(atoms.positions)

            #print("Correcting small precision errors to run in fortran")
            r += boxdim/2
            r = r - boxdim*np.floor(r/boxdim)
            r -= boxdim/2
            #print(np.min(r,axis = 0)/boxdim)
            #print(np.max(r,axis = 0)/boxdim)

            #concatenate array
            inpc = np.vstack((boxdim,r))

            #get BOO through fortran
            os.chdir(fortranPATH)
            #write positions 
            np.savetxt('coord.inp',inpc,fmt = "%s")
            #remove ClusStat first
            try:
                os.remove("ClusStat")
            except OSError:
                pass
            #run BOO code in fortran
            cmd = "./BOest_execute"
            os.system(cmd)

            #check ClusStat"
            if not os.path.exists("ClusStat"):
                os.chdir(workpath)
                print(workpath)
                return "FAILED"

            #read qmavg,qmstar,ClusStat
            if isOutqavg:
                arr = np.loadtxt("qmavg")
            if isOutqstar:
                arrstar = np.loadtxt("qmstar")
            if isOutncon:
                arrncon = np.loadtxt("Nconstar")
            if isOutclus:
                clusstr,self.clusdict = readclus("ClusStat",self.clusdict)
            os.chdir(workpath)

            #update avgarr
            if isOutqavg:
                if frame == 0:
                    avgarr = np.copy(arr)
                else:
                    avgarr = np.concatenate([avgarr,arr])

            #update qmstar if required
            if isOutqstar:
                if frame == 0:
                    avgarrstar = np.copy(arrstar)
                else:
                    avgarrstar = np.concatenate([avgarrstar,arrstar])

            #update ncon
            if isOutncon:
                if frame == 0:
                    avgarrncon = np.copy(arrncon)
                else:
                    avgarrncon = np.concatenate([avgarrncon,arrncon])
            

        #optional plotting , qavg
        if isOutqavg:                    
            #histogram parameters
            if "nbins" in vardict.keys():
                nbins = vardict["nbins"]
            else:
                nbins = 50
            if "qmin" in vardict.keys():
                qmin = vardict["qmin"]
            else:
                qmin = np.min(avgarr)
            if "qmax" in vardict.keys():
                qmax = vardict["qmax"]
            else:
                qmax = np.max(avgarr)
            #bins
            bins = np.linspace(qmin,qmax,nbins+1)
        
            #histogram
            qhist,edges = np.histogram(avgarr,bins,density = True)
            print(edges,nbins,qmin,qmax,bins)

            centers = 0.5*(edges[1:]+ edges[:-1])
        
            #save the avg results
            f = open(outavgfile,'w')
            f.write("# BOO avg estimation : Courtesy of MDAnalysis, FortranOP and OPboo\n")
            f.write("# Inputs\n")
            f.write("# rcutnb : %s\n"%rcutnb)
            f.write("# rcut : %s\n"%rcut)
            f.write("# nbins : %s\n"%nbins)
            f.write("# qmin : %s\n"%qmin)
            f.write("# qmax : %s\n"%qmax)
            f.write("# Mq : %s\n"%Mq)
            f.write("# selstr (selatoms in input dict) : %s\n"%(selstr))        
            f.write("# Outputs\n")
            f.write("# Frames Read :  %s to %s\n"%(N1,N2))        
            f.write("# centers qmavg\n")
            for i in range(len(centers)):
                f.write("%s\t%s\n"%(centers[i],qhist[i]))
            f.close()

            #plotting
            if "imname" in plotdict.keys():
                imname = plotdict["imname"]
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.set_title("Average Q%s"%Mq)    
                ax1.set_xlabel("q")
                ax1.set_ylabel("P(q)")
                ax1.plot(centers,qhist)
                plt.savefig(imname)
                plt.close()
        

        #optional plotting , qstar
        if isOutqstar:
            #histogram parameters
            if "nbinsstar" in vardict.keys():
                nbinsstar = vardict["nbinsstar"]
            else:
                nbinsstar = 50
            if "qstarmin" in vardict.keys():
                qstarmin = vardict["qstarmin"]
            else:
                qstarmin = np.min(avgarrstar)
            if "qstarmax" in vardict.keys():
                qstarmax = vardict["qstarmax"]
            else:
                qstarmax = np.max(avgarrstar)
            #bins
            binsstar = np.linspace(qstarmin,qstarmax,nbinsstar+1)
            #histogram
            qstarhist,edges = np.histogram(avgarrstar,binsstar,density = True)
            print(edges,nbinsstar,qstarmin,qstarmax,binsstar)

            centers = 0.5*(edges[1:]+ edges[:-1])
            
            #save the q star results
            f = open(outstarfile,'w')
            f.write("# BOO star estimation : Courtesy of MDAnalysis, FortranOP and OPboo\n")
            f.write("# Inputs\n")
            f.write("# rcutnb : %s\n"%rcutnb)
            f.write("# rcut : %s\n"%rcut)
            f.write("# nbinsstar : %s\n"%nbinsstar)
            f.write("# qstarmin : %s\n"%qstarmin)
            f.write("# qstarmax : %s\n"%qstarmax)
            f.write("# Mq : %s\n"%Mq)
            f.write("# selstr (selatoms in input dict) : %s\n"%(selstr))        
            f.write("# Outputs\n")
            f.write("# Frames Read :  %s to %s\n"%(N1,N2))        
            f.write("# centers qmstar\n")
            for i in range(len(centers)):
                f.write("%s\t%s\n"%(centers[i],qstarhist[i]))
            f.close()

            if "imnamestar" in plotdict.keys():
                imnamestar = plotdict["imnamestar"]
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.set_title("Q%sstar"%Mq)    
                ax1.set_xlabel("qstar")
                ax1.set_ylabel("P(qstar)")
                ax1.plot(centers,qstarhist)
                plt.savefig(imnamestar)
                plt.close()


        #optional plotting , ncon
        if isOutncon:
            #histogram parameters
            if "nconmin" in vardict.keys():
                nconmin = vardict["nconmin"]
            else:
                nconmin = np.min(avgarrncon)
            if "nconmax" in vardict.keys():
                nconmax = vardict["nconmax"]
            else:
                nconmax = np.max(avgarrncon)
            #bins
            binsncon = np.arange(nconmin,nconmax+2).astype(int)
            #histogram
            nconhist,edges = np.histogram(avgarrncon.astype(int),binsncon,density = True)

            #save the q star results
            f = open(outnconfile,'w')
            f.write("# BOO ncon estimation : Courtesy of MDAnalysis, FortranOP and OPboo\n")
            f.write("# Inputs\n")
            f.write("# rcutnb : %s\n"%rcutnb)
            f.write("# rcut : %s\n"%rcut)
            f.write("# dcstar : %s\n"%dcstar)
            f.write("# nconmin : %s\n"%nconmin)
            f.write("# nconmax : %s\n"%nconmax)
            f.write("# Mq : %s\n"%Mq)
            f.write("# selstr (selatoms in input dict) : %s\n"%(selstr))        
            f.write("# Outputs\n")
            f.write("# Frames Read :  %s to %s\n"%(N1,N2))        
            f.write("# edges qmstar\n")
            for i in range(len(nconhist)):
                f.write("%s\t%s\n"%(edges[i],nconhist[i]))
            f.close()

            if "imnamencon" in plotdict.keys():
                imnamencon = plotdict["imnamencon"]
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.set_title("Ncon%sstar"%Mq)    
                ax1.set_xlabel("Ncon")
                ax1.set_ylabel("P(Ncon)")
                ax1.plot(edges[:-1],nconhist)
                plt.savefig(imnamencon)
                plt.close()
            
        

        #optional plotting , clusfile
        if isOutclus:

            #get average and standard deviation of cluster
            for i in range(len(clusstr)):
                key = clusstr[i]
                value = self.clusdict[key]
                self.clusdict["%s_avg"%key] = np.mean(np.asarray(value))
                self.clusdict["%s_std"%key] = np.std(np.asarray(value))
                
            #save the clus results
            f = open(outclusfile,'w')
            f.write("# BOO cluster size estimation : Courtesy of MDAnalysis, FortranOP and OPboo\n")
            f.write("# Inputs\n")
            f.write("# rcutnb : %s\n"%rcutnb)
            f.write("# rcut : %s\n"%rcut)
            f.write("# qLOW : %s\n"%qLOW)
            f.write("# qHIGH : %s\n"%qHIGH)
            f.write("# dcstar : %s\n"%dcstar)
            f.write("# epsstar : %s\n"%epsstar)
            f.write("# Mq : %s\n"%Mq)
            f.write("# selstr (selatoms in input dict) : %s\n"%(selstr))        
            f.write("# Outputs\n")
            f.write("# Frames Read :  %s to %s\n"%(N1,N2))
            #averages
            for i in range(len(clusstr)):
                key = "%s_avg"%clusstr[i]
                value = self.clusdict[key]
                f.write("# %s : %s\n"%(key,value))
            #standard deviations
            for i in range(len(clusstr)):
                key = "%s_std"%clusstr[i]
                value = self.clusdict[key]
                f.write("# %s : %s\n"%(key,value))
            #string clus
            strclus = ""
            for i in range(len(clusstr)):
                strclus += " %s"%clusstr[i]
            f.write("# frame_id%s\n"%strclus)
            
            #write values
            for i in range(len(frameids)):
                valstr = "%s"%frameids[i]
                for k in range(len(clusstr)):
                    val = self.clusdict[clusstr[k]][i]
                    valstr += " %s"%val
                f.write("%s\n"%valstr)
            f.close()
            
