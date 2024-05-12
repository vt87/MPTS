      implicit real*8(a-h, o-z)
      integer MX,MY,MZ,NA,NbMAX,LIST
      integer NM,NMZ,MAPSZE,NCELL,MAPS,HEAD
      integer neighbor,Nb
      integer neigClus,NbC
      integer Nconstar,epsstar,belongsto
      integer epsMAX
      integer Mq
      integer Nqflag
      integer grClusstar,grClusavg
      double precision dcstar
      double precision qmstar
      complex*8 qm
      double precision qmavg
      double precision qLOW, qHIGH
      double precision r,q,CL
      double precision rcutsq,rcut
      double precision rcutsqnb,rcutnb
      double precision rcutsqcl,rcutcl
      double precision cellix,celliy,celliz
      integer qmavgFLAG
      PARAMETER(NA=3456)
      PARAMETER(rcutnb=3.0d0)
      PARAMETER(rcutcl=3.0d0)
      PARAMETER(rcut=4.0d0)
      PARAMETER(epsstar=8)
      PARAMETER(epsMAX=18)
      PARAMETER(dcstar=0.75d0)
      PARAMETER(NbMAX=500)
      PARAMETER(NM=50000)
      PARAMETER(NMZ=1300000)
      PARAMETER(rcutsqnb=rcutnb**2.0d0)
      PARAMETER(rcutsq=rcut**2.0d0)
      PARAMETER(rcutsqcl=rcutcl**2.0d0)
      PARAMETER(Mq=6)
      PARAMETER(qmavgFLAG=0)
      PARAMETER(qLOW=0.05d0)
      PARAMETER(qHIGH=0.1d0)


