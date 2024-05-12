        common /coords/ r(3,NA),q(3,NA)
	common /boxsize/ CL(3),cellix,celliy,celliz,MX,MY,MZ
	common /var/ LIST(NA)
	common /var2/ neighbor(NbMAX,NA) ,Nb(NA)
	common /var7/ neigClus(NbMAX,NA) ,NbC(NA)
	common /var3/ belongsto(NA) 
	common /var4/ Nconstar(NA)
        common /var5/ grClusstar,grClusavg
	common /var6/ Nqflag(NA)
	common /map/ MAPS(NMZ),HEAD(NM),NCELL,MAPSZE
	common /Qmval/ qm(-Mq:Mq,NA)
	common /Qmavgval/ qmavg(NA)
	common /Qmstar/ qmstar(NA*NA)
