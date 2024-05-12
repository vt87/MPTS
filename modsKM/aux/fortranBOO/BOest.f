      program BOest
      include 'BOmeters.h'
      character*30 RIN,AXESIN
      integer i
      integer lenstar
      integer totOrdstar
      integer totOrdavg
      double precision fracOrdstar
      double precision fracOrdavg
      double precision fracgrCstar
      double precision fracgrCavg
      include 'BOcoords.h'
      RIN = 'coord.inp'
      open (unit=1,file='coord.inp')
		read (1,*) CL(:)
		do i=1,NA
			read(1,*) r(1,i),r(2,i),r(3,i)
		enddo
      close(1)
      q(1,:)=r(1,:)/CL(1)
      q(2,:)=r(2,:)/CL(2)
      q(3,:)=r(3,:)/CL(3)
      
      MX=int(CL(1)/RCUT)
      MY=int(CL(2)/RCUT)
      MZ=int(CL(3)/RCUT)
      NCELL = MX*MY*MZ
      MAPSZE = 26*NCELL
      !write(*,*) NCELL,MAPSZE
      if (MAPSZE >= NMZ) then 
         stop 'MAPSZE is greater than NMZ in parameters file'
      endif
      if (NCELL >= NM) then 
         stop 'NCELL is greater than NM in parameters file'
      endif

      ! neighbors
      call MAPS1()
      call LINKS()
      call getallneighbor()
      call getallneighborClus()

      ! order parameter
      if (qmavgFLAG == 1) then
         call OrderParameterAVGQW()
      else 
         call OrderParameterQW()
      endif
      open(unit = 20, File = "qmavg")
      do i=1,NA
      	write(20,*)  qmavg(i)
      enddo
      close(20)

      ! dm
      lenstar = 0 
      call getqmstar(lenstar)
      Open(unit = 20, File = "qmstar")
      do i = 1,lenstar
         write(20,*) qmstar(i)
      enddo
      close(20)

      ! Ncon based on qmstar
      call getNConnectstar() 
      Open(unit = 20,File = "Nconstar")
      do i = 1,NA
         write(20,*) Nconstar(i)
      enddo
      close(20)

      ! total number of ordered particles based on epsstar
      totOrdstar = 0
      do i = 1,NA
         if (Nconstar(i) .ge. epsstar) then
            if (epsMAX .gt. 0) then
               if (Nconstar(i) .le. epsMAX) then
                  totOrdstar = totOrdstar +  1
               endif
            else
               totOrdstar = totOrdstar +  1
            endif
         endif
      enddo
      fracOrdstar = totOrdstar/real(NA)
      !write(*,*) "Total Ordered star",totOrdstar,fracOrdstar

      ! greatest cluster
      call getGreatestClusstar()
      fracgrCstar = grClusstar/real(NA)
      !write(*,*) "GreatestOrd star",grClusstar,fracgrCstar
      
      ! total number of ordered particles based on qLOW and qHIGH
      totOrdavg = 0
      do i = 1,NA
         Nqflag(i) = 0
      enddo
      do i = 1,NA
         if (qmavg(i).ge.qLOW.and.qmavg(i).le.qHIGH) then
            totOrdavg = totOrdavg +  1
            Nqflag(i) = 1
         endif
      enddo
      fracOrdavg = totOrdavg/real(NA)
      !write(*,*) "Total Ordered avg",totOrdavg,fracOrdavg

      ! greatest cluster avg
      call getGreatestClusavg()
      fracgrCavg = grClusavg/real(NA)
      !write(*,*) "GreatestOrd avg",grClusavg,fracgrCavg


      Open(unit = 20,File = "ClusStat")
      write(20,*) "N_Totstar : ",totOrdstar
      write(20,*) "f_Totstar : ",fracOrdstar
      write(20,*) "N_GRstar : ",grClusstar
      write(20,*) "f_GRstar : ",fracgrCstar
      write(20,*) "N_Totavg : ",totOrdavg
      write(20,*) "f_Totavg : ",fracOrdavg
      write(20,*) "N_GRavg : ",grClusavg
      write(20,*) "f_GRavg : ",fracgrCavg
      close(20)

      END


!-------------------------------------------------------------------------
!	SUBROUTINE OrderParameterAVGQW() - 
!       Sets the order parameter avg values for each particle
!-------------------------------------------------------------------------         
      SUBROUTINE OrderParameterAVGQW()
      include 'BOmeters.h'
      integer i,j,k,m,p
      double precision Theta,CosTheta,Phi
      double precision qij(3), rij(3)
      double precision fac3
      real*8 fac1(-Mq:Mq)
      real*8 fac2(-Mq:-1)
      complex*8 AtomMq
      complex*8 qmloc(-Mq:Mq)
      include 'BOcoords.h'

      !get prefactors
      call Prefactors(fac1,fac2,fac3)
      write(*,*) fac1,fac2,fac3

      ! qm,i
      do i = 1,NA
        do m = -Mq,Mq
          qm(m,i) = Cmplx(0.0)
        enddo

        ! partcle i neighbors
        do m = -Mq,Mq
          qmloc(m) = Cmplx(0.0)
        enddo
        do j = 1,Nb(i)
          qij(:) = q(:,i) - q(:,neighbor(j,i))
          qij(:) = qij(:) - nint(qij(:))
          rij(:)=qij(:)*CL(:)
          Theta = DATan2(Dsqrt(rij(1)**2 + rij(2)**2),rij(3))	
          if(rij(1).eq.0.0 .and. rij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(rij(2),rij(1))
          end if
          CosTheta = cos(Theta)

          qmloc(0) = qmloc(0) + cmplx(Plgndr(Mq,0,CosTheta))
          do m = 1,Mq
            AtomMq = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(Mq,m,costheta))
            qmloc(m)  = qmloc(m) + AtomMq
            qmloc(-m) = qmloc(-m)+ cmplx(fac2(-m))*conjg(AtomMq)          
          enddo      
        enddo
        do m = -Mq,Mq
          qmloc(m) = fac1(m)*qmloc(m)
          qmloc(m) = qmloc(m)/real(Nb(i))
          qm(m,i)  = qm(m,i) +  qmloc(m)
        enddo

        !neighbors of particle i 
        do j = 1,Nb(i)
           k = neighbor(j,i)
           do m = -Mq,Mq
              qmloc(m) = Cmplx(0.0)
           enddo
           do p = 1,Nb(k)
              qij(:) = q(:,k) - q(:,neighbor(p,k))
              qij(:) = qij(:) - nint(qij(:))
              rij(:)=qij(:)*CL(:)
              Theta = DATan2(Dsqrt(rij(1)**2 + rij(2)**2),rij(3))	
              if(rij(1).eq.0.0 .and. rij(2) .eq.0.0) then
                 Phi = 0.0
              else
                 Phi = DATan2(rij(2),rij(1))
              end if
              CosTheta = cos(Theta)

              qmloc(0) = qmloc(0) + cmplx(Plgndr(Mq,0,CosTheta))
              do m = 1,Mq
                 AtomMq = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                *cmplx(Plgndr(Mq,m,costheta))
                 qmloc(m)  = qmloc(m) + AtomMq
                 qmloc(-m) = qmloc(-m)+ cmplx(fac2(-m))*conjg(AtomMq)          
              enddo      
           enddo
           do m = -Mq,Mq
              qmloc(m) = fac1(m)*qmloc(m)
              qmloc(m) = qmloc(m)/real(Nb(k))
              qm(m,i)  = qm(m,i) +  qmloc(m)
           enddo
        enddo         

        !final average
        do m = -Mq,Mq
           qm(m,i) = qm(m,i)/(real(1+Nb(i)))
        enddo              
      enddo

      ! avg qm
      do i =1,NA
      	qmavg(i) = 0
      enddo	
      do i =1,NA
        do m =-Mq,Mq
          qmavg(i) = qmavg(i) + Cabs(qm(m,i))**2        
      	enddo
      	qmavg(i) = fac3*qmavg(i)
      	qmavg(i) = qmavg(i)**0.5
      enddo	      


      END



!-------------------------------------------------------------------------
!	SUBROUTINE OrderParameterQW() - 
!       Sets the order parameter values for each particle
!-------------------------------------------------------------------------         
      SUBROUTINE OrderParameterQW()
      include 'BOmeters.h'
      integer i,j,m
      double precision Theta,CosTheta,Phi
      double precision qij(3), rij(3)
      double precision fac3
      real*8 fac1(-Mq:Mq)
      real*8 fac2(-Mq:-1)
      complex*8 AtomMq
      include 'BOcoords.h'

      !get prefactors
      call Prefactors(fac1,fac2,fac3)
      write(*,*) fac1,fac2,fac3

      ! qm,i
      do i = 1,NA
        do m = -Mq,Mq
          qm(m,i) = Cmplx(0.0)
        enddo

        do j = 1,Nb(i)
          qij(:) = q(:,i) - q(:,neighbor(j,i))
          qij(:) = qij(:) - nint(qij(:))
          rij(:)=qij(:)*CL(:)
          Theta = DATan2(Dsqrt(rij(1)**2 + rij(2)**2),rij(3))	
          if(rij(1).eq.0.0 .and. rij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(rij(2),rij(1))
          end if
          CosTheta = cos(Theta)

          qm(0,i) = qm(0,i) + cmplx(Plgndr(Mq,0,CosTheta))
          do m = 1,Mq
            AtomMq = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(Mq,m,costheta))
            qm(m,i)  = qm(m,i) + AtomMq
            qm(-m,i) = qm(-m,i)+ cmplx(fac2(-m))*conjg(AtomMq)          
          enddo      
        enddo
        do m = -Mq,Mq
          qm(m,i) = fac1(m)*qm(m,i)
          qm(m,i) = qm(m,i)/real(Nb(i))
        enddo
      enddo

      ! avg qm
      do i =1,NA
      	qmavg(i) = 0
      enddo	
      do i =1,NA
        do m =-Mq,Mq
          qmavg(i) = qmavg(i) + Cabs(qm(m,i))**2        
      	enddo
      	qmavg(i) = fac3*qmavg(i)
      	qmavg(i) = qmavg(i)**0.5
      enddo	      

      END

!-------------------------------------------------------------------------
!	SUBROUTINE getqmstar(lenstar) - 
!       Gives me the cross product
!-------------------------------------------------------------------------          
      SUBROUTINE getqmstar(lenstar)
      include 'BOmeters.h'
      integer lenstar,neindex
      double precision sqiqm,sqjqm,rdpqm,idpqm
      double precision scalar      
      integer i,j,m
      include 'BOcoords.h'
      do i =1,NA
         do j = 1,NbC(i)
            neindex = neigClus(j,i)
            sqiqm =0.0
            sqjqm =0.0
            rdpqm =0.0
            idpqm =0.0
            do m =-Mq,Mq
               sqiqm = sqiqm + Cabs(qm(m,i))**2
               sqjqm = sqjqm + Cabs(qm(m,neindex))**2
               rdpqm = rdpqm + real(qm(m,i))*real(qm(m,neindex))
               idpqm = idpqm +  aimag(qm(m,i))*aimag(qm(m,neindex))
            enddo
            scalar = (rdpqm+idpqm)/((sqiqm**0.5)*(sqjqm**0.5))
            lenstar = lenstar+1
            qmstar(lenstar) = scalar
         enddo
      enddo

      END


!-------------------------------------------------------------------------
!	SUBROUTINE getNconnectstar() - 
!       Gives me the number of connections for each particle based on q
!-------------------------------------------------------------------------                
      SUBROUTINE getNConnectstar()
      include 'BOmeters.h'
      integer lenstar,neindex
      double precision sqiqm,sqjqm,rdpqm,idpqm
      double precision scalar      
      integer i,j,m
      include 'BOcoords.h'
      do  i = 1,NA
        Nconstar(i) = 0
      enddo  
      do i =1,NA
         do j = 1,NbC(i)
            neindex = neigClus(j,i)
            sqiqm =0.0
            sqjqm =0.0
            rdpqm =0.0
            idpqm =0.0
            do m =-Mq,Mq
               sqiqm = sqiqm + Cabs(qm(m,i))**2
               sqjqm = sqjqm + Cabs(qm(m,neindex))**2
               rdpqm = rdpqm + real(qm(m,i))*real(qm(m,neindex))
               idpqm = idpqm +  aimag(qm(m,i))*aimag(qm(m,neindex))
            enddo
            scalar = (rdpqm+idpqm)/((sqiqm**0.5)*(sqjqm**0.5))
            if(scalar>dcstar) then
               Nconstar(i)=Nconstar(i) +1
            endif
         enddo
      enddo      
      END


!-------------------------------------------------------------------------
!	SUBROUTINE getGreatestClusstar() - 
!       Returns the Greatest Cluster based on dcstar
!-------------------------------------------------------------------------    
      SUBROUTINE getGreatestClusstar()
      include 'BOmeters.h'
      integer i,NofCluster,hist(NA)
      include 'BOcoords.h'
      EXTERNAL harvestClusstar
      
      grClusstar =0
      NofCluster = 0
      do i=1,NA
      	belongsto(i) = -1
      	hist(i) =0
      enddo
      
      do i=1,NA
      	if(belongsto(i).eq.-1.and.Nconstar(i).ge.epsstar) then
          if (epsMAX .gt. 0) then
             if (Nconstar(i) .le. epsMAX) then
                NofCluster=NofCluster +1 
                belongsto(i)=NofCluster
                call harvestClusstar(i,NofCluster,harvestClusstar)
             endif
           else
              NofCluster=NofCluster +1 
              belongsto(i)=NofCluster
              call harvestClusstar(i,NofCluster,harvestClusstar)
          endif
        endif
      enddo 
      do i =1,NA
	if(belongsto(i) .ne. -1) then
	  hist(belongsto(i)) = hist(belongsto(i)) +1 
        endif
      enddo
      do i=1,NA 
        if(hist(i) > grClusstar) then
         grClusstar = hist(i)
        endif 
      enddo
        
      END	      
 
 
!-------------------------------------------------------------------------
!	SUBROUTINE harvestClusstar() - 
!       **Recursively harvest the Cluster Size** DANGER SUBROUTINE
!-------------------------------------------------------------------------         
      SUBROUTINE harvestClusstar(i,NofCluster,DUMHC)
      include 'BOmeters.h'
      integer i,NofCluster,c,ni
      include 'BOcoords.h'
      EXTERNAL DUMHC
      c=1
      do while (c .le. NbC(i))
      	ni = neigClus(c,i)
        if(belongsto(ni).eq.-1.and.Nconstar(ni).ge.epsstar) then
          if (epsMAX .gt. 0) then
             if (Nconstar(ni) .le. epsMAX) then
                belongsto(ni) = NofCluster
                call DUMHC(ni,NofCluster,DUMHC)
             endif
          else
             belongsto(ni) = NofCluster
             call DUMHC(ni,NofCluster,DUMHC)
          endif
        endif 
        c=c+1 
      enddo
      END


!-------------------------------------------------------------------------
!	SUBROUTINE getGreatestClusavg() - 
!     Returns the Greatest Cluster based on qmavg
!-------------------------------------------------------------------------    
      SUBROUTINE getGreatestClusavg()
      include 'BOmeters.h'
      integer i,NofCluster,hist(NA)
      include 'BOcoords.h'
      EXTERNAL harvestClusavg
      
      grClusavg =0
      NofCluster = 0
      do i=1,NA
      	belongsto(i) = -1
      	hist(i) =0
      enddo

      do i=1,NA
         if(belongsto(i).eq.-1.and.Nqflag(i).eq.1) then
      	  NofCluster=NofCluster +1 
          belongsto(i)=NofCluster
          call harvestClusavg(i,NofCluster,harvestClusavg)
        endif
      enddo 
      do i =1,NA
	if(belongsto(i) .ne. -1) then
	  hist(belongsto(i)) = hist(belongsto(i)) +1 
        endif
      enddo
      do i=1,NA 
        if(hist(i) > grClusavg) then
         grClusavg = hist(i)
        endif 
      enddo
        
      END	      
 
 
!-------------------------------------------------------------------------
!	SUBROUTINE harvestClusavg() - 
!       **Recursively harvest the Cluster Size** DANGER SUBROUTINE
!-------------------------------------------------------------------------         
      SUBROUTINE harvestClusavg(i,NofCluster,DUMHCavg)
      include 'BOmeters.h'
      integer i,NofCluster,c,ni
      include 'BOcoords.h'
      EXTERNAL DUMHCavg
      c=1
      do while (c .le. NbC(i))
      	ni = neigClus(c,i)
        if(belongsto(ni).eq.-1.and.Nqflag(ni).eq.1) then
          belongsto(ni) = NofCluster
          call DUMHCavg(ni,NofCluster,DUMHCavg)
        endif 
        c=c+1 
      enddo
      END


!-------------------------------------------------------------------------
!	SUBROUTINE getallneighbor() - 
!       Gives me the neighborlist of all the particles
!-------------------------------------------------------------------------                   
      SUBROUTINE getallneighbor()
      include 'BOmeters.h'
      double precision qij(3),rij(3)
      double precision rsq
      integer icell,i,j,append,JCELL0,NABOR,JCELL
      include 'BOcoords.h'
      do i=1,NA
        append = 1
        icell = 1 + INT (( q(1,i) + 0.5d0 ) * cellix )
     :          + INT ( ( q(2,i) + 0.5d0 ) * celliy )*MX
     :          + INT ( ( q(3,i) + 0.5d0 ) * celliz )*MX*MY
        j = HEAD(icell)
        do while(j.ne.0)
          if(j.ne.i) then
            qij(:) = q(:,i) - q(:,j)
            qij(:) = qij(:) - nint(qij(:))
            rij(:)=qij(:)*CL(:)
            rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
            if(rsq.le.rcutsqnb) then
              neighbor(append,i) = j  
              append = append+1
            endif
          endif
          j = LIST(j)
         enddo
        JCELL0 = 26*(icell - 1)
        do NABOR =1,26
          JCELL = MAPS(JCELL0 +NABOR)
          j = HEAD(JCELL)
          do while(j.ne.0)
            if(j.ne.i) then
              qij(:) = q(:,i) - q(:,j)
              qij(:) = qij(:) - nint(qij(:))
              rij(:)=qij(:)*CL(:)
              rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
              if(rsq.le.rcutsqnb) then
                neighbor(append,i) = j
                append = append+1
              endif
            endif
            j = LIST(j)
           enddo                                
        enddo
        Nb(i) = append-1
        if (Nb(i) >= NbMAX) then
           write(*,*) "Neighbors > NbMAX for particle",i
           stop 'EXITING'
        endif
      enddo
           
      END


!-------------------------------------------------------------------------
!	SUBROUTINE getallneighbor() - 
!       Gives me the neighborlist of all the particles
!-------------------------------------------------------------------------                   
      SUBROUTINE getallneighborClus()
      include 'BOmeters.h'
      double precision qij(3),rij(3)
      double precision rsq
      integer icell,i,j,append,JCELL0,NABOR,JCELL
      include 'BOcoords.h'
      do i=1,NA
        append = 1
        icell = 1 + INT (( q(1,i) + 0.5d0 ) * cellix )
     :          + INT ( ( q(2,i) + 0.5d0 ) * celliy )*MX
     :          + INT ( ( q(3,i) + 0.5d0 ) * celliz )*MX*MY
        j = HEAD(icell)
        do while(j.ne.0)
          if(j.ne.i) then
            qij(:) = q(:,i) - q(:,j)
            qij(:) = qij(:) - nint(qij(:))
            rij(:)=qij(:)*CL(:)
            rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
            if(rsq.le.rcutsqcl) then
              neigClus(append,i) = j  
              append = append+1
            endif
          endif
          j = LIST(j)
         enddo
        JCELL0 = 26*(icell - 1)
        do NABOR =1,26
          JCELL = MAPS(JCELL0 +NABOR)
          j = HEAD(JCELL)
          do while(j.ne.0)
            if(j.ne.i) then
              qij(:) = q(:,i) - q(:,j)
              qij(:) = qij(:) - nint(qij(:))
              rij(:)=qij(:)*CL(:)
              rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
              if(rsq.le.rcutsqcl) then
                neigClus(append,i) = j
                append = append+1
              endif
            endif
            j = LIST(j)
           enddo                                
        enddo
        NbC(i) = append-1
        if (NbC(i) >= NbMAX) then
           write(*,*) "Neighbors > NbMAX for particle",i
           stop 'EXITING'
        endif
      enddo
           
      END
            

!-------------------------------------------------------------------------
!	SUBROUTINE PreFactors(fac1,fac2,fac3) - 
!-------------------------------------------------------------------------  
      SUBROUTINE PreFactors(fac1,fac2,fac3)
      include 'BOmeters.h'
      integer i,j
      double precision fac3
      real*8 fac1(-Mq:Mq)
      real*8 fac2(-Mq:-1)
      real*8 fac81(-8:8)
      real*8 fac82(-8:-1)
      real*8 fac61(-6:6)
      real*8 fac62(-6:-1)
      real*8 fac41(-4:4)
      real*8 fac42(-4:-1)
      real*8 fac21(-2:2)
      real*8 fac22(-2:-1)
      include 'BOcoords.h'
           
      if (Mq == 8) then

         ! fac1
         fac81(-8) = 2.5427855325e-07
         fac81(-7) = 1.0171142130e-06
         fac81(-6) = 5.5709639801e-06
         fac81(-5) = 3.6103972995e-05
         fac81(-4) = 2.6034945177e-04
         fac81(-3) = 2.0166581818e-03
         fac81(-2) = 1.6383408518e-02
         fac81(-1) = 1.3707343005e-01
         fac81(0)  = 1.1631066229
         fac81(1)  = 1.3707343005e-01
         fac81(2)  = 1.6383408518e-02
         fac81(3)  = 2.0166581818e-03
         fac81(4)  = 2.6034945177e-04
         fac81(5)  = 3.6103972995e-05
         fac81(6)  = 5.5709639801e-06
         fac81(7)  = 1.0171142130e-06
         fac81(8)  = 2.5427855325e-07
         do i=-Mq,Mq
            fac1(i) = fac81(i)
         enddo

         ! fac2 
         fac82(-8) = 1.
         fac82(-7) = -1.
         fac82(-6) = 1.
         fac82(-5) = -1.
         fac82(-4) = 1.
         fac82(-3) = -1.
         fac82(-2) = 1.
         fac82(-1) = -1.
         do i=-Mq,-1
            fac2(i) = fac82(i)
         enddo
         
         ! fac3
         fac3 = 0.7391983

      else if (Mq == 6) then
         ! fac1
         fac61(-6) = 4.647273819e-5
         fac61(-5) = 1.609862874e-4
         fac61(-4) = 7.550926198e-4
         fac61(-3) = 0.004135812609
         fac61(-2) = 0.02481487565
         fac61(-1) = 0.15694306
         fac61(0)  = 1.0170172
         fac61(1)  = 0.15694306
         fac61(2)  = 0.02481487565
         fac61(3)  = 0.004135812609
         fac61(4)  = 7.550926198e-4
         fac61(5)  = 1.609862874e-4
         fac61(6)  = 4.647273819e-5
         do i=-Mq,Mq
            fac1(i) = fac61(i)
         enddo

         ! fac2 
         fac62(-6) = 1.
         fac62(-5) = -1.
         fac62(-4) = 1.
         fac62(-3) = -1.
         fac62(-2) = 1.
         fac62(-1) = -1.
         do i=-Mq,-1
            fac2(i) = fac62(i)
         enddo

         ! fac3
         fac3 = 0.9666439

      else if (Mq == 4) then

         ! fac1
         fac41(-4) = 0.00421459707
         fac41(-3) = 0.01192068067
         fac41(-2) = 0.04460301029
         fac41(-1) = 0.18923493915
         fac41(0)  = 0.846283753
         fac41(1)  = 0.18923493915
         fac41(2)  = 0.04460301029
         fac41(3)  = 0.01192068067
         fac41(4)  = 0.00421459707
         do i=-Mq,Mq
            fac1(i) = fac41(i)
         enddo

         ! fac2 
         fac42(-4) = 1.
         fac42(-3) = -1.
         fac42(-2) = 1.
         fac42(-1) = -1.
         do i=-Mq,-1
            fac2(i) = fac42(i)
         enddo

         ! fac3
         fac3 = 1.3962634

      else if (Mq == 2) then

         ! fac1
         fac21(-2) = 1.2875806734e-01
         fac21(-1) = 2.5751613468e-01
         fac21(0)  = 6.3078313051e-01
         fac21(1)  = 2.5751613468e-01
         fac21(2)  = 1.2875806734e-01
         do i=-Mq,Mq
            fac1(i) = fac21(i)
         enddo

         ! fac2 
         fac22(-2) = 1.
         fac22(-1) = -1.
         do i=-Mq,-1
            fac2(i) = fac22(i)
         enddo

         ! fac3
         fac3 = 2.5132741

      else
         stop 'Invalid Mq'
      endif

      END
  
!-------------------------------------------------------------------------
!	SUBROUTINE Plgndr(l,m,x) - 
!       It gives me the value of legendre polynomial
!-------------------------------------------------------------------------  
C***************************************************************************  
        FUNCTION    Plgndr( l, m, x)
        INTEGER     l, m
        REAL*8      plgndr, x
        
*   **  Computes the associated Legendre polynomial P m l (x). 
*   **  m and l are integers 0 <= m <= l, while x lies in the range - 1 <= x <= 1.
        
        INTEGER     i, ll
        REAL*8      fact, pll, pmm, pmmp1, somx2
        
        pmm=1. 
*   **  Compute Pmm .

        IF (m.gt.0) THEN
            somx2=sqrt((1.-x)*(1.+x))
            fact=1.
            DO i=1,m
                pmm=-pmm*fact*somx2
                fact=fact+2.
            ENDDO
        ENDIF
        IF (l.eq.m) THEN
            plgndr=pmm
        ELSE
            pmmp1=x*(2*m+1)*pmm 
*   **  Compute P m m+1 .
        
            IF (l.eq.m+1) THEN
                plgndr=pmmp1
            ELSE 
*   **  Compute P m l , l >m+ 1.
                DO ll=m+2,l
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/float(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                ENDDO 
                plgndr=pll
            ENDIF
        ENDIF
        
        RETURN
        END
c************       
  
  
  
  
!-------------------------------------------------------------------------
!	SUBROUTINE MAPS1() - It builds the neighborList for each cell
!-------------------------------------------------------------------------	
	SUBROUTINE MAPS1()

	  include 'BOmeters.h'
        INTEGER     IX, IY, IZ, IMAP, ICELL
		
	  include 'BOcoords.h'
C    *******************************************************************

C    ** STATEMENT FUNCTION TO GIVE CELL INDEX **

        ICELL ( IX, IY, IZ) = 1 + MOD ( IX - 1 + MX, MX )
     :                          + MOD ( IY - 1 + MY, MY) * MX
     :                          + MOD ( IZ - 1 + MZ, MZ ) 
     :	* MX * MY

C    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **

        DO 50 IZ = 1, MZ

           DO 40 IY = 1, MY

              DO 30 IX = 1, MX

                 IMAP = ( ICELL ( IX, IY, IZ ) - 1 ) * 26

                 MAPS( IMAP + 1  ) = ICELL( IX + 1, IY    , IZ     )
                 MAPS( IMAP + 2  ) = ICELL( IX + 1, IY + 1, IZ     )
                 MAPS( IMAP + 3  ) = ICELL( IX    , IY + 1, IZ     )
                 MAPS( IMAP + 4  ) = ICELL( IX - 1, IY + 1, IZ     )
                 MAPS( IMAP + 5  ) = ICELL( IX + 1, IY    , IZ - 1 )
                 MAPS( IMAP + 6  ) = ICELL( IX + 1, IY + 1, IZ - 1 )
                 MAPS( IMAP + 7  ) = ICELL( IX    , IY + 1, IZ - 1 )
                 MAPS( IMAP + 8  ) = ICELL( IX - 1, IY + 1, IZ - 1 )
                 MAPS( IMAP + 9  ) = ICELL( IX + 1, IY    , IZ + 1 )
                 MAPS( IMAP + 10 ) = ICELL( IX + 1, IY + 1, IZ + 1 )
                 MAPS( IMAP + 11 ) = ICELL( IX    , IY + 1, IZ + 1 )
                 MAPS( IMAP + 12 ) = ICELL( IX - 1, IY + 1, IZ + 1 )
                 MAPS( IMAP + 13 ) = ICELL( IX    , IY    , IZ + 1 )
		 MAPS( IMAP + 14 ) = ICELL( IX - 1, IY    , IZ     )
                 MAPS( IMAP + 15 ) = ICELL( IX - 1, IY - 1, IZ     )
                 MAPS( IMAP + 16 ) = ICELL( IX    , IY - 1, IZ     )
                 MAPS( IMAP + 17 ) = ICELL( IX + 1, IY - 1, IZ     )
                 MAPS( IMAP + 18 ) = ICELL( IX - 1, IY    , IZ + 1 )
                 MAPS( IMAP + 19 ) = ICELL( IX - 1, IY - 1, IZ + 1 )
		 MAPS( IMAP + 20 ) = ICELL( IX    , IY - 1, IZ + 1 )
                 MAPS( IMAP + 21 ) = ICELL( IX + 1, IY - 1, IZ + 1 )
                 MAPS( IMAP + 22 ) = ICELL( IX - 1, IY    , IZ - 1 )
                 MAPS( IMAP + 23 ) = ICELL( IX - 1, IY - 1, IZ - 1 )
                 MAPS( IMAP + 24 ) = ICELL( IX    , IY - 1, IZ - 1 )
                 MAPS( IMAP + 25 ) = ICELL( IX + 1, IY - 1, IZ - 1 )
                 MAPS( IMAP + 26 ) = ICELL( IX    , IY    , IZ - 1 )

30            CONTINUE

40         CONTINUE

50      CONTINUE

        RETURN
        END
	
	
!-------------------------------------------------------------------------
!	SUBROUTINE LINKS() - Bulid the Head and List Array
!-------------------------------------------------------------------------	
	
        SUBROUTINE LINKS ()
	  include 'BOmeters.h'
        INTEGER     ICELL, I
	  include 'BOcoords.h'
        DO 10 ICELL = 1, NCELL
           HEAD(ICELL) = 0
10      CONTINUE

        CELLIX = REAL ( MX )
        CELLIY = REAL ( MY )
        CELLIZ = REAL ( MZ )

C    ** SORT ALL ATOMS **

        DO 20 I = 1, NA

           ICELL = 1 + INT ( ( q(1,I) + 0.5d0 ) * CELLIX )
     :               + INT ( ( q(2,I) + 0.5d0 ) * CELLIY ) * MX
     :               + INT ( ( q(3,I) + 0.5d0 ) * CELLIZ ) *MY*MX
		 if (q(1,I).eq.0.5d0) ICELL=ICELL-1
		 if (q(2,I).eq.0.5d0) ICELL=ICELL-MX
		 if (q(3,I).eq.0.5d0) ICELL=ICELL-MX*MY

 
           LIST(I)     = HEAD(ICELL)
           HEAD(ICELL) = I

20      CONTINUE

        RETURN
        END
        
        
