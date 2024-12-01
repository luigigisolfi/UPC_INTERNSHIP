       SUBROUTINE TPX(N,VI,IPO,IND,IFI,PES,FU,CO,FNO,CNO,TOL,CAMP)
C********************************************************************
C  PARALLEL SHOOTING MINIMUM NORM CORRECTION FORWARD OR BACKWARD
C  IN TIME, SOLVING THE OVERDET. SYSTEM (DF,-I)CO=-F BY TXOLESKY.
C
C  NOTE: THIS ROUTINE DEALS WITH THE SYSTEM OF EQUATIONS RELATED 
C        TO THE NEWTON PROCEDURE IN THE FORM (DF,-I)CO=-F
C        INSTEAD OF (-I,DF)CO=-F AS IN OTHER ROUTINES OF
C        PARALLEL SHOOTING.
C
C  INPUT VALUES:
C  
C   N    NUMBER OF NODES OF THE PARALLEL SHOOTING (ALL OF THEM)
C   VI   VI(7*N) CONTAINING SEQUENTIALLY THE INITIAL VALUES OF
C        TIME, POS AND VEL FOR EACH ONE OF THE NODES. IT MUST
C        HAPPEN THAT IF I.GT.J THEN THE TIME VALUE OF NODE I .GT.
C        TIME VALUE OF NODE J.
C   IND  +1 OR -1 ACCORDING WHETHER PARALLEL SHOOTING MUST BE
C        PERFORMED FORWARD (from node 1 to 2, from 2 to 3, ...
C        from n-1 to n), OR BACKWARD (from node 2 to 1, from 3 
C        to 2,...,from n to n-1) IN TIME.
C   IPO  NUMBER INDICATING THE TYPE OF THE PROPAGATION IMPLEMENTED
C        IN ROUTINE PROTPX AND USED IN THE FUNCTION WE WANT TO
C        MAKE ZERO. BUT AT THIS MOMENT ONLY IPO=1 IS ALLOWED,
C        SINCE THE EPOCHS OF THE NODES ARE NOT CHANGING.
C        (SEE THIS ROUTINE BELOW FOR MORE DETAILS)
C   IFI  IFI(7). IN THIS VERSION OF THE PARALLEL SHOOTING THIS
C        VECTOR IS NOT USED.
C   PES  PES(7*N) WEIGHT VECTOR FOR THE LEAST SQUARES PROCEDURE.
C   TOL  THRESHOLD ERROR VALUE FOR THE INTEGRATOR RK78.
C   CAMP NAME OF THE ROUTINE CONTAINING THE VECTORFIELD. THE
C        VECTORFIELD IS ASSUMED GIVEN AS IN ROUTINE DERIVP
C        FOR EXAMPLE.
C
C  OUTPUT VALUES:
C  
C   FU   FU(7*(N-1)) THE VALUES OF THE FUNCTION DEFINING THE
C        PARALLEL SHOOTING WHICH WE WANT TO MAKE ZERO.
C        NOTE: WHEN IND.GT.0 THE VALUES OF FU ARE TIME, POS AND
C              VEL AT NODES 2,3,...,N-1,N. 
C              WHEN IND.LT.0 THE VALUES OF FU ARE TIME, POS AND
C              VEL AT NODES N-1,N-2,N-3,...,1.
C   CO   CO(7*N) THE VECTOR CONTAINING THE CORRECTIONS. THE
C        NEW GUESS MUST BE COMPUTED OUTSIDE THIS ROUTINE DOING
C        VI=VI+CO. THEN YOU CAN CALL AGAIN THIS ROUTINE FOR THE
C        NEXT ITERATION. WE REMEMBER THAT THE CORRECTIONS OF
C        THE EPOCH COMPONENTS ARE ALWAYS ZERO.
C  FNO   THE 2-NORM OF THE FUNCTION. (NO WEIGHTS)
C  CNO   THE 2-NORM OF THE CORRECTION (COMPUTED WITH OR WITHOUT 
C        WEIGHTS ACCORDING TO PARAMETER IPES IN ROUTINE ESCITX)
C
C NOTES: -TROUGH CHANNEL NCR AND IN FILE ARXIT (SEE PARAMETERS)
C         THE ROUTINE WRITES INFORMATION ABOUT THE PARALLEL SHOOT.
C        -DECOMENTING THE LINE "CALL TESSOL(..)" YOU HAVE INFORMATION
C         ABOUT THE OBTAINED CORRECTION.
C********************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*64 ARXIT
       PARAMETER (NMA=30000,NCR=8,ARXIT='tpx.itr')
       DIMENSION DF(7,7*(NMA-1)),FU(7*(N-1)),CO(7*N),VI(7*N)
       DIMENSION PES(7*N),X(48),IFI(7),AA(6,6),AAT(6,6)
       DIMENSION A(6,6*NMA-6),B(6*NMA-6),DX(6*NMA),W(6*NMA)
       DIMENSION EL(6,6*NMA-6),D(6,6*NMA-6),U(6*NMA-6),V(6*NMA-6)
       DIMENSION Y(6*NMA-6)
       EXTERNAL CAMP
       NBM=NMA-1
       IF (IND.NE.1.AND.IND.NE.-1) STOP
       IF (IPO.NE.1) THEN
       WRITE (*,*) 'TPX. ONLY THE VALUE IPO=1 IS IMPLEMENTED'
       WRITE (*,*) '     YOU GAVE IPO=',IPO
       STOP
       ENDIF
       IF (N.GT.NMA) THEN
       WRITE (*,*) 'TPX. VALUE OF N TOO LARGE: ',N
       WRITE (*,*) '     CHANGE PARAMETER NMA: ',NMA
       STOP
       ENDIF
C PUTTING THE WEIGHTS IN THE RIGHT ORDER IF NEEDED
       IF (IND.LT.0) THEN
       DO 2 I=1,N/2
       L1=7*(N-I)
       L2=7*(I-1)
       DO 3 J=1,7
       AUX=PES(L1+J)
       PES(L1+J)=PES(L2+J)
       PES(L2+J)=AUX
3      CONTINUE
2      CONTINUE
       ENDIF
C STARTING MAIN PART OF THE PROGRAM
       II=1
       IF=N-1
       IF (IND.LT.0) THEN
       II=2
       IF=N
       ENDIF
C INTEGRATION OF EACH INITIAL CONDITION
       DO 10 NI=II,IF
       DO 11 I=7,48
       X(I)=0.D0
11     CONTINUE
       DO 12 I=1,6
       X(I)=VI(7*NI-6+I)
       X(7*I+6)=1.D0
12     CONTINUE
       T=VI(7*NI-6)
       CALL PROTPX(IPO,T,X,VI,NI,IND,TOL,CAMP)
C STORING FUNCTION AND DIFFERENTIAL
C La matriu es guarda diferent que en altres metodes de tir 
C on es feia   NFIL=7*(N-NI)-6 i la FU esta canv de signe
       ITF=7*NI+1
       NFIL=7*(NI-1)+1
       IF (IND.LT.0) THEN
       ITF=7*NI-13
       NFIL=7*(N-NI)+1
       ENDIF
       FU(NFIL)=VI(ITF)-T
       DF(1,NFIL)=1.D0
       DO 15 I=1,6
       FU(NFIL+I)=VI(ITF+I)-X(I)
       DF(I+1,NFIL)=0.D0
       DF(1,NFIL+I)=X(6+I)
       NC=6*(I+1)
       DO 16 J=1,6
       DF(I+1,NFIL+J)=X(NC+J)
16     CONTINUE
15     CONTINUE
10     CONTINUE
C POSEM LES COSES EN MATRIUS 6*6 APTES PER A RESTPX
       NB=N-1
       AUX=1.D0
       DO 17 J=1,NB
       L1=7*(J-1)
       L2=6*(J-1)-1
       DO 19 K=2,7
       DO 21 I=1,6
       A(I,L2+K)=DF(I+1,L1+K)
21     CONTINUE
19     CONTINUE
       CALL GETM(NB,A,J,AA)
       CALL TRAN(AA,AAT)
       CALL PUTM(NB,AAT,J,AUX,A)
17     CONTINUE
       DO 25 J=1,NB
       L1=7*(J-1)
       L2=6*(J-1)-1
       DO 27 K=2,7
       B(L2+K)=FU(L1+K)
       W(L2+K)=PES(L1+K)
27     CONTINUE
25     CONTINUE
       DO 29 K=2,7
       W(6*(N-1)+K-1)=PES(7*(N-1)+K)
29     CONTINUE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c aixo pot servir per escriure les matrius en binari
c       write (*,*) 'COMENCO A ESCRIURE ....'
c       open(4,file='df.bin',form='unformatted')
c       write(4) N
c       do 30 i=1,7*(n-1)
c       write(4) fu(i)
c       do 31 j=1,7
c       write(4) df(j,i)
c31     continue
c30     continue
c       close(4)
c       write (*,*) 'JA L''HE ESCRIT !!!!'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C COMPUTING THE SOLUTION OF (DF,-I)CO=-F. FIRST OF ALL THE
C SOLUTION, DX, CORRESPONDING TO THE 6*6 PARTS IS COMPUTED
C AND THEN STORED IN THE FINAL ONE, CO.
       write (*,*) 'tpx. passo a calcular la solucio del sistema'
       CALL RESTPX(NB,A,B,DX,W,NBM,EL,D,U,V,Y)
       DO 50 J=1,N
       L1=7*(J-1)
       L2=6*(J-1)-1
       CO(L1+1)=0.D0
       DO 51 K=2,7
       CO(L1+K)=DX(L2+K)
51     CONTINUE
50     CONTINUE
C TESTINGS OF THE SOLUTIONS IF WANTED.
C      CALL TESSOL(DF,FU,CO,NMA,N,6)
C SETTING THE CORRECTION AND WEIGHTS IN THE INITIAL ORDER 
C WHEN IND.LT.0
       IF (IND.LT.0) THEN
       DO 102 I=1,N/2
       L1=7*(N-I)
       L2=7*(I-1)
       DO 103 J=1,7
       AUX=CO(L1+J)
       CO(L1+J)=CO(L2+J)
       CO(L2+J)=AUX
       AUX=PES(L1+J)
       PES(L1+J)=PES(L2+J)
       PES(L2+J)=AUX
103    CONTINUE
102    CONTINUE
       ENDIF
C WRITTING INFORMATION ABOUT THE COMPUTATIONS, IT COMPUTES
C AS WELL THE NORMS OF THE FUNCTION AND OF THE CORRECTION.
c+++++++++++++++++++++++++++++++++++++++++++
c aquesta rutina escriu un arxiu on la correccio es veu en
c coordenades adimensionals. Ara pero no es pot fer servir ja
c que correspon a altres verions de tir paral.lel i cal 
c modificar-la una mica (arxiu sdfr.f)
c       CALL CORRN(VI,PA,CS,CO,NMA,N)
c+++++++++++++++++++++++++++++++++++++++++++
       CALL ESCITX(N,FU,CO,PES,NMA,FNO,CNO,IND,NCR,ARXIT)
       RETURN
       END



       SUBROUTINE ESCITX(N,FU,CO,PES,NMA,FNO,CNO,IND,NCR,ARXIT)
C********************************************************************
C WRITTING INFORMATION OF THE PARALLEL SHOOTING THROUGH CHANNEL
C NCR IN FILE ARXIT AND THROUGH THE STD OUTPUT.
C PARAMETER IPES 1,0 CONTROLS WHETHER THE NORM OF THE CORRECTION 
C MUST BE COMPUTED WITH (1) OR WITHOUT (0) WEIGHTS.
C********************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (IPES=0)
       CHARACTER*64 ARXIT
       DIMENSION FU(7*(N-1)),CO(7*N),XMAC(3),XMIC(3),XNC(3),V(3)
       DIMENSION XMAF(3),XMIF(3),XNF(3),PES(7*N)
       OPEN (NCR,FILE=ARXIT)
       WRITE (NCR,*) '*****  INFORMATION FROM ROUTINE TPX  *****'
       WRITE (NCR,*) 'N and IND: ',N,IND
       DO 1 I=1,3
       XNF(I)=0.D0
       XNC(I)=0.D0
       XMAF(I)=0.D0
       XMAC(I)=0.D0
       XMIF(I)=1.D64
       XMIC(I)=1.D64
1      CONTINUE
       FNO=0.D0
       CNO=0.D0
       DO 5 I=1,N
       L=7*(I-1)
       IF (IPES.EQ.1) THEN
       V(1)=CO(L+1)*CO(L+1)*PES(L+1)
       V(2)=CO(L+2)*CO(L+2)*PES(L+2)+CO(L+3)*CO(L+3)*PES(L+3)+
     $      CO(L+4)*CO(L+4)*PES(L+4)
       V(3)=CO(L+5)*CO(L+5)*PES(L+5)+CO(L+6)*CO(L+6)*PES(L+6)+
     $      CO(L+7)*CO(L+7)*PES(L+7)
       ELSE
       V(1)=CO(L+1)*CO(L+1)
       V(2)=CO(L+2)*CO(L+2)+CO(L+3)*CO(L+3)+CO(L+4)*CO(L+4)
       V(3)=CO(L+5)*CO(L+5)+CO(L+6)*CO(L+6)+CO(L+7)*CO(L+7)
       ENDIF
       DO 6 K=1,3
       IF (V(K).GT.XMAC(K)) XMAC(K)=V(K)
       IF (V(K).LT.XMIC(K)) XMIC(K)=V(K)
       XNC(K)=XNC(K)+V(K)
       CNO=CNO+V(K)
6      CONTINUE
       IF (I.EQ.N) GOTO 5
       V(1)=FU(L+1)*FU(L+1)
       V(2)=FU(L+2)*FU(L+2)+FU(L+3)*FU(L+3)+FU(L+4)*FU(L+4)
       V(3)=FU(L+5)*FU(L+5)+FU(L+6)*FU(L+6)+FU(L+7)*FU(L+7)
       DO 7 K=1,3
       IF (V(K).GT.XMAF(K)) XMAF(K)=V(K)
       IF (V(K).LT.XMIF(K)) XMIF(K)=V(K)
       XNF(K)=XNF(K)+V(K)
       FNO=FNO+V(K)
7      CONTINUE
5      CONTINUE
       DO 8 I=1,3
       XNF(I)=DSQRT(XNF(I))
       XNC(I)=DSQRT(XNC(I))
       XMAF(I)=DSQRT(XMAF(I))
       XMAC(I)=DSQRT(XMAC(I))
       XMIF(I)=DSQRT(XMIF(I))
       XMIC(I)=DSQRT(XMIC(I))
8      CONTINUE
       FNO=DSQRT(FNO)
       CNO=DSQRT(CNO)
       WRITE (*,*) 'ROUTINE ESCITX'
       WRITE (*,*) 'NORM OF THE FUNCTION:   ',FNO
       WRITE (*,*) 'NORM in TIME, POS and VEL, MAX and MIN VALUES:'
       WRITE (*,110) (XNF(I),I=1,3)
       WRITE (*,110) (XMAF(I),I=1,3)
       WRITE (*,110) (XMIF(I),I=1,3)
       WRITE (NCR,*) 'NORM OF THE FUNCTION:   ',FNO
       WRITE (NCR,*) 'NORM in TIME, POS and VEL, MAX and MIN VALUES:'
       WRITE (NCR,110) (XNF(I),I=1,3)
       WRITE (NCR,110) (XMAF(I),I=1,3)
       WRITE (NCR,110) (XMIF(I),I=1,3)
       WRITE (*,*) 'PARAMETER IPES=', IPES
       WRITE (NCR,*) 'PARAMETER IPES=', IPES
       IF (IPES.EQ.1) THEN
       WRITE (*,*) ' CORRECTION NORM COMPUTED   WITH  WEIGHTS'
       WRITE (NCR,*) ' CORRECTION NORM COMPUTED   WITH  WEIGHTS'
       ELSE
       WRITE (*,*) ' CORRECTION NORM COMPUTED   WITHOUT  WEIGHTS'
       WRITE (NCR,*) ' CORRECTION NORM COMPUTED   WITHOUT  WEIGHTS'
       ENDIF
       WRITE (*,*) 'NORM OF THE CORRECTION: ',CNO
       WRITE (*,*) 'NORM in TIME, POS and VEL, MAX and MIN VALUES:'
       WRITE (*,110) (XNC(I),I=1,3)
       WRITE (*,110) (XMAC(I),I=1,3)
       WRITE (*,110) (XMIC(I),I=1,3)
       WRITE (NCR,*) 'NORM OF THE CORRECTION: ',CNO
       WRITE (NCR,*) 'NORM in TIME, POS and VEL, MAX and MIN VALUES:'
       WRITE (NCR,110) (XNC(I),I=1,3)
       WRITE (NCR,110) (XMAC(I),I=1,3)
       WRITE (NCR,110) (XMIC(I),I=1,3)
110    FORMAT(1X,3E18.7)
C WRITTING THE FUNCTION
       WRITE (NCR,*) 'VALUES OF THE FUNCTION FU AT THE NODE:'
       DO 10 L=1,N-1
       IN=L+1
       IF (IND.LT.0) IN=N-L
       WRITE (NCR,100) IN,(FU(7*(L-1)+K),K=1,7)
10     CONTINUE
C WRITTING THE CORRECTION
       WRITE (NCR,*) 'VALUES OF THE CORRECTION CO AT THE NODE:'
       DO 14 L=1,N
       WRITE (NCR,100) L,(CO(7*(L-1)+K),K=1,7)
14     CONTINUE
100    FORMAT (I4,7E12.4)
       CLOSE(NCR)
       RETURN
       END


       SUBROUTINE TESSOL(DF,FU,CO,NMA,N,NCR)
C*******************************************************************
C THIS ROUTINE DOES A TEST OF THE CORRECTION (-I,DF)CO=-F. 
C*******************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION DF(7,7*(NMA-1)),FU(7*(N-1)),CO(7*NMA)
       WRITE (NCR,*) 'TESTING CORRECTION AS A SOL. OF (-I,DF)CO=-F: '
       VMX=0.D0
       ERR=0.D0
       DO 20 I=1,N-1
       LL=7*(I-1)
       DO 23 J=1,7
       VAL=-FU(LL+J)-CO(LL+J+7)
       IF (DABS(CO(LL+J)).GT.VMX) THEN
       VMX=DABS(CO(LL+J))
       IMAX=LL+J
       ENDIF
       DO 25 K=1,7
       VAL=VAL+DF(K,LL+J)*CO(LL+K)
25     CONTINUE
       IF (DABS(VAL).GT.ERR) THEN
       ERR=DABS(VAL)
       IEN=LL+J
       ENDIF
23     CONTINUE
       WRITE (4,100) I,(CO(LL+K),K=1,7)
20     CONTINUE
       WRITE (4,100) N,(CO(7*(N-1)+K),K=1,7)
100    FORMAT(I5,7E24.16)
       WRITE (NCR,*) '*** MAX ERROR IN CORREC.: ',IEN,ERR
       DO 27 I=7*(N-1)+1,7*N
       IF (DABS(CO(I)).GT.VMX) THEN
       VMX=DABS(CO(I))
       IMAX=I
       ENDIF
27     CONTINUE
       WRITE (NCR,*) 'MAX VALUE OF VECTOR: ',IMAX,VMX
       WRITE (NCR,*)
       CLOSE(4)
       RETURN
       END



       SUBROUTINE PROTPX(IPO,T,X,VI,NI,IND,TOL,CAMP)
C***************************************************************
C PROPAGATION OF THE INITIAL CONDITION.
C IPO INDICATES THE TYPE OF PROPAGATION. THE POSSIBILITIES
C ARE:
C  IPO=1  THE INITIAL CONDITION IS INTEGRATED (FORWARD OR
C         BACKWARDS ACCORDING TO IND) UP TO THE EPOCH OF THE 
C         NEXT TARGET POINT IN THE PARALLEL SHOOTING.
C  IPO=2  SAME AS BEFORE BUT WHEN REACHING THE TARGET EPOCH THE
C         INTEGRATION IS FOLLOWED AND STOPPED WHEN THE NEAREST
C         LOCAL MINIMUM OF THE DISTANCE TO THE NEXT TARGET POINT
C         IN THE PARALLEL SHOOTING IS FOUND. THIS SECOND
C         INTEGRATION IS PERFORMED FORWARD OR BACKWARDS IN TIME
C         INDEPENDENTLY OF IND AND JUST TRYING TO FIND THE MIN.
C***************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NEV=48,IPOM=2,EFM=1.D-3,EDTM=1.D-10)
       DIMENSION X(*),VI(*),R(13,NEV),B(NEV),F(NEV)
       EXTERNAL CAMP
       HMIN=1.D-4
       HMAX=1.D0
       H=IND*5.D-2
       IF (IPO.EQ.2) GOTO 100
       IF (IPO.GT.2) THEN
          WRITE (*,*) 'PROTPX. ONLY ',IPOM,' TYPES'
          WRITE (*,*) '         OF COMPUTATIONS IPO: ',IPO
          STOP
       ENDIF
C TYPE 1. FROM EPOCH TO EPOCH.
       ITF=7*NI+1
       IF (IND.LT.0) ITF=7*NI-13
       TF=VI(ITF)
10     CALL RK78 (T,X,NEV,H,HMIN,HMAX,TOL,R,B,F,CAMP)
       IF ((TF-T)*IND.GT.0.D0) GOTO 10
       H=TF-T
       CALL RK78 (T,X,NEV,H,HMIN,HMAX,TOL,R,B,F,CAMP)
       RETURN
C TYPE 2. FROM EPOCH TO NEAREST DISTANCE NEAR EPOCH.
100    ITF=7*NI+1
       IF (IND.LT.0) ITF=7*NI-13
       IXF=ITF+1
       TI=T
       TF=VI(ITF)
101    CALL RK78 (T,X,NEV,H,HMIN,HMAX,TOL,R,B,F,CAMP)
       IF ((TF-T)*IND.GT.0.D0) GOTO 101
       D0=(X(1)-VI(IXF))**2+(X(2)-VI(IXF+1))**2+(X(3)-VI(IXF+2))**2
       ITER=0
110    IF (ITER.GT.5) WRITE (*,*) 'PROTPX. NI, DISTANCE: ',NI,DSQRT(D0)
       ITER=ITER+1
       FF=(X(1)-VI(IXF))*X(4)+(X(2)-VI(IXF+1))*X(5)+
     $    (X(3)-VI(IXF+2))*X(6)
       CALL CAMP(T,X,NEV,F)
       DFF=(X(1)-VI(IXF))*F(4)+(X(2)-VI(IXF+1))*F(5)+
     $     (X(3)-VI(IXF+2))*F(6)+X(4)*X(4)+X(5)*X(5)+X(6)*X(6)
       H=-FF/DFF
       IF (ITER.GT.5) WRITE (*,*)' ITER,DT,FUN and DER: ',ITER,H,FF,DFF
       TAF=T+H
105    CALL RK78 (T,X,NEV,H,HMIN,HMAX,TOL,R,B,F,CAMP)
       IF (DABS(TAF-T).GT.HMIN) GOTO 105
       IF (DABS(TAF-T).GT.1.D-10) THEN
       H=TAF-T
       CALL RK78 (T,X,NEV,H,HMIN,HMAX,TOL,R,B,F,CAMP)
       ENDIF
       D0=(X(1)-VI(IXF))**2+(X(2)-VI(IXF+1))**2+(X(3)-VI(IXF+2))**2
       IF (DABS(FF).GT.EFM.AND.DABS(FF/DFF).GT.EDTM) GOTO 110
       IF (ITER.GT.5)WRITE (*,*)'PROTPX. NI, LAST DIST: ',NI,DSQRT(D0)
       IF ((T-TI)*IND.LT.0.D0) THEN
       WRITE (*,*) 'PROTPX. WARNING ! MINIMUM BACKWARDS:'
       WRITE (*,*) 'IND, TI, TMIN: ',IND,TI,T
       ENDIF
       RETURN
       END



      subroutine restpx (nb,a,b,x,w,nbm,el,d,u,v,y)
      implicit real*8(a-h,o-z)
      dimension a(6,6*nbm),w(6*(nbm+1)),b(6*nbm),x(6*(nbm+1))
      dimension el(6,nbm*6),d(6,nbm*6),u(nbm*6),v(nbm*6),y(nbm*6)
      dimension xmax(6),nmax(6)
      if(nb.gt.nbm) then
	 write(*,*) ' A restpx: nb =',nb,'es mes gran que nbm =',nbm
	 stop
      endif
      call solveit(nb,a,w,b,x,nbm,el,d,u,v,y)
      nblo=0
      do i=1,6
	 xmax(i)=0
      enddo
      do 9 i=1,6*(nb+1)
	 ii=mod(i-1,6)+1
	 if(dabs(x(i)).gt.xmax(ii))then
	    xmax(ii)=dabs(x(i))
	    nmax(ii)=nblo+1
         endif
9     continue
10    format(i6,d24.16)
      write(*,*)' correccions maximes i blocs'
      do i=1,6
	 write(*,*) nmax(i),xmax(i)
      enddo
c     call check(nb,a,w,b,x,nbm)
      close(2)
      return 
      end

      subroutine solveit(nb,a,w,b,x,nbm,el,d,u,v,y)
      implicit real*8(a-h,o-z)
      dimension a(6,6*nbm),w(6*(nbm+1)),b(nbm*6),x(6*(nbm+1)),
     1   el(6,nbm*6),d(6,nbm*6),u(nbm*6),v(nbm*6),y(nbm*6),
     2   aa(6,6),wwi(6,6),aat(6,6),dd(6,6),ell(6,6),ellt(6,6),
     3   uu(6),vv(6),yy(6),xx(6),bb(6),
     4   em1(6,6),em2(6,6),em3(6,6)
      call getm(nb,a,1,aa)
      call getiw(nb+1,w,1,wwi)
      call tran(aa,aat)
      call prom(aa,wwi,em1)
      call prom(em1,aat,em3)
      call getiw(nb+1,w,2,wwi)
      call summ(em3,wwi,1.d0,dd)
      call putm(nb,dd,1,1.d0,d)
      do 1 n=2,nb
         call matin(dd,em1,isin)
	 if(isin.eq.1)then
	    write(*,*)' el error es per n=',n,' en crida a matin'
	    stop
         endif
         call getm(nb,a,n,aa)
         call getiw(nb+1,w,n,wwi)
         call prom(aa,wwi,em2)
         call prom(em2,em1,ell)
         call putm(nb,ell,n,-1.d0,el)
         call tran(aa,aat)
         call prom(em2,aat,em3)
         call getiw(nb+1,w,n+1,wwi)
         call summ(em3,wwi,1.d0,em1)
         call prom(ell,dd,em2)
         call tran(ell,ellt)
         call prom(em2,ellt,em3)
         call summ(em1,em3,-1.d0,dd)
         call putm(nb,dd,n,1.d0,d)
1     continue
      if (n.lt.100) then
C Basicament evitem trams massa curts. De fet pero no es necessari.
      write (*,*) 'solveit. Aqui diu que pari si n<100. n=',n
      stop
      endif
      call getv(nb,b,1,bb)
      call sumv(bb,bb,0.d0,uu)
      call putv(nb,uu,1,1.d0,u)
      do 2 n=2,nb
         call getv(nb,b,n,vv)
         call getm(nb,el,n,ell)
         call promv(ell,uu,yy)
         call sumv(vv,yy,-1.d0,uu)
         call putv(nb,uu,n,1.d0,u)
2     continue
      do 3 n=1,nb
         call getm(nb,d,n,dd)
         call getv(nb,u,n,uu)
         call solv(dd,uu,vv,isin)
	 if(isin.eq.1)then
	    write(*,*)' el error es per n=',n,' en crida a solv'
	    stop
         endif
         call putv(nb,vv,n,1.d0,v)
3     continue
      call putv(nb,vv,nb,1.d0,y)
      do 4 n=nb-1,1,-1
         call getv(nb,v,n,vv)
         call getm(nb,el,n+1,ell)
         call getv(nb,y,n+1,uu)
         call tran(ell,ellt)
         call promv(ellt,uu,xx)
         call sumv(vv,xx,-1.d0,yy)
         call putv(nb,yy,n,1.d0,y)
4     continue
      call getm(nb,a,1,aa)
      call tran(aa,aat)
      call getv(nb,y,1,yy)
      call promv(aat,yy,xx)
      call putv(nb+1,xx,1,1.d0,x)
      do 5 n=2,nb
         call getm(nb,a,n,aa)
         call tran(aa,aat)
         call getv(nb,y,n,yy)
         call promv(aat,yy,uu)
         call getv(nb,y,n-1,yy)
         call sumv(uu,yy,-1.d0,xx)
         call putv(nb+1,xx,n,1.d0,x)
5     continue
      call getv(nb,y,nb,yy)
      call putv(nb+1,yy,nb+1,-1.d0,x)
      do 6 i=1,6*(nb+1)
6     x(i)=x(i)/w(i)
      return
      end

      subroutine check(nb,a,w,b,x,nbm)
      implicit real*8(a-h,o-z)
      dimension a(6,6*nbm),w(6*(nbm+1)),b(nbm*6),x(6*(nbm+1)),
     1   aa(6,6),bb(6),xx(6),yy(6),zz(6)
      do 1 n=1,nb
	 call getm(nb,a,n,aa)
         call getv(nb+1,x,n,xx)
	 call promv(aa,xx,yy)
	 call getv(nb+1,x,n+1,xx)
         call sumv(yy,xx,-1.d0,zz) 
         call getv(nb,b,n,bb)
	 call sumv(zz,bb,-1.d0,yy)
	 call norm(yy,vno)
         if(vno.gt.1.d-12)write(*,*)' n,vno=',n,vno
1     continue
      sq=0
      do 2 i=1,6*(nb+1)
2     sq=sq+w(i)*x(i)*x(i)
      write(*,*)' suma de quadrats',sq
      return
      end

      subroutine norm(v,vno)
      implicit real*8(a-h,o-z)
      dimension v(6)
      vno=0
      do 1 i=1,6
1     vno=vno+v(i)*v(i)
      vno=dsqrt(vno)
      return
      end

      subroutine getm(nb,a,n,aa)
      implicit real*8(a-h,o-z)
      dimension a(6,6*nb),aa(6,6)
      n6=n*6-6
      do 1 j=1,6
         jj=j+n6
         do 2 i=1,6
2        aa(i,j)=a(i,jj)
1     continue
      return
      end

      subroutine getiw(nbm1,w,n,ww)
      implicit real*8(a-h,o-z)
      dimension w(6*nbm1),ww(6,6)
      do 1 j=1,6
         do 2 i=1,6
2        ww(i,j)=0
1     continue
      n6=n*6-6
      do 3 j=1,6
3     ww(j,j)=1/w(j+n6)
      return
      end

      subroutine prom(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(6,6),b(6,6),c(6,6)
      do 1 j=1,6
         do 2 i=1,6
            c(i,j)=0
            do 3 k=1,6
3           c(i,j)=c(i,j)+a(i,k)*b(k,j)
2        continue
1     continue
      return
      end

      subroutine tran(a,at)
      implicit real*8(a-h,o-z)
      dimension a(6,6),at(6,6)
      do 1 j=1,6
         do 2 i=1,6
2        at(j,i)=a(i,j)
1     continue
      return
      end

      subroutine putm(nb,aa,n,coef,a)
      implicit real*8(a-h,o-z)
      dimension a(6,6*nb),aa(6,6)
      n6=n*6-6
      do 1 j=1,6
         do 2 i=1,6
2        a(i,j+n6)=coef*aa(i,j)
1     continue
      return
      end

      subroutine matin(a,ain,isin)
      implicit real*8(a-h,o-z)
      dimension a(6,6),ac(6,6),ain(6,6),b(6),ipivot(6)
      call copy(a,ac)
      call matinv(ac,ain,b,6,ipivot,isin)
      return
      end

      subroutine copy(a,ac)
      implicit real*8(a-h,o-z)
      dimension a(6,6),ac(6,6)
      do 1 j=1,6
         do 2 i=1,6
2        ac(i,j)=a(i,j)
1     continue
      return
      end

      subroutine summ(a,b,coef,c)
      implicit real*8(a-h,o-z)
      dimension a(6,6),b(6,6),c(6,6)
      do 1 j=1,6
         do 2 i=1,6
2        c(i,j)=a(i,j)+coef*b(i,j)
1     continue
      return
      end

      subroutine getv(nbm1,b,n,bb)
      implicit real*8(a-h,o-z)
      dimension b(6*nbm1),bb(6)
      ii=6*n-6
      do 1 i=1,6
1     bb(i)=b(i+ii)
      return
      end

      subroutine sumv(a,b,coef,c)
      implicit real*8(a-h,o-z)
      dimension a(6),b(6),c(6)
      do 1 i=1,6
1     c(i)=a(i)+coef*b(i)
      return
      end

      subroutine putv(nbm1,uu,n,coef,u)
      implicit real*8(a-h,o-z)
      dimension u(6*nbm1),uu(6)
      ii=6*n-6
      do 1 i=1,6
1     u(i+ii)=coef*uu(i)
      return
      end

      subroutine promv(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(6,6),b(6),c(6)
      do 1 i=1,6
         c(i)=0
         do 2 j=1,6
2        c(i)=c(i)+a(i,j)*b(j)
1     continue
      return
      end

      subroutine solv(a,b,x,isin)
      implicit real*8(a-h,o-z)
      dimension a(6,6),ac(6,6),b(6),x(6),d(6),ipivot(6)
      call copy(a,ac)
      call solin(ac,b,x,d,ipivot,6,isin)
      return
      end

      subroutine matinv(a,ainv,b,n,ipivot,isin)
      implicit real*8(a-h,o-z)
      dimension a(n*n),ainv(n*n),b(n),ipivot(n)
      isin=0
      call factor(a,ipivot,b,n,iflag)
      if(iflag.eq.2)then
         write(*,*)'singular matrix'
         isin=1
	 return
      endif
      do 1 i=1,n
         b(i)=0
1     continue
      ibeg=1
      do 2 j=1,n
         b(j)=1
         call subst(a,b,ainv(ibeg),ipivot,n)
         b(j)=0
         ibeg=ibeg+n
2     continue
      return
      end

      subroutine factor(a,ipivot,d,n,iflag)
      implicit real*8(a-h,o-z)
      dimension a(n,n),ipivot(n),d(n)
      iflag=1
      do 1 i=1,n
         ipivot(i)=i
         rowmax=0
         do 2 j=1,n
            rowmax=dmax1(rowmax,dabs(a(i,j)))
2        continue
         if(rowmax.le.0)then
            iflag=2
            return
         endif
         d(i)=rowmax
1     continue
      nm1=n-1
      if(nm1.eq.0)return
      do 3 k=1,nm1
         j=k
         kp1=k+1
         ip=ipivot(k)
         colmax=dabs(a(ip,k))/d(ip)
         do 4 i=kp1,n
            ip=ipivot(i)
            awikov=dabs(a(ip,k))/d(ip)
            if(awikov.le.colmax)go to 4
            colmax=awikov
            j=i
4        continue
         if(colmax.le.0)then
            iflag=2
            return
         endif
         ipk=ipivot(j)
         ipivot(j)=ipivot(k)
         ipivot(k)=ipk
         do 5 i=kp1,n
            ip=ipivot(i)
            a(ip,k)=a(ip,k)/a(ipk,k)
            ratio=-a(ip,k)
            do 6 j=kp1,n
               a(ip,j)=ratio*a(ipk,j)+a(ip,j)
6           continue
5        continue
3     continue
      return
      end

      subroutine subst(a,b,x,ipivot,n)
      implicit real*8(a-h,o-z)
      dimension a(n,n),b(n),x(n),ipivot(n)
      if(n.gt.1)go to 1
      x(1)=b(1)/a(1,1)
      return
1     ip=ipivot(1)
      x(1)=b(ip)
      do 2 k=2,n
         ip=ipivot(k)
         km1=k-1
         sum=0
         do 3 j=1,km1
            sum=a(ip,j)*x(j)+sum
3        continue
         x(k)=b(ip)-sum
2     continue
      x(n)=x(n)/a(ip,n)
      k=n
      do 4 np1mk=2,n
         kp1=k
         k=k-1
         ip=ipivot(k)
         sum=0
         do 5 j=kp1,n
            sum=a(ip,j)*x(j)+sum
5        continue
         x(k)=(x(k)-sum)/a(ip,k)
4     continue
      return
      end

      subroutine solin(a,b,x,d,ipivot,n,isin)
      implicit real*8(a-h,o-z)
      dimension a(n,n),b(n),x(n),d(n),ipivot(n)
      isin=0
      call factor(a,ipivot,d,n,iflag)
      if(iflag.eq.2)then
         write(*,*)' iflag=',iflag
         isin=1
	 return
      endif
      call subst(a,b,x,ipivot,n)
      return
      end
