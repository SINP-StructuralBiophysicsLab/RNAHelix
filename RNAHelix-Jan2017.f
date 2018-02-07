C         PROGRAM TO GENERATE A BENT OR STRAIGHT HELIX
C   RNAHelix version 1.0 -- D. Bhattacharyya, M. Bansal, S. Halder, D. Mukherje
C   Latest modification on Sept. 27, 2016
C

      COMMON /COORDS/XTM1(10000),YTM1(10000),ZTM1(10000),NRES(10000),
     1        ATOM(10000,2),ISTS(5000),IENS(5000)
      COMMON /GENE/ATMA(40,2),CA(3,40),ATMB(40,2),CB(3,40),ATMC(40,2),
     1        CC(3,40),ATMD(40,2),CD(3,40)
      DOUBLE PRECISION THXB,THXT,ROLLB,ROLLT,THYB,BUCK1,BUCK2,
     1                 TILTSAVE,ROLLSAVE,PROPSAVE,BUCKSAVE
      DIMENSION TILTSAVE(16),ROLLSAVE(16),TWISTSAV(16),DISPSAVE(16),
     1          SLIPSAVE(16),HSAVE(16),PROPSAVE(16),BUCKSAVE(16)
      DIMENSION XA(40),YA(40),ZA(40),XB(40),YB(40),ZB(40)
      DIMENSION XC(40),YC(40),ZC(40),XD(40),YD(40),ZD(40),BASTEP(2,16)
      DIMENSION X1(40),Y1(40),Z1(40),X2(40),Y2(40),Z2(40),X3(40),Y3(40),
     1          Z3(40),X4(40),Y4(40),Z4(40)
      DIMENSION BTEM(3),RMAT(3,3),BASTYP(4),YAXIS(3)
      CHARACTER*4 ATOMA(40,2),ATOMB(40,2),ATOMC(40,2),ATOMD(40,2),ATOM
      CHARACTER*4 ATMA,ATMB,ATMC,ATMD,dumm
      CHARACTER*1 STEP(2),BASTEP,BASTYP,B1,B2,BASESEQ(5000),SYMSEQ(5000)
     1     ,edge1,edge2,orient,BASEPERT(5000)
      CHARACTER*3 atmat(40),atmbt(40)
      CHARACTER*8 BASEPR,valp
      CHARACTER*80 LINE,outname,val1(10)
        external matmul

      DATA NAA,NAT,NAG,NAC/18,18,19,16/
      DATA YAXIS/0.000000,1.000000,0.0/
      DATA BASTYP/'A','T','G','C'/
      nout=0
!***********************************************************************
      write(*,*) '              ______________________________'
      write(*,*) '              |                            |'
      write(*,*) '              |          RNAHElix          |'
      write(*,*) '              |        Version -1.0        |'
      write(*,*) '              |____________________________|'
      write(*,*) ''
      write(*,*) '    Written by D. Bhattacharyya, M. Bansal, S. Haldar,
     1 D. Mukherjee'
      write(*,*) ''
      write(*,*) 'Computational Science Division, Saha Institute of',
     1' Nuclear Physics, INDIA. '
      write(*,*) ''
      write(*,*) '                 Last Update May 19, 2016'
      write(*,*) ''
      write(*,*)
     1'  This program is free software; you can redistribute it and/or'
      write(*,*)'  modify it under the terms of the  GNU General Public
     1License as ',
     2 '  published by the Free Software oundation; either version of '
     3 ,'  the License, or (at your option) any later version.'
C************************************************************************
C------------------------------------------------------------------------
      NYOPT=0
      incr=0
      outname='allatoms.pdb'
      CONV = 180.0/3.14159
      iagc= iargc()
c        write(*,*) 'No. of Arguments',iagc
      if(iagc.gt.1) then
        do i=1,iagc
          call getarg(i,val1(i))
        enddo
        if(iagc.eq.2) then
          if (val1(1)(1:1).eq.'-') then
            if(val1(2)(2:4).eq.'inp'.or.val1(i)(2:4).eq.'par'.or.
     1          val1(i)(2:4).eq.'loc') then
              line=val1(i+1)
c        write(6,*) 'Input filename',line
              open(unit=2,file=line,status='OLD',ERR=399)
            elseif(val1(2)(2:4).eq.'out'.or.val1(i)(2:4).eq.'pdb') then
              outname=val1(i+1)
              nout=1
c        write(6,*) 'Output filename',outname
            endif
          endif
        elseif(iagc.eq.4) then
          i=1
          do while(i.le.iagc)
            if (val1(i)(1:1).eq.'-') then
              if(val1(i)(2:4).eq.'inp'.or.val1(i)(2:4).eq.'par') then
                line=val1(i+1)
c        write(6,*) 'Input filename',line
                open(unit=2,file=line,status='OLD',ERR=399)
                i=i+2
              elseif(val1(i)(2:4).eq.'out'.or.val1(i)(2:4).eq.'pdb')then
                outname=val1(i+1)
                nout=1
c        write(6,*) 'Output filename',outname
                i=i+2
              endif
            endif
          enddo

        endif
c        write(6,*) 'File for Input parameter',line
c        write(6,*) 'File for Output coordinates',outname
      elseif(iagc.eq.1) then
c        write(*,*) 'No. of Arguments',iagc,line
        call getarg(1,line)
        open(unit=2,file=line,status='OLD',ERR=399)
      else
        OPEN(2,FILE='parameter.loc',status='OLD',ERR=399)
      endif
      nn=index(outname,' ')
      write(*,*) ''
      write(*,*) '        Coordinates of the generated model are stored
     1in : ',outname(1:nn)
      write(*,*) ''
      write(*,*) '          ==================OOO================='
      write(*,*) ''
c      OPEN(9,FILE='centers.bpc')
      open(11,file='double.hlx')
      KOUNT = 0
      FLENGTH = 0.0
      KO = 0
      READ(2,*) NBASE,NREPET,NYOPT

      IF(NREPET.NE.0) THEN
        OPEN(3,FILE='sequence.ply')
        DO 904 KK=1,16
          READ(2,712)BASEPR,WTILT,WROLL,TWI,SLX,SLY,SLZ
712       FORMAT(A8,2X,6F8.2)
          IF(TWI.EQ.0.0E0) TWI = 36.0
          if(slz.eq.0.0e0) slz = 3.4
          IF (WTILT.EQ.0.00) WTILT = 0.0001
          IF (WROLL.eq.0.00) WROLL = 0.0001

          ST = SIN(WTILT/(2.0 * CONV))
          SR = SIN(WROLL/(2.0 * CONV))
          SG = SQRT(ST * ST + SR * SR)
          GAMAB2 = ASIN(SG)
          CG = COS(GAMAB2)
          STW = SIN(TWI/(2.0 * CONV))
          STWIST = SQRT(STW * STW * CG * CG + SG * SG)
          HTWIST = ASIN(STWIST) * CONV * 2.0
          IF(TWI.LT.0.0) HTWIST = - HTWIST
          
          CALL FUNCT(WTILT,WROLL,HTWIST,THXB,ROLLB)
          CALL SLIDCA(THXB,ROLLB,HTWIST,SLX,SLY,SLZ,DISP,SLIP,HT)

          BASTEP(1,KK) = BASEPR(1:1)
          BASTEP(2,KK) = BASEPR(3:3)

          TILTSAVE(KK) = THXB
          ROLLSAVE(KK) = ROLLB
          TWISTSAV(KK) = HTWIST
          DISPSAVE(KK) = DISP
          SLIPSAVE(KK) = SLIP
          HSAVE(KK) = HT
      write(*,713)BASEPR(1:1),BASEPR(3:3),THXB,ROLLB,HTWIST,DISP,SLIP,HT
 904    CONTINUE
          
        IBS = 0
        DO 9050 KK=1,1000
          READ(3,9051,END=9052) LINE
          DO 9053 K=1,80
            IF(IBS.LE.(NBASE+2)) THEN
              IF(LINE(K:K).EQ.'A'.OR.LINE(K:K).EQ.'a') THEN
                IBS = IBS + 1
                BASESEQ(IBS)='A'
              END IF
              IF(LINE(K:K).EQ.'G'.OR.LINE(K:K).EQ.'g') THEN
                IBS = IBS + 1
                BASESEQ(IBS)='G'
              END IF
              IF(LINE(K:K).EQ.'C'.OR.LINE(K:K).EQ.'c') THEN 
                IBS = IBS + 1
                BASESEQ(IBS)='C'
              END IF
              IF(LINE(K:K).EQ.'T'.OR.LINE(K:K).EQ.'t') THEN
                IBS = IBS + 1
                BASESEQ(IBS)='T'
              END IF
            END IF
 9053     CONTINUE
 9050   CONTINUE
 9051   FORMAT(A80)
 9052   CONTINUE
 9054   FORMAT(1X,60A1)
      END IF       

      DO 10 KOUNT=1,NBASE
        IF(NREPET.EQ.0) THEN

          THXB = 0.0D0
          ROLLB = 0.0D0
          READ(2,711)BASEPR,WTILT,WROLL,TWI,SLX,SLY,SLZ,
     1         BUCK1s,openan,thybs,stag1,shear1,strc1
711       FORMAT(A8,2X,12F8.2)

!===========Added to fetch basepairing type====================DM

          BASESEQ(KOUNT)=BASEPR(1:1)
          IF(BASEPR(5:5) .EQ. 'W' .OR. BASEPR(5:5) .EQ. 'w') THEN
            BASEPERT(KOUNT)= '|'
          END IF
          IF(BASEPR(5:5) .EQ. 'H' .OR. BASEPR(5:5) .EQ. 'h') THEN
            BASEPERT(KOUNT)= '*'
          END IF
          IF(BASEPR(5:5) .EQ. 'S' .OR. BASEPR(5:5) .EQ. 's') THEN
            BASEPERT(KOUNT)= '-'
          END IF
          BASESEQ(KOUNT+NBASE)=BASEPR(3:3)
          IF(BASEPR(7:7) .EQ. 'W' .OR. BASEPR(7:7) .EQ. 'w') THEN
            BASEPERT(KOUNT+NBASE)= '|'
          END IF
          IF(BASEPR(7:7) .EQ. 'H' .OR. BASEPR(7:7) .EQ. 'h') THEN
            BASEPERT(KOUNT+NBASE)= '*'
          END IF
          IF(BASEPR(7:7) .EQ. 'S' .OR. BASEPR(7:7) .EQ. 's') THEN
            BASEPERT(KOUNT+NBASE)= '-'
          END IF
!==============================================================DM
          IF(TWI.EQ.0.0E0) TWI = 36.0
          if(slz.eq.0.0e0) slz = 3.4
          IF (WTILT.EQ.0.00) WTILT = 0.0001
          IF (WROLL.EQ.0.00) WROLL = 0.0001

          ST = SIN(WTILT/(2.0 * CONV))
          SR = SIN(WROLL/(2.0 * CONV))
          SG = SQRT(ST * ST + SR * SR)
          GAMAB2 = ASIN(SG)
          CG = COS(GAMAB2)
          STW = SIN(TWI/(2.0 * CONV))
          STWIST = SQRT(STW * STW * CG * CG + SG * SG)
          HTWIST = ASIN(STWIST) * CONV * 2.0
          IF(TWI.LT.0.0) HTWIST = - HTWIST

          CALL FUNCT(WTILT,WROLL,HTWIST,THXB,ROLLB)
          CALL SLIDCA(THXB,ROLLB,HTWIST,SLX,SLY,SLZ,DISP,SLIP,HT)
c          write(*,*)WTILT,WROLL,TWI,SLX,SLY,SLZ
c          write(*,*)THXB,ROLLB,HTWIST,DISP,SLIP,HT
 
          B1 = BASEPR(1:1)
          B2 = BASEPR(3:3)
          edge1=BASEPR(5:5)
          edge2=BASEPR(7:7)
          orient=BASEPR(8:8)
          DO 906 KB=1,4
            IF(B1.EQ.BASTYP(KB)) NTY1 = KB
            IF(B2.EQ.BASTYP(KB)) NTY2 = KB
 906      CONTINUE
          NTY3 = NTY1
          NTY4 = NTY2
        ELSE
          STEP(1) = BASESEQ(KOUNT)
          STEP(2) = BASESEQ(KOUNT+1)
          DO 907 KK=1,16
      IF(STEP(1).EQ.BASTEP(1,KK).AND.STEP(2).EQ.BASTEP(2,KK)) GO TO 200
 907      CONTINUE
 200      CONTINUE

          THXB = TILTSAVE(KK)
          ROLLB = ROLLSAVE(KK)
          HTWIST = TWISTSAV(KK)
          DISP = DISPSAVE(KK)
          SLIP = SLIPSAVE(KK)
          HT = HSAVE(KK)
          write(*,713)STEP(1),STEP(2),THXB,ROLLB,HTWIST,DISP,SLIP,HT 
713       FORMAT(A1,1X,A1,2X,6F8.2)

          THYBS = -10.0
          BUCK1S = 0.001
          OPENAN = 0.001
          STAG1 = 0.001
          SHEAR1 = 0.001
          STRC1 = 2.80
           
          DO 908 KB=1,4
            IF(STEP(1).EQ.BASTYP(KB)) NTY1 = KB
 908      CONTINUE
          NT1 = NTY1 + 1
          NTY2 = NT1 - 2 * MOD(NT1,2)
          NTY3 = NTY1
          NTY4 = NTY2
          b1 = BASTYP(NTY1)
          b2 = BASTYP(NTY2)
          edge1 = 'W'
          edge2 = 'W'
          orient = 'C'
        END IF

        THXT = THXB
        ROLLT = ROLLB
  5     FORMAT(I3)
190     FORMAT(1X,A4,3F10.4,I3,F8.3)
        call genbp(b1,b2,edge1,edge2,orient,buck1s,openan,thybs,
     1  stag1,shear1,strc1,x1,y1,z1,atmat,x2,y2,z2,atmbt,na,nb,nyopt)
        NATOM = NA + NB
        KO = KO + NATOM

C  All the previously generated atoms are given the same amount of tip
C  rotation now.

        FLENGTH = FLENGTH + HT
        NT = (KO - NATOM)
        CALL ROTMAD(0.0,1.0,0.0,ROLLB,RMAT)
        DO 909 K=1,NT
          CALL MATMUL(RMAT,XTM1(K),YTM1(K),ZTM1(K),BTEM)
          XTM1(K) = BTEM(1)
          YTM1(K) = BTEM(2)
          ZTM1(K) = BTEM(3)
 909    CONTINUE
        do k=1,na
           call matmul(rmat,x1(k),y1(k),z1(k),btem)
           x1(k)=btem(1)
           y1(k)=btem(2)
           z1(k)=btem(3)
        enddo
        do k=1,nb
          call matmul(rmat,x2(k),y2(k),z2(k),btem)
          x2(k)=btem(1)
          y2(k)=btem(2)
          z2(k)=btem(3)
        enddo

C  All the previously generated atoms are given the same amount of
C  inclination rotation now.  The displacement parameters are also
C  applied here.

      CALL ROTMAD(1.0,0.0,0.0,THXB,RMAT)
      DO 910 K=1,NT
         CALL MATMUL(RMAT,XTM1(K),YTM1(K),ZTM1(K),BTEM)
          IF(HTWIST.GT.0.0) THEN
           XTM1(K) = BTEM(1) + DISP
           YTM1(K) = BTEM(2) + SLIP
          ELSE
            XTM1(K) = BTEM(1) - DISP
            YTM1(K) = BTEM(2) - SLIP
          END IF
         ZTM1(K) = BTEM(3)
 910    CONTINUE
        do k=1,na
          call matmul(rmat,x1(k),y1(k),z1(k),btem)
          if(htwist.gt.0.0) then
            x1(k)=btem(1) + disp
            y1(k)=btem(2) + slip
          else
            x1(k)=btem(1) - disp
            y1(k)=btem(2) - slip
          endif
          z1(k)=btem(3)
        enddo
        do k=1,nb
          call matmul(rmat,x2(k),y2(k),z2(k),btem)
          if(htwist.gt.0.0) then
            x2(k)=btem(1) + disp
            y2(k)=btem(2) + slip
          else
            x2(k)=btem(1) - disp
            y2(k)=btem(2) - slip
          endif
          z2(k)=btem(3)
        enddo
        ISTS(2*KOUNT-1) = NT + 1

C  The newly generated atoms are added to the previously generated set.

      DO 911 K=1,NA
        NAK = NT + K
        XTM1(NAK) = X1(K)
        YTM1(NAK) = Y1(K)
        ZTM1(NAK) = Z1(K)
        atom(nak,1) = atmat(k)//' '
        atom(nak,2) = '   '//B1
        NRES(NAK) = KOUNT * 2 -1
 911  CONTINUE
      IENS(2*KOUNT-1) = NT + NA
      ISTS(2*KOUNT) = NT + NA + 1
      NT = NT + NA
      DO 912 K =1,NB
        NAK = NT + K
        XTM1(NAK) = X2(K)
        YTM1(NAK) = Y2(K)
        ZTM1(NAK) = Z2(K)
        atom(nak,1) = atmbt(k)//' '
        atom(nak,2) = '   '//b2
        NRES(NAK) = KOUNT * 2 
 912  CONTINUE
      IENS(2*KOUNT) = NT + NB

C  Helical twist in the form of (-H,-T) has been applied to keep the
C  top base-pair neearly oriented along the global coordinate system.  All
C  the atoms are translated to the  origin at the same time.

      CALL ROTMAT(0.,0.,1.,-HTWIST,RMAT)
      DO 913 K=1,KO
        CALL MATMUL(RMAT,XTM1(K),YTM1(K),ZTM1(K),BTEM)
         IF(HTWIST.GT.0.0) THEN
          XTM1(K) = BTEM(1) - DISP
          YTM1(K) = BTEM(2) - SLIP
         ELSE
           XTM1(K) = BTEM(1) + DISP
           YTM1(K) = BTEM(2) + SLIP
         END IF
        ZTM1(K) = BTEM(3) - HT
 913  CONTINUE
193   FORMAT('ATOM',9X,A4,A3,2X,I4,4X,3F8.3)

C  Inclination of the last basepair is subtracted now from all the
C  atoms, including those generated now, to bring the long  axis of the
C  last basepair in the X-Y plane. 

        CALL ROTMAD(1.0,0.0,0.0,-THXT,RMAT)
        DO 914 I=1,KO
          CALL MATMUL(RMAT,XTM1(I),YTM1(I),ZTM1(I),BTEM)
          XTM1(I) = BTEM(1)
          YTM1(I) = BTEM(2)
          ZTM1(I) = BTEM(3)
 914    CONTINUE

C  Tip of the last basepar is subtracted now from all the atoms to
C  bring the last basepair in the X-Y plane.

        CALL ROTMAD(0.0,1.0,0.0,-ROLLT,RMAT)
        DO 915 KJ=1,KO
          CALL MATMUL(RMAT,XTM1(KJ),YTM1(KJ),ZTM1(KJ),BTEM)
          XTM1(KJ) = BTEM(1)
          YTM1(KJ) = BTEM(2)
          ZTM1(KJ) = BTEM(3)
 915    CONTINUE

C ---------------------------------------------------------------------
C            THE ABOVE PART OF THE PROGRAM ROTATES THE SET BY EXACTLY
C       THE OPPOSITE ROTATIONS APPLIED TO THE TWO BASE-PAIRS IN "DNPAIR"
C ----------------------------------------------------------------------

 10   CONTINUE


C      CALL FINDVEC(NBASE,EEDIST)
!      WRITE(6,*) 'Path length =',FLENGTH
C      RATIO = EEDIST / FLENGTH
c      WRITE(9,'(/'' Bending ratio "d/l" is ='',f10.5/)') RATIO
c      WRITE(6,'(/'' Bending ratio "d/l" is ='',f10.5/)') RATIO
c      WRITE(11,'(''DNA of sequence:'',20A1,''...'')')(BASESEQ(M),M=1,20)
      WRITE(11,'(''Sequence of Generated Nucleotides:'')')
      WRITE(11,*) '5',"'",'-',(BASESEQ(M),M=1,NBASE),'-3',"'"
      WRITE(11,*) '   ',(BASEPERT(i),i=1,NBASE)
      imb=NBASE+NBASE
      WRITE(11,*) '   ',(BASEPERT(j),j=NBASE+1,imb)
      WRITE(11,*) '3',"'",'-',(BASESEQ(M),M=NBASE+1,imb),'-5',"'",'\n'

      NUM = 0
      DO 918 I=1,10000
        IF(ATOM(I,1).EQ.'C1'' ') THEN
          NUM = NUM + 1
          IF(ATOM(I,2).EQ.'   A') SYMSEQ(NUM) = 'A'
          IF(ATOM(I,2).EQ.'   T') SYMSEQ(NUM) = 'T'
          IF(ATOM(I,2).EQ.'   G') SYMSEQ(NUM) = 'G'
          IF(ATOM(I,2).EQ.'   U') SYMSEQ(NUM) = 'U'
          IF(ATOM(I,2).EQ.'   C') SYMSEQ(NUM) = 'C'
          WRITE(11,'(3F10.3,3X,A1)') XTM1(I),YTM1(I),ZTM1(I),SYMSEQ(NUM)
        END IF
 918    CONTINUE
        DO 20 I=1,NUM,2
 20     CONTINUE
        numold=0
        if(nout.eq.0) then
          OPEN(7,FILE='allatoms.pdb')
        else
          open(unit=7,file=outname)
        endif
        DO 916 K=1,KO
           NRESODD=(NRES(K)/2) * 2
        IF(NRES(K).NE.NRESODD) THEN
          if(nres(k).ne.numold) then
            numold=nres(k)
            incr=incr+1
          endif
        dumm=atom(k,2)
        WRITE(7,193) ATOM(K,1),dumm(2:4),incr,XTM1(K),YTM1(K),ZTM1(K)
        END IF
 916    CONTINUE
        DO 917 K=KO,1,-1
         NRESEVN = (NRES(K)/2) * 2
         IF(NRES(K).EQ.NRESEVN) THEN
           if(nres(k).ne.numold) then
             numold=nres(k)
             incr=incr+1
           endif
        dumm=atom(k,2)
        WRITE(7,193) ATOM(K,1),dumm(2:4),incr,XTM1(K),YTM1(K),ZTM1(K)
        END IF
 917    CONTINUE         
        write(7,73)'END'
73      format(a3)
        CLOSE (UNIT=7)
	stop
399     write(*,398) 
398    format('Please create parameter.loc file and then run RNAHelix')
      STOP
      END

C***********************************************************************
      SUBROUTINE ROTMAD(EL,EM,EN,THETA,R)
  
C  This subroutine calculates the matrix elements corresponding to
C  rotation about the vector (EL,EM,EN) through an angle THETA.

      DOUBLE PRECISION THETA,CST,SNT,CSTF,THETR
      DIMENSION R(3,3)
      CON=3.14195/180.0
      THETR=THETA*CON
      CST=DCOS(THETR)
      SNT=DSIN(THETR)
      CSTF=1.0D0-CST
      R(1,1)=EL*EL*CSTF+CST
      R(2,2)=EM*EM*CSTF+CST
      R(3,3)=EN*EN*CSTF+CST
      AA=EL*EM*CSTF
      BB=EN*SNT
      R(1,2)=AA-BB
      R(2,1)=AA+BB
      AA=EL*EN*CSTF
      BB=EM*SNT
      R(1,3)=AA+BB
      R(3,1)=AA-BB
      AA=EM*EN*CSTF
      BB=EL*SNT
      R(2,3)=AA-BB
      R(3,2)=AA+BB
      RETURN
      END
C***********************************************************************

      SUBROUTINE FUNCT(WTLT,WRLL,TWIST,TILT,ROLL)

C  Calulation of basepair inclination and basepair tip angles from the
C  local doublet parameters, tilt, roll and twist using the  the
C  expression derived. 

      DOUBLE PRECISION SR,ST,RT1,RT4,BRACET,TANT2,TL,RL,X,Y,X1,Y1
      DOUBLE PRECISION WTILT,WROLL,TWI,TILT,ROLL,TANT22
      CONV=180.0/3.1415926
      IF(WTLT.EQ.0.00) WTLT = 0.0001
      IF(WRLL.EQ.0.00) WRLL = 0.0001
      WTILT = WTLT/CONV
      WROLL = WRLL/CONV
      TWI   = TWIST/CONV
      SR =  DSIN(WROLL/2.0D0)
      ST =  DSIN(WTILT/2.0D0)
      SR = SR * SR
      ST = ST * ST
      RT1 = 1.0 - SR - ST
      RT4 = SR * ST * 4.0
      BRACET = RT1 -  DSQRT(RT1 * RT1 - RT4)
      TANT2 = DTAN(TWI/2.0D0)
      TANT22 = TANT2 * TANT2
      Y1 = BRACET/(2.0 * ST * TANT22)
      IF(Y1.LT.0.0D0) THEN
        Y = 0.0D0
      ELSE
          Y = DSQRT(Y1)
      END IF
      IF(DABS(Y).GT.1.0D0) THEN
          TL = 0.0D0
      ELSE
           TL = DASIN(Y) * CONV
      END IF
      TILT =  TL * 1.0D0
      X1 = BRACET/(2.0 * SR * TANT22)
      IF(X1.LT.0.0D0) THEN
         X = 0.0D0
      ELSE
         X = DSQRT(X1)
      END IF
      IF(DABS(X).GT.1.0D0) THEN
         RL = 0.0D0
      ELSE
           RL = DASIN(X) * CONV
      END IF
      ROLL =  RL * (-1.0D0)

      IF(WTILT.LT.0.0D0) ROLL = -ROLL
      IF(WROLL.LT.0.0D0) TILT = -TILT

      RETURN
      END
C***********************************************************************

      SUBROUTINE ROTMAT(EL,EM,EN,THETA,R)
      DIMENSION R(3,3)
      CON=3.14195/180.0
      THETR=THETA*CON
      CST= COS(THETR)
      SNT= SIN(THETR)
      CSTF=1.0-CST
      R(1,1)=EL*EL*CSTF+CST
      R(2,2)=EM*EM*CSTF+CST
      R(3,3)=EN*EN*CSTF+CST
      AA=EL*EM*CSTF
      BB=EN*SNT
      R(1,2)=AA-BB
      R(2,1)=AA+BB
      AA=EL*EN*CSTF
      BB=EM*SNT
      R(1,3)=AA+BB
      R(3,1)=AA-BB
      AA=EM*EN*CSTF
      BB=EL*SNT
      R(2,3)=AA-BB
      R(3,2)=AA+BB
      RETURN
      END
C***********************************************************************
      SUBROUTINE SLIDCA(TIL,ROL,TWIS,SLIDEX,SLIDEY,SLIDEZ,D1,D2,D3)

C  This subroutine calculates the local helical displacement dx,dy and
C  dz from the local displacements, shift (Dx), slide (Dy) and rise (Dz).

      DIMENSION N(3),M(3),A(3,3)
      DOUBLE PRECISION  TIL,ROL
      CONV = 3.1415926/180.000
      ROLL = ROL * CONV
      IF(TWIS.LT.0.0) THEN
        TILT = -1.0D0 * TIL * CONV
      ELSE
        TILT = TIL * CONV
      END IF
      TWIST = TWIS * CONV
      CY = COS(-ROLL)
      CX = COS(TILT)
      SY = SIN(-ROLL)
      SX = SIN(TILT)
      CT = COS(TWIST)
      ST = SIN(TWIST)
      A(1,1) = 2.0 * CX * ST
      A(1,2) = 0.0
      A(1,3) = 2.0 * SX
      A(2,1) = 0.0
      A(2,2) = -2.0 * CY * ST
      A(2,3) = 2.0 * SY
      A(3,1) = -2.0 * CY * SX * ST
      A(3,2) = 2.0 * CX * SY * ST
      A(3,3) = CX * CY * (1.0 + CT)
      CALL MINV(A,3,DXX,N,M)
      DX = 1.0 + CT + SX * SX * ( 1.0 - CT)
      DX = SQRT(2.0 * DX)
      DY = 1.0 + CT + SY * SY * (1.0 - CT)
      DY = SQRT(2.0 * DY)
      B1 = SLIDEY * DY
      B2 = SLIDEX * DX
      B3 = (SLIDEZ * DX * DY)/2
      D1 = A(1,1) * B1 + A(1,2) * B2 + A(1,3) * B3
      D2 = A(2,1) * B1 + A(2,2) * B2 + A(2,3) * B3
      D3 = A(3,1) * B1 + A(3,2) * B2 + A(3,3) * B3
      RETURN
      END

C***********************************************************************

      SUBROUTINE FINDVEC(NTIMES,DIST)
  
C  This subroutine calculates the basepair centers of all the
C  generated basepairs and the end-to-end distance, i.e. distance
C  between the first and last basepair centers.

        COMMON /COORDS/XTM1(10000),YTM1(10000),ZTM1(10000),NRES(10000),
     1        ATOM(10000,2),ISTS(5000),IENS(5000)
      DIMENSION ORIG(3), FIRST(3),LAST(3)
      CHARACTER*4 ATOM
      REAL LAST

        DO  201 KKTIME =1,NTIMES
         NF = KKTIME * 2 - 1
         NS = KKTIME * 2

      IF(ATOM(ISTS(NF),2).EQ.'   G'.OR.ATOM(ISTS(NF),2).EQ.'   A') THEN
             DO 52 J=ISTS(NF),IENS(NF)
                IF(ATOM(J,1).EQ.'C8  ') NC8 = J
 52           CONTINUE
           ELSE
              DO 53 J=ISTS(NF),IENS(NF)
                IF(ATOM(J,1).EQ.'C6  ') NC6 = J
  53          CONTINUE
          END IF
      IF(ATOM(ISTS(NS),2).EQ.'   G'.OR.ATOM(ISTS(NS),2).EQ.'   A') THEN
             DO 50 I=ISTS(NS),IENS(NS)
                IF(ATOM(I,1).EQ.'C8  ') NC8 = I
 50          CONTINUE
      ELSE
            DO 51 I=ISTS(NS),IENS(NS)
              IF(ATOM(I,1).EQ.'C6  ') NC6 = I
  51          CONTINUE
      END IF
C
  1     FORMAT(1X,A4,3F10.4,I4)
           ORIG(1) = (XTM1(NC8) + XTM1(NC6)) * 0.5
           ORIG(2) = (YTM1(NC8) + YTM1(NC6)) * 0.5
           ORIG(3) = (ZTM1(NC8) + ZTM1(NC6)) * 0.5
            IF(KKTIME.EQ.1) THEN
               FIRST(1) = ORIG(1)
        	FIRST(2) = ORIG(2)
        	FIRST(3) = ORIG(3)
            END IF
            IF(KKTIME.EQ.NTIMES) THEN
              LAST(1) = ORIG(1)
              LAST(2) = ORIG(2)
              LAST(3) = ORIG(3)
            END IF
C
c            WRITE(9,6) (ORIG(KJ),KJ=1,3)
 201    CONTINUE
 	DIST = SQRT( (FIRST(1)-LAST(1))**2 + (FIRST(2)-LAST(2))**2 +
     1              (FIRST(3)-LAST(3))**2)
c       WRITE(9,'('' End-to-end Distance is'',f10.3)') DIST
  6     FORMAT(3F12.4)
      RETURN
      END

C***********************************************************************
 
C***********************************************************************
      subroutine genbp(base1,base2,edge1,edge2,orient,buck,openan,
     1  prop,stag,shear,strch,xt1,yt1,zt1,atpt1,xt2,yt2,zt2,atpt2,
     2  nat1,nat2,nyopt)
      
      external MATMUL
      character*3 basepair,type,atmpr(80),residue(15,40),
     1  respr(40),atpt1(40),atpt2(40)
      character*1 base1,base2,edge1,edge2,orient
      character*4 atomnm(15,40)
      character*80 line
      real mid_c1c1(3)
      double precision thxb, rollb,exl,exm,exn,eyl,eym,eyn,ezl,ezm,ezn,
     1  el1,em1,en1,el2,em2,en2,yaxis(3),ezl12,ezm12,ezn12,y_c1c1(3),
     2  el11,em11,en11,el22,em22,en22,exl1,exm1,exn1,eyl1,eym1,eyn1,
     3  ezl1,ezm1,ezn1,ydir(3),ydir1(3)
      dimension xyzbase(15,40,3), rotmats(15,3,3) ,transl(15,3),
     1  natmbase(15),xyzpr(80,3),b(3),rmat(3,3),xt1(40),xt2(40),yt1(40),
     2  yt2(40),zt1(40),zt2(40),xc11(3),xc12(3),np1(9999),np2(9999),
     3  xcor1(3,9999),xcor2(3,9999),xcor11(3,9999),xcor22(3,9999)
      dimension yo(3),yp(3)

      conv = 180.0/3.14159
      open(unit=1,file='/usr/local/bin/DataSet.dat',STATUS='OLD',ERR=98)

      jcontr=0
      do while(jcontr.eq.0)
        read(1,1,end=99) line
        if(line(1:6).eq.'COMPND') then
          read(line,3) natombs
          if(line(9:11).eq.'ADE') then
            if(line(13:13).eq.'W'.or.line(13:13).eq.'+'.or.line(13:13)
     1   .eq.'w') then
              k=1
              natmbase(k)=natombs
            elseif(line(13:13).eq.'H'.or.line(13:13).eq.'h'.or.
     1    line(13:13).eq.'g') then
              k=2
              natmbase(k)=natombs
            elseif(line(13:13).eq.'S'.or.line(13:13).eq.'s'.or.
     1    line(13:13).eq.'z') then
              k=3
              natmbase(k)=natombs
            endif
          elseif(line(9:11).eq.'GUA') then
            if(line(13:13).eq.'W'.or.line(13:13).eq.'w') then
              k=4
              natmbase(k)=natombs
            elseif(line(13:13).eq.'H'.or.line(13:13).eq.'h'.or.
     1    line(13:13).eq.'g') then
              k=5
              natmbase(k)=natombs
            elseif(line(13:13).eq.'S'.or.line(13:13).eq.'s'.or.
     1    line(13:13).eq.'z') then
              k=6
              natmbase(k)=natombs
            endif
          elseif(line(9:11).eq.'CYT') then
            if(line(13:13).eq.'W'.or.line(13:13).eq.'+'.or.line(13:13)
     1   .eq.'w') then
              k=7
              natmbase(k)=natombs
            elseif(line(13:13).eq.'H'.or.line(13:13).eq.'h'.or.
     1   line(13:13).eq.'g') then
              k=8
              natmbase(k)=natombs
            elseif(line(13:13).eq.'S'.or.line(13:13).eq.'s'.or.
     1   line(13:13).eq.'z') then
              k=9
              natmbase(k)=natombs
            endif
          elseif(line(9:11).eq.'URA') then
            if(line(13:13).eq.'W'.or.line(13:13).eq.'w') then
              k=10
              natmbase(k)=natombs
            elseif(line(13:13).eq.'H'.or.line(13:13).eq.'h'.or.
     1   line(13:13).eq.'g') then
              k=11
              natmbase(k)=natombs
            elseif(line(13:13).eq.'S'.or.line(13:13).eq.'s'.or.
     1   line(13:13).eq.'z') then
              k=12
              natmbase(k)=natombs
            endif
          elseif(line(9:11).eq.'THY') then
            if(line(13:13).eq.'W'.or.line(13:13).eq.'w') then
              k=13
              natmbase(k)=natombs
            elseif(line(13:13).eq.'H') then
              k=14
              natmbase(k)=natombs
            elseif(line(13:13).eq.'S'.or.line(13:13).eq.'s') then
              k=15
              natmbase(k)=natombs
            endif
          endif 
          read(1,2) (transl(k,l),l=1,3)
          read(1,2) ((rotmats(k,l,m),m=1,3),l=1,3)
          do i=1,natombs
            read(1,4) atomnm(k,i),residue(k,i),(xyzbase(k,i,l),l=1,3)
          enddo
        endif
      enddo
99    continue
27    format(a1,1x,a1,1x,a1,1x,a1,1x,a1,6f8.2)

      if(base1.eq.'A') ib1 = 1
      if(base1.eq.'G') ib1 = 2
      if(base1.eq.'C') ib1 = 3
      if(base1.eq.'U') ib1 = 4
      if(base1.eq.'T') ib1 = 5
      if((edge1.eq.'W').or.(edge1.eq.'w').or.edge1.eq.'+') ie1 = 1
      if((edge1.eq.'H').or.(edge1.eq.'h')) ie1 = 2
      if((edge1.eq.'S').or.(edge1.eq.'s')) ie1 = 3
      if(base2.eq.'A') ib2 = 1
      if(base2.eq.'G') ib2 = 2
      if(base2.eq.'C') ib2 = 3
      if(base2.eq.'U') ib2 = 4
      if(base2.eq.'T') ib2 = 5
      if((edge2.eq.'W').or.(edge2.eq.'w').or.edge2.eq.'+') ie2 = 1
      if((edge2.eq.'H').or.(edge2.eq.'h')) ie2 = 2
      if((edge2.eq.'S').or.(edge2.eq.'s')) ie2 = 3
      nbseg1 = (ib1-1)*3+ie1
      nbseg2 = (ib2-1)*3+ie2

      do i=1,natmbase(nbseg1)
        xyzpr(i,1) = xyzbase(nbseg1,i,1)
        xyzpr(i,2) = xyzbase(nbseg1,i,2)
        xyzpr(i,3) = xyzbase(nbseg1,i,3)
        atmpr(i) = atomnm(nbseg1,i)
        respr(i) = residue(nbseg1,i)
      enddo
      num=natmbase(nbseg1)
      do i=1,natmbase(nbseg2)
        xyzpr(i+num,1)=xyzbase(nbseg2,i,1)
        xyzpr(i+num,2)=-xyzbase(nbseg2,i,2)
        xyzpr(i+num,3)=-xyzbase(nbseg2,i,3)
        atmpr(i+num) = atomnm(nbseg2,i)
        respr(i+num) = residue(nbseg2,i)
      enddo
      if(orient.eq.'T') then
        do i=1,natmbase(nbseg2)
          xyzpr(i+num,1)=-xyzpr(i+num,1)
          xyzpr(i+num,3)=-xyzpr(i+num,3)
        enddo
      endif
      num=num+natmbase(nbseg2)

      k=1
!      write(*,*) 'No. of atoms in base pair',num

C Assignment of Idealized Coordinates of both the bases COMPLETE
c=======================================================================
c Conversion of BL parameters in LOCAL frame to LOCAL HELICAL frame

      if(buck.eq.0.0) buck=0.001
      if(openan.eq.0.0) openan=0.001
      if(prop.eq.0.0) prop=0.001
      if(stag.eq.0.0) stag=0.001
      if(shear.eq.0.0) shear=0.001
      if(prop.gt.0.0)buck = -buck   ! probably relevent for Cis bps only
      if(orient.eq.'T') then
        buck=-buck
        openan=-openan
        stag=-stag
        shear=-shear
      endif

      ST = SIN(buck/(2.0 * CONV))
      SR = SIN(openan/(2.0 * CONV))
      SG = SQRT(ST * ST + SR * SR)
      GAMAB2 = ASIN(SG)
      CG = COS(GAMAB2)
      STW = SIN(prop/(2.0 * CONV))
      STWIST = SQRT(STW * STW * CG * CG + SG * SG)
      HTWIST = ASIN(STWIST) * CONV * 2.0
      CALL FUNCT(buck,openan,HTWIST,THXB,ROLLB)
      if(prop.lt.0.0) then
        htwist=-htwist
        thxb = -thxb
        CALL SLIDCA(-ROLLB,-THXB,HTWIST,stag,shear,strch,DISP,SLIP,HT)
      else
        CALL SLIDCA(ROLLB,-THXB,HTWIST,stag,shear,strch,DISP,SLIP,HT)
      endif

c   Rotation is applied corresponding to OPEN angle

      angle = rollb *(1.0)
      call rotmat(0.0,0.0,1.0,angle,rmat)
      do i=1,natmbase(nbseg1)
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
      enddo   
      call rotmat(0.0,0.0,1.0,angle,rmat)
      do i=natmbase(nbseg1)+1,num
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
      enddo   

c   Rotation is applied corresponding to BUCKLE

      angle = thxb*(-1.0) 
      call rotmat(1.0,0.0,0.0,angle,rmat)
      do i=1,natmbase(nbseg1)
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
      enddo
      call rotmat(1.0,0.0,0.0,angle,rmat)
      do i=natmbase(nbseg1)+1,num
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
      enddo

c   Translational parameters are applied here 

      do i=1,natmbase(nbseg1)
        xyzpr(i,3) = xyzpr(i,3) + disp
        xyzpr(i,1) = xyzpr(i,1) + slip
      enddo
      do i=natmbase(nbseg1)+1,num
        xyzpr(i,3) = xyzpr(i,3) + disp
        xyzpr(i,1) = xyzpr(i,1) + slip
      enddo

c Rotation & Translation applied corresponding to PROPELLER & STRETCH

      angle = htwist/2.0
      call rotmat(0.0,1.0,0.0,angle,rmat)
      do i=1,natmbase(nbseg1)
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
        xyzpr(i,2)=xyzpr(i,2)+ht/2.0
      enddo
      angle=-angle
      call rotmat(0.0,1.0,0.0,angle,rmat)
      do i=natmbase(nbseg1)+1,num
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
        xyzpr(i,2)=xyzpr(i,2)-ht/2.0
      enddo

C     Reverse Translation

      do i=1,num
        xyzpr(i,3)=xyzpr(i,3)-disp
        xyzpr(i,1)=xyzpr(i,1)-slip
      enddo

C     Reverse Rotation

      angle=thxb*1.0
      call rotmat(1.0,0.0,0.0,angle,rmat)
      do i=1,num
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
      enddo

      angle=-rollb*1.0
      call rotmat(0.0,0.0,1.0,angle,rmat)
      nsecond=0
      do i=1,num
        call matmul(rmat,xyzpr(i,1),xyzpr(i,2),xyzpr(i,3),b)
        do k=1,3
          xyzpr(i,k)=b(k)
        enddo
        if(i.le.natmbase(nbseg1)) then
          np1(i)=0
          xt1(i)=xyzpr(i,1)
          yt1(i)=xyzpr(i,2)
          zt1(i)=xyzpr(i,3)
          atpt1(i)=atmpr(i)
          do nz=1,3
            xcor1(nz,i)=xyzpr(i,nz)
          enddo
          if(nyopt.gt.0)then
            if(atpt1(i).eq.'C1''') nc11=i 
          else
            if(((respr(i).eq.'ADE').or.(respr(i).eq.'GUA'))
     1      .and.(atpt1(i).eq.'C8'))nc11=i
            if(((respr(i).eq.'CYT').or.(respr(i).eq.'URA')
     1      .or.(respr(i).eq.'THY')).and.(atpt1(i).eq.'C6'))nc11=i
          endif
          if((atpt1(i).ne.'C1''').and.(atpt1(i)(1:1).ne.'H'))np1(i)=i 
        else
          nsecond=nsecond+1
          np2(nsecond)=0
          xt2(nsecond)=xyzpr(i,1)
          yt2(nsecond)=xyzpr(i,2)
          zt2(nsecond)=xyzpr(i,3)
          atpt2(nsecond)=atmpr(i)
          do nz=1,3
            xcor2(nz,nsecond)=xyzpr(i,nz)
          enddo
          if(nyopt.gt.0)then
            if(atpt2(nsecond).eq.'C1''') nc12=i
          else
            if(((respr(i).eq.'ADE').or.(respr(i).eq.'GUA'))
     1      .and.(atpt2(nsecond).eq.'C8')) nc12=i
            if(((respr(i).eq.'CYT').or.(respr(i).eq.'URA')
     1      .or.(respr(i).eq.'THY')).and.(atpt2(nsecond).eq.'C6'))nc12=i
          endif
          if((atpt2(nsecond).ne.'C1''').and.(atpt2(nsecond)(1:1).ne.
     1    'H')) np2(nsecond)=nsecond
        endif
      enddo
      nat1=natmbase(nbseg1)
      nat2=natmbase(nbseg2)

C Determining Direction of Major Groove

      do i=1,nat1
        if((base1.eq.'A').or.(base1.eq.'G'))then
          if(edge1.eq.'W')then
            if(atpt1(i).eq.'N1')then
              yo(1)=xt1(i)
              yo(2)=yt1(i)
              yo(3)=zt1(i)
            elseif((atpt1(i).eq.'N6').or.(atpt1(i).eq.'O6'))then
              yp(1)=xt1(i)
              yp(2)=yt1(i)
              yp(3)=zt1(i)
            endif
          elseif(edge1.eq.'H'.or.edge1.eq.'h')then
            if(atpt1(i).eq.'N7')then
              yo(1)=xt1(i)
              yo(2)=yt1(i)
              yo(3)=zt1(i)
            elseif((atpt1(i).eq.'N6').or.(atpt1(i).eq.'O6'))then
              yp(1)=xt1(i)
              yp(2)=yt1(i)
              yp(3)=zt1(i)
            endif
          elseif(edge1.eq.'S'.or.edge1.eq.'s')then
            if(atpt1(i).eq.'C1*')then
              yo(1)=xt1(i)
              yo(2)=yt1(i)
              yo(3)=zt1(i)
            elseif(atpt1(i).eq.'N3')then
              yp(1)=xt1(i)
              yp(2)=yt1(i)
              yp(3)=zt1(i)
            endif
          endif
        elseif((base1.eq.'C').or.(base1.eq.'U').or.(base1.eq.'T'))then
          if(edge1.eq.'W')then
            if(atpt1(i).eq.'N3')then
              yo(1)=xt1(i)
              yo(2)=yt1(i)
              yo(3)=zt1(i)
            elseif((atpt1(i).eq.'N4').or.(atpt1(i).eq.'O4'))then
              yp(1)=xt1(i)
              yp(2)=yt1(i)
              yp(3)=zt1(i)
            endif
          elseif(edge1.eq.'H'.or.edge1.eq.'h')then
            if(atpt1(i).eq.'C5')then
              yo(1)=xt1(i)
              yo(2)=yt1(i)
              yo(3)=zt1(i)
            elseif((atpt1(i).eq.'N4').or.(atpt1(i).eq.'O4'))then
              yp(1)=xt1(i)
              yp(2)=yt1(i)
              yp(3)=zt1(i)
            endif
          elseif(edge1.eq.'S'.or.edge1.eq.'s')then
            if(atpt1(i).eq.'C1*')then
              yo(1)=xt1(i)
              yo(2)=yt1(i)
              yo(3)=zt1(i)
            elseif(atpt1(i).eq.'O2')then
              yp(1)=xt1(i)
              yp(2)=yt1(i)
              yp(3)=zt1(i)
            endif
          endif
        endif
      enddo
      ydir1(1)=yp(1)-yo(1)
      ydir1(2)=yp(2)-yo(2)
      ydir1(3)=yp(3)-yo(3)
      if(orient.eq.'T')then
        ydir1(1)=-ydir1(1)
        ydir1(2)=-ydir1(2)
        ydir1(3)=-ydir1(3)
      endif
      call normal(ydir1(1),ydir1(2),ydir1(3),ydir(1),ydir(2),ydir(3))

C Y-Axis Determination Before Transformation

      do k=1,3
        xc11(k)=xyzpr(nc11,k)
        xc12(k)=xyzpr(nc12,k)
        yaxis(k)=xc11(k)-xc12(k)
      enddo
      call normal(yaxis(1),yaxis(2),yaxis(3),eyl,eym,eyn)

C Z-Axis Determination Before Transformation

      call planeb(xcor1,nat1,np1,1,dva,el1,em1,en1,DETR,SDV4)
      call planeb(xcor2,nat2,np2,1,dva,el2,em2,en2,DETR,SDV4)
     
      ang=el1*el2+em1*em2+en1*en2
154   format(a30,3f8.3)
      if(ang.gt.0.0) then
        ezl12=el1+el2
        ezm12=em1+em2
        ezn12=en1+en2
      else
        ezl12=el1-el2
        ezm12=em1-em2
        ezn12=en1-en2
      endif
      call normal(ezl12,ezm12,ezn12,ezl,ezm,ezn)

C X-Axis Determination Before Transformation

      call cross(eyl,eym,eyn,ezl,ezm,ezn,exl1,exm1,exn1)
      call normal(exl1,exm1,exn1,exl,exm,exn)
      angle=exl*ydir(1)+exm*ydir(2)+exn*ydir(3)
      if(angle.lt.0)then
        exl=-exl
        exm=-exm
        exn=-exn
        ezl=-ezl
        ezm=-ezm
        ezn=-ezn
      endif
      call cross(exl,exm,exn,eyl,eym,eyn,ezl1,ezm1,ezn1)
      call normal(ezl1,ezm1,ezn1,ezl,ezm,ezn)
c      write(54,*)base1,base2,edge1,edge2,orient
c      write(54,*)'angle =',angle,'ang = ',ang
c      write(54,154)'Xaxis before transformation= ',exl,exm,exn
c      write(54,154)'Yaxis before transformation= ',eyl,eym,eyn
c      write(54,154)'Zaxis before transformation= ',ezl,ezm,ezn
c      write(54,*)'   '

C TRANSFORMATION

      rmat(1,1)=exl*1.0
      rmat(2,1)=eyl*1.0
      rmat(3,1)=ezl*1.0
      rmat(1,2)=exm*1.0
      rmat(2,2)=eym*1.0
      rmat(3,2)=ezm*1.0
      rmat(1,3)=exn*1.0
      rmat(2,3)=eyn*1.0
      rmat(3,3)=ezn*1.0
      do ny=1,3
        mid_c1c1(ny)=(xc11(ny)+xc12(ny))/2
      enddo
      do i=1,nat1
        xt1(i)=xt1(i)-mid_c1c1(1)
        yt1(i)=yt1(i)-mid_c1c1(2)
        zt1(i)=zt1(i)-mid_c1c1(3)
        call matmul(rmat,xt1(i),yt1(i),zt1(i),b)
        xt1(i)=b(1)
        yt1(i)=b(2)
        zt1(i)=b(3)
        xcor11(1,i)=xt1(i)
        xcor11(2,i)=yt1(i)
        xcor11(3,i)=zt1(i)
c        write(15,6)i,atpt1(i),residue(nbseg1,i),1,xt1(i),yt1(i),zt1(i)
        if(nyopt.gt.0)then
          if(atpt1(i).eq.'C1''') nc11=i 
        else
          if(((residue(nbseg1,i).eq.'ADE').or.(residue(nbseg1,i)
     1    .eq.'GUA')).and.(atpt1(i).eq.'C8')) nc11=i
          if(((residue(nbseg1,i).eq.'CYT').or.(residue(nbseg1,i)
     1    .eq.'URA').or.(residue(nbseg1,i).eq.'THY')).and.
     2    (atpt1(i).eq.'C6')) nc11=i
        endif
C
C Dirty fix of U:A W:H T getting flipped... MBU, 2015 and SINP Sept. 2016
C
       if((base2.eq.'A'.or.base2.eq.'G').and.(edge1.eq.'W'.or.edge1.eq.
     1   'H').and.(edge2.eq.'H'.or.edge2.eq.'W').and.orient.eq.'T') then
C         write(6,*) base1,base2,edge1,edge2,orient
         xt1(i)=-xt1(i)
         zt1(i)=-zt1(i)
       endif
      enddo
      do i=1,nat2
        xt2(i)=xt2(i)-mid_c1c1(1)
        yt2(i)=yt2(i)-mid_c1c1(2)
        zt2(i)=zt2(i)-mid_c1c1(3)
        call matmul(rmat,xt2(i),yt2(i),zt2(i),b)
        xt2(i)=b(1)
        yt2(i)=b(2)
        zt2(i)=b(3)
        xcor22(1,i)=xt2(i)
        xcor22(2,i)=yt2(i)
        xcor22(3,i)=zt2(i)
C        write(15,6)i,atpt2(i),residue(nbseg2,i),2,xt2(i),yt2(i),zt2(i)
        if(nyopt.gt.0)then
          if(atpt2(i).eq.'C1''') nc12=i
        else
          if(((residue(nbseg2,i).eq.'ADE').or.(residue(nbseg2,i)
     1    .eq.'GUA')).and.(atpt2(i).eq.'C8')) nc12=i
          if(((residue(nbseg2,i).eq.'CYT').or.(residue(nbseg2,i)
     1    .eq.'URA').or.(residue(nbseg2,i).eq.'THY')).and.
     2    (atpt2(i).eq.'C6')) nc12=i
        endif
C
C Dirty fix of U:A W:H T getting flipped... MBU, 2015 and SINP Sept 2016
C
       if((base2.eq.'A'.or.base2.eq.'G').and.(edge1.eq.'W'.or.edge1.eq.
     1  'H').and.(edge2.eq.'H'.or.edge2.eq.'W').and.orient.eq.'T') then
C         write(6,*) base1,base2,edge1,edge2,orient
         xt2(i)=-xt2(i)
         zt2(i)=-zt2(i)
       endif
      enddo
c      write(*,*)base1,base2,edge1,edge2,orient,'done'

C Y-axis Determination after Transformation

      do ny=1,3
c        write(*,*)'ny:  ',ny,nc11,nc12
        xc11(ny)=xcor11(ny,nc11)
        xc12(ny)=xcor22(ny,nc12)
        y_c1c1(ny)=xc11(ny)-xc12(ny)
      enddo
      call normal(y_c1c1(1),y_c1c1(2),y_c1c1(3),eyl,eym,eyn)

C Z-axis Determination after Transformation

      call planeb(xcor11,nat1,np1,0,dva,el11,em11,en11,DETR,SDV4)
      call planeb(xcor22,nat2,np2,0,dva,el22,em22,en22,DETR,SDV4)
      ang=el11*el22+em11*em22+en11*en22
      if(ang.gt.0.0)then
         ezl12=el11+el22
         ezm12=em11+em22
         ezn12=en11+en22
      else
         ezl12=el11-el22
         ezm12=em11-em22
         ezn12=en11-en22
      endif
      call normal(ezl12,ezm12,ezn12,ezl,ezm,ezn)

C X-axis Determination after Transformation

      call cross(eyl,eym,eyn,ezl,ezm,ezn,exl1,exm1,exn1)
      call normal(exl1,exm1,exn1,exl,exm,exn)
      call cross(exl,exm,exn,eyl,eym,eyn,ezl1,ezm1,ezn1)
      call normal(ezl1,ezm1,ezn1,ezl,ezm,ezn)
c      write(54,*)'ang= ',ang
c      write(54,154)'Xaxis after transformation= ',exl,exm,exn
c      write(54,154)'Yaxis after transformation= ',eyl,eym,eyn
c      write(54,154)'Zaxis after transformation= ',ezl,ezm,ezn
c      write(54,*)'   '


      close(unit=1)
1     format(a80)
2     format(23x,3f10.7)
3     format(18x,i3)
4     format(13x,a4,1x,a3,9x,3f8.3)
5     format(a1,1x,a1,1x,a1)
6     format('ATOM',2x,i5,1x,a4,1x,a3,5x,i1,4x,3f8.3)
7     format('ATOM   37     X    A             1.840   0.000  10.000')
8     format('ATOM   38     O    A             0.000   0.000  10.000')
9     format('ATOM   39     Y    A             0.000   1.840  10.000')
10    format('ATOM   40     Z    A             0.000   0.000  11.840')

      return
98    continue
        write(6,*)'Please copy DataSet.dat file at /usr/local/bin before
     1 runing the program'
        stop
      end
c********************************************************************
      SUBROUTINE NORMAL(A,B,C,X,Y,Z)
      DOUBLE PRECISION A,B,C, X,Y,Z, ALEN,ALENGT
C
      ALEN =  A * A + B * B + C * C
C
      ALENGT = DSQRT(ALEN)
C
      IF (ABS(ALENGT) .LT. 1.0E-6) THEN
         X=0.D0
         Y=0.D0
         Z=0.D0
      ELSE
         X = A/ALENGT
         Y = B/ALENGT
         Z = C/ALENGT
      ENDIF
C
      RETURN
      END

C***********************************************************************
      SUBROUTINE PLANEB(X,NTOM,NX,NPLN,DV,EL,EM,EN,DET,SDV)
        parameter ( nrs = 9999 )
      DOUBLE PRECISION EL,EM,EN,Sl,Sm,Sn,y12(3),y13(3)
      DIMENSION X(3,nrs),NX(nrs),DV(nrs),Y(3,nrs),T(3,3),B(3),A(3)
      DIMENSION XY(3),XZ(3),VEC1(3),VEC2(3)
      external matmul
      CHARACTER*1 ansnrm
      INTEGER XY,XZ
      ansnrm='Y'

      DMOVE = 0.0
      ntim = 0
      NCOUNT = 0
      DO 100 J1=1,NTOM
        NX1=NX(J1)
        IF(NX1.NE.0) THEN
          NCOUNT = NCOUNT + 1
          DO 2000 K=1,3
            Y(K,NCOUNT)=X(K,NX1)
 2000     CONTINUE
        END IF
 100  CONTINUE
      if(ncount.ge.4) then
        NTOM1 = NCOUNT
      else
        write(40,*) ncount,ntom,nx1,' nx=',nx
        write(6,801) 'in planeb',ncount
        stop 1
      endif
      if(ansnrm.ne.'Y') then
        do while(ntim.le.5)
          DO 101 J1=1,3
            B(J1)=0.0
            DO 101 J2=1,3
              T(J1,J2)=0.0
 101      CONTINUE
          DO 102 J1=1,NTOM1
            T(1,1)=T(1,1)+Y(1,J1)*Y(1,J1)
            T(1,2)=T(1,2)+Y(1,J1)*Y(2,J1)
            T(1,3)=T(1,3)+Y(1,J1)*Y(3,J1)
            T(2,2)=T(2,2)+Y(2,J1)*Y(2,J1)
            T(2,3)=T(2,3)+Y(2,J1)*Y(3,J1)
            T(3,3)=T(3,3)+Y(3,J1)*Y(3,J1)
            B(1)=B(1)+Y(1,J1)
            B(2)=B(2)+Y(2,J1)
            B(3)=B(3)+Y(3,J1)
 102      CONTINUE
          DO 103 J1=1,3
            DO 103 J2=J1,3
              T(J2,J1)=T(J1,J2)
 103      CONTINUE
          CALL MINV(T,3,DX,XY,XZ)
          IF(DX.NE.0.0E0) THEN
            B1 = B(1)
            B2 = B(2)
            B3 = B(3)
            CALL MATMUL(T,B1,B2,B3,A)
            DET=1.0/(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))
            DET=SQRT(DET)
            EL=A(1)*DET
            EM=A(2)*DET
            EN=A(3)*DET
            SDV = 0.0
            DO 104 J1=1,NTOM1
              DV(J1)=EL*Y(1,J1)+EM*Y(2,J1)+EN*Y(3,J1)-DET
              SDV = SDV + DV(J1) * DV(J1)
 104        CONTINUE
            SDV = SQRT(SDV)
            IF(NPLN.NE.0.AND.SDV.GT.0.75) THEN
              DMOVE = DMOVE + 25.0
              ntim = ntim + 1
              do j1=1,ntom1
                do k=1,3
                  y(k,j1)=y(k,j1)+dmove
                enddo
              enddo
            else
              return
            endif
          endif
        enddo

      else

C LEAST SQUARE FIT did not work for the set of coordinates.  Going back to
C finding best plane through cross-product method and averaging.

           do k=1,3
              y12(k)=y(k,1) - y(k,2)
              y13(k)=y(k,1) - y(k,3)
           enddo

           call cross(y12(1),y12(2),y12(3),y13(1),y13(2),y13(3),
     1 sl,sm,sn)
           do 105 i=2,ntom1-2
             do 106 k=1,3
               y12(k) = y(k,i) - y(k,i+1)
               y13(k) = y(k,i) - y(k,i+2)
 106         continue
             call cross(y12(1),y12(2),y12(3),y13(1),y13(2),y13(3),
     1  el,em,en)
             ang = (sl*el+sm*em+sn*en)/(sqrt(sl*sl+sm*sm+sn*sn)*
     1  sqrt(el*el+em*em+en*en))
             IF(ang.lt.0.0) then
                el = -el
                em = -em
                en = -en
             endif
             sl = sl+el
             sm = sm+em
             sn = sn+en
             call normal(sl,sm,sn,el,em,en)
             sl=el
             sm=em
             sn=en
 105       continue

          endif
      RETURN
801    format('BAD Selection of BASE, only',I3,' atoms found in PLANEB')
      END

c********************************************************************
      SUBROUTINE CROSS(A,B,C,A1,B1,C1,ALPHA,BETA,GAMA)

      DOUBLE PRECISION A,B,C,ALPHA,BETA,GAMA, A1,B1,C1
      ALPHA=B*C1-C*B1
      BETA=C*A1-A*C1
      GAMA=A*B1-B*A1
      AMOD=DSQRT(ALPHA*ALPHA+BETA*BETA+GAMA*GAMA)
      ALPHA=ALPHA/AMOD
      BETA=BETA/AMOD
      GAMA=GAMA/AMOD

      RETURN
      END

C***********************************************************************

      SUBROUTINE MINV (A,N,D,L,M)

C     ----- STANDARD IBM MATRIX INVERSION ROUTINE -----

      DIMENSION A(1),L(1),M(1)

C     ----- SEARCH FOR LARGEST ELEMENT -----

      D = 1.E0
      NK = -N
      DO 80 K = 1,N
         NK = NK+N
         L(K) = K
         M(K) = K
         KK = NK+K
         BIGA = A(KK)
         DO 20 J = K,N
         IZ = N*(J-1)
            DO 20 I = K,N
              IJ = IZ+I
   10         IF( ABS(BIGA)- ABS(A(IJ))) 15,20,20
   15         BIGA = A(IJ)
              L(K) = I
              M(K) = J
   20    CONTINUE

C     ----- INTERCHANGE ROWS -----

         J = L(K)
         IF(J-K) 35,35,25
   25    KI = K-N
         DO 30 I = 1,N
           KI = KI+N
           HOLD = -A(KI)
           JI = KI-K+J
           A(KI) = A(JI)
   30    A(JI) = HOLD

C     ----- INTERCHANGE COLUMNS -----

   35     I = M(K)
          IF(I-K) 45,45,38
   38     JP = N*(I-1)
          DO 40 J = 1,N
            JK = NK+J
            JI = JP+J
        HOLD = -A(JK)
          A(JK) = A(JI)
   40     A(JI) = HOLD

C     ----- DIVIDE COLUMN BY MINUS PIVOT -----

   45     IF(BIGA) 48,46,48
   46     D = 0.E0
          GO TO 150
   48     DO 55 I = 1,N
            IF(I-K) 50,55,50
   50       IK = NK+I
            A(IK) = A(IK)/(-BIGA)
   55     CONTINUE

C     ----- REDUCE MATRIX -----

          DO 65 I = 1,N
             IK = NK+I
             HOLD = A(IK)
        IJ = I-N
          DO 65 J = 1,N
           IJ = IJ+N
          IF(I-K) 60,65,60
   60      IF(J-K) 62,65,62
   62      KJ = IJ-I+K
             A(IJ) = HOLD*A(KJ)+A(IJ)
   65    CONTINUE

C     ----- DIVIDE ROW BY PIVOT -----

         KJ = K-N
         DO 75 J = 1,N
           KJ = KJ+N
           IF(J-K) 70,75,70
   70   A(KJ) = A(KJ)/BIGA
   75    CONTINUE

C     ----- PRODUCT OF PIVOTS -----

         D = D*BIGA

C     ----- REPLACE PIVOT BY RECIPROCAL -----

         A(KK) = 1.E0/BIGA
   80 CONTINUE

C     ----- FINAL ROW AND COLUMN INTERCHANGE -----

      K = N
  100 K = (K-1)
      IF(K) 150,150,105
  105 I = L(K)
      IF(I-K) 120,120,108
  108 JQ = N*(K-1)
      JR = N*(I-1)
      DO 110 J = 1,N
         JK = JQ+J
         HOLD = A(JK)
         JI = JR+J
         A(JK) = -A(JI)
  110 A(JI) = HOLD
  120 J = M(K)
      IF(J-K) 100,100,125
  125 KI = K-N
      DO 130 I = 1,N
         KI = KI+N
         HOLD = A(KI)
         JI = KI-K+J
         A(KI) = -A(JI)
  130 A(JI) = HOLD
      GOTO 100
  150 RETURN
      END

      SUBROUTINE MATMUL(R,A1,A2,A3,B)

C Matrix multiplication subroutine.

      DIMENSION B(3),R(3,3)
      DO 923 I=1,3
         B(I) = R(I,1) * A1 + R(I,2) * A2 + R(I,3) * A3
 923    CONTINUE
      RETURN
      END

