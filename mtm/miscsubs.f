c
c $Id: miscsubs.f,v 1.1 1997/01/10 01:14:14 weibel Exp $
c
c     calculate cos function by recursion relations for efficiency
c
      SUBROUTINE COSSIN(N,C,S,C0,S0,F,DT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(1),S(1)
      DATA PI/3.14159265358979D0/
      C(1)=C0
      S(1)=S0
      CS=DCOS(2.D0*PI*F*DT)
      SN=DSIN(2.D0*PI*F*DT)
      DO 100 I=2,N
      C(I)=C(I-1)*CS-S(I-1)*SN
  100 S(I)=C(I-1)*SN+S(I-1)*CS
      RETURN
      END
c
      SUBROUTINE FFT2(ar,ai,N)
c
c     Jeff Park's fft
c
c     2 real input arrays rather than complex 
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension ar(1),ai(1)
      mex=dlog(dble(float(n)))/.693147d0
      nv2=N/2
      nm1=N-1
      j=1
      do 7 i=1,nm1
      if(i .ge. j) go to 5
      tr=ar(j)
      ti=ai(j)
      ar(j)=ar(i)
      ai(j)=ai(i)
      ar(i)=tr
      ai(i)=ti
   5  k=nv2
   6  if(k .ge. j) go to 7
      j=j-k
      k=k/2
      go to 6
   7  j=j+k
      pi=3.14159265358979d0
      do 20 l=1,mex
      le=2**l
      le1=le/2
      ur=1.0
      ui=0.
      wr=dcos(pi/le1 )
      wi=dsin (pi/le1)
      do 20 j=1,le1
      do 10 i=j,N,le
      ip=i+le1
      tr=ar(ip)*ur - ai(ip)*ui
      ti=ai(ip)*ur + ar(ip)*ui
      ar(ip)=ar(i)-tr
      ai(ip)=ai(i) - ti
      ar(i)=ar(i)+tr
      ai(i)=ai(i)+ti
  10  continue
      utemp=ur
      ur=ur*wr - ui*wi
      ui=ui*wr + utemp*wi
  20  continue
      return
      end
c
c     subroutines TRIDIB, TINVIT from DEISPK
c  

      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(M)
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0D0
      DO 40 I = 1, N
         X1 = U
         U = 0.0D0
         IF (I .NE. N) U = DABS(E(I+1))
         XU = DMIN1(D(I)-(X1+U),XU)
         X0 = DMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         TST1 = DABS(D(I)) + DABS(D(I-1))
         TST2 = TST1 + DABS(E(I))
         IF (TST2 .GT. TST1) GO TO 40
   20    E2(I) = 0.0D0
   40 CONTINUE
      X1 = N
      X1 = X1 * EPSLON(DMAX1(DABS(XU),DABS(X0)))
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5D0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0D0
      DO 120 Q = P, N
         X1 = U
         U = 0.0D0
         V = 0.0D0
         IF (Q .EQ. N) GO TO 110
         U = DABS(E(Q+1))
         V = E2(Q+1)
  110    XU = DMIN1(D(Q)-(X1+U),XU)
         X0 = DMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0D0) GO TO 140
  120 CONTINUE
  140 X1 = EPSLON(DMAX1(DABS(XU),DABS(X0)))
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
      X0 = UB
      ISTURM = 5
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
      K = M2
  250    XU = LB
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
  300    X1 = (XU + X0) * 0.5D0
         IF ((X0 - XU) .LE. DABS(EPS1)) GO TO 420
         TST1 = 2.0D0 * (DABS(XU) + DABS(X0))
         TST2 = TST1 + (X0 - XU)
         IF (TST2 .EQ. TST1) GO TO 420
  320    S = P - 1
         U = 1.0D0
         DO 340 I = P, Q
            IF (U .NE. 0.0D0) GO TO 325
            V = DABS(E(I)) / EPSLON(1.0D0)
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0D0) S = S + 1
  340    CONTINUE
         GO TO (60,80,200,220,360), ISTURM
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
      END
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,EPSLON,
     X       PYTHAG
      INTEGER IND(M)
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0D0 - E2(1)
      Q = 0
  100 P = Q + 1
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140
  120 CONTINUE
  140 TAG = TAG + 1
      S = 0
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
         XU = 1.0D0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0D0
         GO TO 870
  490    NORM = DABS(D(P))
         IP = P + 1
         DO 500 I = IP, Q
  500    NORM = DMAX1(NORM, DABS(D(I))+DABS(E(I)))
         EPS2 = 1.0D-3 * NORM
         EPS3 = EPSLON(NORM)
         UK = Q - P + 1
         EPS4 = UK * EPS3
         UK = EPS4 / DSQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3
  520    V = 0.0D0
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0D0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0D0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
         IF (U .EQ. 0.0D0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0D0
         RV3(Q) = 0.0D0
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
         IF (GROUP .EQ. 0) GO TO 700
         J = R
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0D0
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
  680    CONTINUE
  700    NORM = 0.0D0
         DO 720 I = P, Q
  720    NORM = NORM + DABS(RV6(I))
         IF (NORM .GE. 1.0D0) GO TO 840
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0D0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
  780    DO 820 I = IP, Q
            U = RV6(I)
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
         ITS = ITS + 1
         GO TO 600
  830    IERR = -R
         XU = 0.0D0
         GO TO 870
  840    U = 0.0D0
         DO 860 I = P, Q
  860    U = PYTHAG(U,RV6(I))
         XU = 1.0D0 / U
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0D0
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
         X0 = X1
  920 CONTINUE
      IF (Q .LT. N) GO TO 100
 1001 RETURN
      END
      DOUBLE PRECISION FUNCTION EPSLON (X)
      DOUBLE PRECISION X
      DOUBLE PRECISION A,B,C,EPS
      A = 4.0D0/3.0D0
   10 B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS .EQ. 0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
c================================================================
c   routines from zlinpack that are used in envelope estimation
c================================================================

      SUBROUTINE ZQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)                   
      INTEGER LDX,N,P,JOB                                               
      INTEGER JPVT(1)                                                   
      COMPLEX*16 X(LDX,1),QRAUX(1),WORK(1)                              
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU                                 
      DOUBLE PRECISION MAXNRM,DZNRM2,TT                                 
      COMPLEX*16 ZDOTC,NRMXL,T                                          
      LOGICAL NEGJ,SWAPJ                                                
      COMPLEX*16 CSIGN,ZDUM,ZDUM1,ZDUM2                                 
      DOUBLE PRECISION CABS1                                            
      DOUBLE PRECISION DREAL,DIMAG                                      
      COMPLEX*16 ZDUMR,ZDUMI                                            
      DREAL(ZDUMR) = ZDUMR                                              
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI                               
      CSIGN(ZDUM1,ZDUM2) = CDABS(ZDUM1)*(ZDUM2/CDABS(ZDUM2))            
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))               
      PL = 1                                                            
      PU = 0                                                            
      IF (JOB .EQ. 0) GO TO 60                                          
         DO 20 J = 1, P                                                 
            SWAPJ = JPVT(J) .GT. 0                                      
            NEGJ = JPVT(J) .LT. 0                                       
            JPVT(J) = J                                                 
            IF (NEGJ) JPVT(J) = -J                                      
            IF (.NOT.SWAPJ) GO TO 10                                    
               IF (J .NE. PL) CALL ZSWAP(N,X(1,PL),1,X(1,J),1)          
               JPVT(J) = JPVT(PL)                                       
               JPVT(PL) = J                                             
               PL = PL + 1                                              
   10       CONTINUE                                                    
   20    CONTINUE                                                       
         PU = P                                                         
         DO 50 JJ = 1, P                                                
            J = P - JJ + 1                                              
            IF (JPVT(J) .GE. 0) GO TO 40                                
               JPVT(J) = -JPVT(J)                                       
               IF (J .EQ. PU) GO TO 30                                  
                  CALL ZSWAP(N,X(1,PU),1,X(1,J),1)                      
                  JP = JPVT(PU)                                         
                  JPVT(PU) = JPVT(J)                                    
                  JPVT(J) = JP                                          
   30          CONTINUE                                                 
               PU = PU - 1                                              
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
      IF (PU .LT. PL) GO TO 80                                          
      DO 70 J = PL, PU                                                  
         QRAUX(J) = DCMPLX(DZNRM2(N,X(1,J),1),0.0D0)                    
         WORK(J) = QRAUX(J)                                             
   70 CONTINUE                                                          
   80 CONTINUE                                                          
      LUP = MIN0(N,P)                                                   
      DO 200 L = 1, LUP                                                 
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120                        
            MAXNRM = 0.0D0                                              
            MAXJ = L                                                    
            DO 100 J = L, PU                                            
               IF (DREAL(QRAUX(J)) .LE. MAXNRM) GO TO 90                
                  MAXNRM = DREAL(QRAUX(J))                              
                  MAXJ = J                                              
   90          CONTINUE                                                 
  100       CONTINUE                                                    
            IF (MAXJ .EQ. L) GO TO 110                                  
               CALL ZSWAP(N,X(1,L),1,X(1,MAXJ),1)                       
               QRAUX(MAXJ) = QRAUX(L)                                   
               WORK(MAXJ) = WORK(L)                                     
               JP = JPVT(MAXJ)                                          
               JPVT(MAXJ) = JPVT(L)                                     
               JPVT(L) = JP                                             
  110       CONTINUE                                                    
  120    CONTINUE                                                       
         QRAUX(L) = (0.0D0,0.0D0)                                       
         IF (L .EQ. N) GO TO 190                                        
            NRMXL = DCMPLX(DZNRM2(N-L+1,X(L,L),1),0.0D0)                
            IF (CABS1(NRMXL) .EQ. 0.0D0) GO TO 180                      
               IF (CABS1(X(L,L)) .NE. 0.0D0)                            
     *            NRMXL = CSIGN(NRMXL,X(L,L))                           
               CALL ZSCAL(N-L+1,(1.0D0,0.0D0)/NRMXL,X(L,L),1)           
               X(L,L) = (1.0D0,0.0D0) + X(L,L)                          
               LP1 = L + 1                                              
               IF (P .LT. LP1) GO TO 170                                
               DO 160 J = LP1, P                                        
                  T = -ZDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)            
                  CALL ZAXPY(N-L+1,T,X(L,L),1,X(L,J),1)                 
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150               
                  IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 150             
                     TT = 1.0D0 - (CDABS(X(L,J))/DREAL(QRAUX(J)))**2    
                     TT = DMAX1(TT,0.0D0)                               
                     T = DCMPLX(TT,0.0D0)                               
                     TT = 1.0D0                                         
     *                    + 0.05D0*TT                                   
     *                      *(DREAL(QRAUX(J))/DREAL(WORK(J)))**2        
                     IF (TT .EQ. 1.0D0) GO TO 130                       
                        QRAUX(J) = QRAUX(J)*CDSQRT(T)                   
                     GO TO 140                                          
  130                CONTINUE                                           
                        QRAUX(J) = DCMPLX(DZNRM2(N-L,X(L+1,J),1),0.0D0) 
                        WORK(J) = QRAUX(J)                              
  140                CONTINUE                                           
  150             CONTINUE                                              
  160          CONTINUE                                                 
  170          CONTINUE                                                 
               QRAUX(L) = X(L,L)                                        
               X(L,L) = -NRMXL                                          
  180       CONTINUE                                                    
  190    CONTINUE                                                       
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE ZQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)      
      INTEGER LDX,N,K,JOB,INFO                                          
      COMPLEX*16 X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1),XB(1)  
      INTEGER I,J,JJ,JU,KP1                                             
      COMPLEX*16 ZDOTC,T,TEMP                                           
      LOGICAL CB,CQY,CQTY,CR,CXB                                        
      COMPLEX*16 ZDUM                                                   
      DOUBLE PRECISION CABS1                                            
      DOUBLE PRECISION DREAL,DIMAG                                      
      COMPLEX*16 ZDUMR,ZDUMI                                            
      DREAL(ZDUMR) = ZDUMR                                              
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI                               
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))               
      INFO = 0                                                          
      CQY = JOB/10000 .NE. 0                                            
      CQTY = MOD(JOB,10000) .NE. 0                                      
      CB = MOD(JOB,1000)/100 .NE. 0                                     
      CR = MOD(JOB,100)/10 .NE. 0                                       
      CXB = MOD(JOB,10) .NE. 0                                          
      JU = MIN0(K,N-1)                                                  
      IF (JU .NE. 0) GO TO 40                                           
         IF (CQY) QY(1) = Y(1)                                          
         IF (CQTY) QTY(1) = Y(1)                                        
         IF (CXB) XB(1) = Y(1)                                          
         IF (.NOT.CB) GO TO 30                                          
            IF (CABS1(X(1,1)) .NE. 0.0D0) GO TO 10                      
               INFO = 1                                                 
            GO TO 20                                                    
   10       CONTINUE                                                    
               B(1) = Y(1)/X(1,1)                                       
   20       CONTINUE                                                    
   30    CONTINUE                                                       
         IF (CR) RSD(1) = (0.0D0,0.0D0)                                 
      GO TO 250                                                         
   40 CONTINUE                                                          
         IF (CQY) CALL ZCOPY(N,Y,1,QY,1)                                
         IF (CQTY) CALL ZCOPY(N,Y,1,QTY,1)                              
         IF (.NOT.CQY) GO TO 70                                         
            DO 60 JJ = 1, JU                                            
               J = JU - JJ + 1                                          
               IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 50                 
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  T = -ZDOTC(N-J+1,X(J,J),1,QY(J),1)/X(J,J)             
                  CALL ZAXPY(N-J+1,T,X(J,J),1,QY(J),1)                  
                  X(J,J) = TEMP                                         
   50          CONTINUE                                                 
   60       CONTINUE                                                    
   70    CONTINUE                                                       
         IF (.NOT.CQTY) GO TO 100                                       
            DO 90 J = 1, JU                                             
               IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 80                 
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  T = -ZDOTC(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)            
                  CALL ZAXPY(N-J+1,T,X(J,J),1,QTY(J),1)                 
                  X(J,J) = TEMP                                         
   80          CONTINUE                                                 
   90       CONTINUE                                                    
  100    CONTINUE                                                       
         IF (CB) CALL ZCOPY(K,QTY,1,B,1)                                
         KP1 = K + 1                                                    
         IF (CXB) CALL ZCOPY(K,QTY,1,XB,1)                              
         IF (CR .AND. K .LT. N) CALL ZCOPY(N-K,QTY(KP1),1,RSD(KP1),1)   
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120                        
            DO 110 I = KP1, N                                           
               XB(I) = (0.0D0,0.0D0)                                    
  110       CONTINUE                                                    
  120    CONTINUE                                                       
         IF (.NOT.CR) GO TO 140                                         
            DO 130 I = 1, K                                             
               RSD(I) = (0.0D0,0.0D0)                                   
  130       CONTINUE                                                    
  140    CONTINUE                                                       
         IF (.NOT.CB) GO TO 190                                         
            DO 170 JJ = 1, K                                            
               J = K - JJ + 1                                           
               IF (CABS1(X(J,J)) .NE. 0.0D0) GO TO 150                  
                  INFO = J                                              
                  GO TO 180                                             
  150          CONTINUE                                                 
               B(J) = B(J)/X(J,J)                                       
               IF (J .EQ. 1) GO TO 160                                  
                  T = -B(J)                                             
                  CALL ZAXPY(J-1,T,X(1,J),1,B,1)                        
  160          CONTINUE                                                 
  170       CONTINUE                                                    
  180       CONTINUE                                                    
  190    CONTINUE                                                       
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240                          
            DO 230 JJ = 1, JU                                           
               J = JU - JJ + 1                                          
               IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 220                
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  IF (.NOT.CR) GO TO 200                                
                     T = -ZDOTC(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)         
                     CALL ZAXPY(N-J+1,T,X(J,J),1,RSD(J),1)              
  200             CONTINUE                                              
                  IF (.NOT.CXB) GO TO 210                               
                     T = -ZDOTC(N-J+1,X(J,J),1,XB(J),1)/X(J,J)          
                     CALL ZAXPY(N-J+1,T,X(J,J),1,XB(J),1)               
  210             CONTINUE                                              
                  X(J,J) = TEMP                                         
  220          CONTINUE                                                 
  230       CONTINUE                                                    
  240    CONTINUE                                                       
  250 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)                            
      COMPLEX*16 ZX(1),ZY(1),ZA                                         
      DOUBLE PRECISION DCABS1                                           
      IF(N.LE.0)RETURN                                                  
      IF (DCABS1(ZA) .EQ. 0.0D0) RETURN                                 
      IF (INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                              
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        ZY(IY) = ZY(IY) + ZA*ZX(IX)                                     
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
   20 DO 30 I = 1,N                                                     
        ZY(I) = ZY(I) + ZA*ZX(I)                                        
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  ZCOPY(N,ZX,INCX,ZY,INCY)                              
      COMPLEX*16 ZX(1),ZY(1)                                            
      INTEGER I,INCX,INCY,IX,IY,N                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        ZY(IY) = ZX(IX)                                                 
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
   20 DO 30 I = 1,N                                                     
        ZY(I) = ZX(I)                                                   
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  ZSCAL(N,ZA,ZX,INCX)                                   
      COMPLEX*16 ZA,ZX(1)                                               
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
      IX = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      DO 10 I = 1,N                                                     
        ZX(IX) = ZA*ZX(IX)                                              
        IX = IX + INCX                                                  
   10 CONTINUE                                                          
      RETURN                                                            
   20 DO 30 I = 1,N                                                     
        ZX(I) = ZA*ZX(I)                                                
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  ZSWAP (N,ZX,INCX,ZY,INCY)                             
      COMPLEX*16 ZX(1),ZY(1),ZTEMP                                      
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        ZTEMP = ZX(IX)                                                  
        ZX(IX) = ZY(IY)                                                 
        ZY(IY) = ZTEMP                                                  
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
   20 DO 30 I = 1,N                                                     
        ZTEMP = ZX(I)                                                   
        ZX(I) = ZY(I)                                                   
        ZY(I) = ZTEMP                                                   
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      DOUBLE PRECISION FUNCTION DZNRM2( N, ZX, INCX)                    
      LOGICAL IMAG, SCALE                                               
      INTEGER          NEXT                                             
      DOUBLE PRECISION CUTLO, CUTHI, HITEST, SUM, XMAX, ABSX, ZERO, ONE 
      COMPLEX*16      ZX(1)                                             
      DOUBLE PRECISION DREAL,DIMAG                                      
      COMPLEX*16 ZDUMR,ZDUMI                                            
      DREAL(ZDUMR) = ZDUMR                                              
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI                               
      DATA         ZERO, ONE /0.0D0, 1.0D0/                             
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /                        
      IF(N .GT. 0) GO TO 10                                             
         DZNRM2  = ZERO                                                 
         GO TO 300                                                      
   10 ASSIGN 30 TO NEXT                                                 
      SUM = ZERO                                                        
      NN = N * INCX                                                     
      DO 210 I=1,NN,INCX                                                
         ABSX = DABS(DREAL(ZX(I)))                                      
         IMAG = .FALSE.                                                 
         GO TO NEXT,(30, 50, 70, 90, 110)                               
   30 IF( ABSX .GT. CUTLO) GO TO 85                                     
      ASSIGN 50 TO NEXT                                                 
      SCALE = .FALSE.                                                   
   50 IF( ABSX .EQ. ZERO) GO TO 200                                     
      IF( ABSX .GT. CUTLO) GO TO 85                                     
      ASSIGN 70 TO NEXT                                                 
      GO TO 105                                                         
  100 ASSIGN 110 TO NEXT                                                
      SUM = (SUM / ABSX) / ABSX                                         
  105 SCALE = .TRUE.                                                    
      XMAX = ABSX                                                       
      GO TO 115                                                         
   70 IF( ABSX .GT. CUTLO ) GO TO 75                                    
  110 IF( ABSX .LE. XMAX ) GO TO 115                                    
         SUM = ONE + SUM * (XMAX / ABSX)**2                             
         XMAX = ABSX                                                    
         GO TO 200                                                      
  115 SUM = SUM + (ABSX/XMAX)**2                                        
      GO TO 200                                                         
   75 SUM = (SUM * XMAX) * XMAX                                         
   85 ASSIGN 90 TO NEXT                                                 
      SCALE = .FALSE.                                                   
      HITEST = CUTHI/FLOAT( N )                                         
   90 IF(ABSX .GE. HITEST) GO TO 100                                    
         SUM = SUM + ABSX**2                                            
  200 CONTINUE                                                          
      IF(IMAG) GO TO 210                                                
         ABSX = DABS(DIMAG(ZX(I)))                                      
         IMAG = .TRUE.                                                  
      GO TO NEXT,(  50, 70, 90, 110 )                                   
  210 CONTINUE                                                          
      DZNRM2 = DSQRT(SUM)                                               
      IF(SCALE) DZNRM2 = DZNRM2 * XMAX                                  
  300 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      COMPLEX FUNCTION ZDOTC*16(N,ZX,INCX,ZY,INCY)                      
      COMPLEX*16 ZX(1),ZY(1),ZTEMP                                      
      ZTEMP = (0.0D0,0.0D0)                                             
      ZDOTC = (0.0D0,0.0D0)                                             
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)                           
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      ZDOTC = ZTEMP                                                     
      RETURN                                                            
   20 DO 30 I = 1,N                                                     
        ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)                             
   30 CONTINUE                                                          
      ZDOTC = ZTEMP                                                     
      RETURN                                                            
      END                                                               
      DOUBLE PRECISION FUNCTION DCABS1(Z)                               
      COMPLEX*16 Z,ZZ                                                   
      DOUBLE PRECISION T(2)                                             
      EQUIVALENCE (ZZ,T(1))                                             
      ZZ = Z                                                            
      DCABS1 = DABS(T(1)) + DABS(T(2))                                  
      RETURN                                                            
      END                                                               
