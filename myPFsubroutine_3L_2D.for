C ======================================================================
C ABAQUS user subroutine for phase field fracture
C Three layered system:
C     1st - Phase-field element - UEL
C     2nd - Plane-strain element - UEL
C     3rd - Visualization - UMAT
C ======================================================================
C 1st layer: User Subroutine UEL for phase-field element
C     Type 1: C2P4 phase-field rectangular element
C     Type 2: C2P3 phase-field triangle element
C ======================================================================
C 2nd layer: User subroutine UEL for plane-strain element
C     Type 21: C2D4 isotropic energy degradation
C     Type 22: C2D4 anisotropic energy degradation by
C                 spectral decomposition
C     Type 23: C2D4 anisotropic energy degradation by
C                 volumetric-deviatoric decomposition
C
C     Type 51: C2D3 isotropic energy degradation
C     Type 52: C2D3 anisotropic energy degradation by
C                 spectral decomposition
C     Type 53: C2D3 anisotropic energy degradation by
C                 volumetric-deviatoric decomposition
C ======================================================================
C Material properties to be given through the input file (*.inp), are
C
C 1st layer: Phase-field element
C PROPS(1) = Crack scale parameter (lc)
C PROPS(2) = Crack surface energy (Gc)
C PROPS(3) = thickness of the element (t)
C
C 2nd layer: Plane-strain element
C for brittle fracture (Type 21-23 & 51-53):
C PROPS(1) = Young's modulus (E)
C PROPS(2) = Poisson's ratio (nu)
C PROPS(3) = thickness of the element (t)
C PROPS(4) = Parameter k (stabilization of the stiffness matrix)
C
C ======================================================================
C ---- Used variables ---------------------
C N_ELEM - number of elements used in the model divided
C            by 3 - (N_phase+N_stress+N_UMAT)/3 (to be changed for each model)
C NSTVTO - solution dependent variables for the phase field element
C            (phase, energy history)
C NSTVTT - solution dependent variables for the strain element
C            (displacements, strains, stresses, elastic stresses,
C             energies, phase)
C NSTV - overall solution dependent variables (NSTVTO+NSTVTT+2), where
C           the additional 2 variables are the: time and iteration number
C ======================================================================
C ---- COMMON blocks ---------------------
C USRVAR(N_ELEM,NSTV,4) - user defined variables at each integration pt
C 2nd layer: Dsiplacement (stress-strain) element
C for Type 21 & 51:
C     SDV1-SDV2     Displacements: u_x, u_y
C     SDV3-SDV4     Axial strains: epsilon_x, epsilon_y
C     SDV5          Engineering shear strains: gamma_xy
C     SDV6-SDV7     Axial stresses: sigma_x, sigma_y
C     SDV8          Shear stresses: sigma_xy
C     SDV9-SDV10    Undegenerated axial stresses: sigma_x0, sigma_y0
C     SDV11         Undegenerated shear stresses: sigma_xy0
C     SDV12         Strain energy: psi
C     SDV13         Undegenerated strain energy: psi_0
C     SDV14         Phase-field: d
C for Type 22 & 52:
C     SDV1-SDV2     Displacements: u_x, u_y
C     SDV3-SDV4     Axial strains: epsilon_x, epsilon_y
C     SDV5          Engineering shear strains: gamma_xy
C     SDV6-SDV7     Axial stresses: sigma_x, sigma_y
C     SDV8          Shear stresses: sigma_xy
C     SDV9-SDV10    Principal strains: epsilon_1, epsilon_2
C     SDV11         1st invariant of strain: volumetric strain
C     SDV12         Strain energy: psi
C     SDV13         Undegenerated strain energy: psi_0^+
C     SDV14         Phase-field: d
C for Type 23 & 53:
C     SDV1-SDV2     Displacements: u_x, u_y
C     SDV3-SDV4     Axial strains: epsilon_x, epsilon_y
C     SDV5          Engineering shear strains: gamma_xy
C     SDV6-SDV7     Axial stresses: sigma_x, sigma_y
C     SDV8          Shear stresses: sigma_xy
C     SDV9-SDV11    Deviatoric strains
C     SDV12         Strain energy: psi
C     SDV13         Undegenerated strain energy: psi_0^+
C     SDV14         Phase-field: d
C
C 1st layer: Phase-field element
C     SDV15         Phase-field: d
C     SDV16         History field: H
C     SDV17         Time
C     SDV18         Iteration number
C ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3 NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4 PERIOD)
C     ==================================================================
      INCLUDE 'ABA_PARAM.INC' ! implicit real(a-h o-z)
C     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-8,FOUR=4.D0,RP25 = 0.25D0,HALF=0.5D0,SIX=6.D0,
     2 N_ELEM=4434,NSTVTO=2,NSTVTT=14,NSTV=18)
C     ==================================================================
C     Initialization for all the element types
C     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1 SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2 U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3 PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5 JPROPS(*)

      INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY,IZ

      REAL*8 AINTW(4) ! weights of integration points
      REAL*8 XII(4,2) ! local coordinates of all integration points
      REAL*8 XI(2) ! local coordinates of a integration point
      REAL*8 dNdxi(NNODE,2) ! local derivatives of the shape functions
      REAL*8 VJACOB(2,2) ! Jacobian matrix of coordinates
      REAL*8 dNdx(NNODE,2) ! global derivatives of the shape functions
      REAL*8 VJABOBINV(2,2) ! inverse of vjacob
      REAL*8 AN(4) ! shape function of each node
      REAL*8 BP(2,NDOFEL) ! B matrix of the phase-field
      REAL*8 DP(2) ! Gradient of phase-field variable
      REAL*8 SDV(NSTV) ! state dependent variables
      REAL*8 BB(3,NDOFEL) ! B matrix of the displacement
      REAL*8 CMAT(3,3) ! materials stiffness matrix
      REAL*8 EPS(3) ! strains
      REAL*8 STRESS(3) ! stresses
      REAL*8 VNI(2,NDOFEL) ! shape function of the displacement
      REAL*8 ULOC(2) ! displacements

      REAL*8 DTM ! determinant of Jacobian matrix
      REAL*8 THCK ! thickness of the element (t)
      REAL*8 HIST ! elastic strain energy history
      REAL*8 CLPAR ! Crack scale parameter (lc)
      REAL*8 GCPAR ! Crack surface energy (Gc)
      REAL*8 EMOD,ENU,PARK,ENG
      REAL*8 EBULK3,EBULK,EG,EG2,EG3,ELAM

C     for type 22
      REAL*8 EIGV(2) ! The three principal strains
      REAL*8 ALPHAI(2) ! alphai = 1 for principal strain >= 0
      REAL*8 DEIGDEPS(2,3)
      REAL*8 DEIG2DEPS(3,3,2)
      REAL*8 CLMAT(2,2) ! L = principal sigma / principal eps
      REAL*8 EPSC(3)

      COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT

C     ==================================================================
C     ******************************************************************
C     Constructing Phase-field element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2P4 phase-field rectangular element
C     ------------------------------------------------------------------
      IF (JTYPE.EQ.ONE) THEN ! C2P4 phase-field rectangular element
C     ==================================================================
C     Time an iteration variables
C     ==================================================================
        TIMEZ=USRVAR(JELEM,17,1)
        IF (TIMEZ.LT.TIME(2)) THEN
          USRVAR(JELEM,17,1)=TIME(2)
          USRVAR(JELEM,18,1)=ZERO
        ELSE
          USRVAR(JELEM,18,1)=USRVAR(JELEM,18,1)+ONE
        ENDIF
        STEPITER=USRVAR(JELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        CLPAR=PROPS(1)
        GCPAR =PROPS(2)
        THCK = PROPS(3)
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = FOUR
        DO I=1,INNODE
          AINTW(I) = ONE
        END DO
C
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initializing solution dependent variables (phase,history)
          DO I=1,NSTVTO
            SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUN(AN,dNdxi,XI)
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          DO INODE=1,NNODE
            BP(1,INODE)=dNdx(INODE,1)
            BP(2,INODE)=dNdx(INODE,2)
          END DO
C
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          PHASE=ZERO
          DPHASE=ZERO
          DO I=1,NDOFEL
            PHASE=PHASE+AN(I)*U(I)
          END DO
          DO I=1,NDOFEL
            DPHASE=DPHASE+AN(I)*DU(I,1)
          END DO
C
          IF (STEPITER.EQ.ZERO) THEN
            SDV(1)=PHASE-DPHASE
          ELSE
            SDV(1)=PHASE
          ENDIF
C         ! phase-field gradient
          DO I=1,2
            DP(I)=ZERO
          END DO
          DO I=1,2
            DO J=1,NNODE
              DP(I)=DP(I)+BP(I,J)*U(J)
            END DO
          END DO
C
C     ==================================================================
C     Calculating elastic ENERGY history
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            ENGN=USRVAR(JELEM,13,INPT)
          ELSE
            ENGN=USRVAR(JELEM,16,INPT)
          ENDIF
C
          HISTN=USRVAR(JELEM,16,INPT)
          IF (ENGN.GT.HISTN) THEN
            HIST=ENGN
          ELSE
            HIST=HISTN
          ENDIF
          SDV(2)=HIST
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
          DO I=1,NNODE
            DO K=1,NNODE
              DO J=1,2
                AMATRX(I,K)=AMATRX(I,K)+BP(J,I)*BP(J,K)*DTM*
     1            THCK*GCPAR*CLPAR*AINTW(INPT)
              END DO
              AMATRX(I,K)=AMATRX(I,K)+AN(I)*AN(K)*DTM*THCK*
     1          AINTW(INPT)*(GCPAR/CLPAR+TWO*HIST)
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO I=1,NDOFEL
            DO J=1,2
              RHS(I,1)=RHS(I,1)-BP(J,I)*DP(J)*GCPAR*CLPAR*
     1          AINTW(INPT)*DTM*THCK
            END DO
            RHS(I,1)=RHS(I,1)-AN(I)*AINTW(INPT)*DTM*THCK*
     1        ((GCPAR/CLPAR+TWO*HIST)*PHASE-TWO*HIST)
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTO
            SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM,I+NSTVTT,INPT)=SVARS(NSTVTO*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points
C
C     ==================================================================
C     ******************************************************************
C     Constructing Phase-field element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2P3 phase-field triangle element
C     ------------------------------------------------------------------
      ELSEIF (JTYPE.EQ.TWO) THEN ! C2P3 phase-field triangle element
C     ==================================================================
C     Time an iteration variables
C     ==================================================================
        TIMEZ=USRVAR(JELEM,17,1)
        IF (TIMEZ.LT.TIME(2)) THEN
          USRVAR(JELEM,17,1)=TIME(2)
          USRVAR(JELEM,18,1)=ZERO
        ELSE
          USRVAR(JELEM,18,1)=USRVAR(JELEM,18,1)+ONE
        ENDIF
        STEPITER=USRVAR(JELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        CLPAR=PROPS(1)
        GCPAR =PROPS(2)
        THCK = PROPS(3)
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = ONE
        AINTW(1) = HALF
C
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initializing solution dependent variables (phase,history)
          DO I=1,NSTVTO
            SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUNT(AN,dNdxi,XI)
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          DO INODE=1,NNODE
            BP(1,INODE)=dNdx(INODE,1)
            BP(2,INODE)=dNdx(INODE,2)
          END DO
C
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          PHASE=ZERO
          DPHASE=ZERO
          DO I=1,NDOFEL
            PHASE=PHASE+AN(I)*U(I)
          END DO
          DO I=1,NDOFEL
            DPHASE=DPHASE+AN(I)*DU(I,1)
          END DO
C
          IF (STEPITER.EQ.ZERO) THEN
            SDV(1)=PHASE-DPHASE
          ELSE
            SDV(1)=PHASE
          ENDIF
C         ! phase-field gradient
          DO I=1,2
            DP(I)=ZERO
          END DO
          DO I=1,2
            DO J=1,NNODE
              DP(I)=DP(I)+BP(I,J)*U(J)
            END DO
          END DO
C
C     ==================================================================
C     Calculating elastic ENERGY history
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            ENGN=USRVAR(JELEM,13,INPT)
          ELSE
            ENGN=USRVAR(JELEM,16,INPT)
          ENDIF
C
          HISTN=USRVAR(JELEM,16,INPT)
          IF (ENGN.GT.HISTN) THEN
            HIST=ENGN
          ELSE
            HIST=HISTN
          ENDIF
          SDV(2)=HIST
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
          DO I=1,NNODE
            DO K=1,NNODE
              DO J=1,2
                AMATRX(I,K)=AMATRX(I,K)+BP(J,I)*BP(J,K)*DTM*
     1            THCK*GCPAR*CLPAR*AINTW(INPT)
              END DO
              AMATRX(I,K)=AMATRX(I,K)+AN(I)*AN(K)*DTM*THCK*
     1          AINTW(INPT)*(GCPAR/CLPAR+TWO*HIST)
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO I=1,NDOFEL
            DO J=1,2
              RHS(I,1)=RHS(I,1)-BP(J,I)*DP(J)*GCPAR*CLPAR*
     1          AINTW(INPT)*DTM*THCK
            END DO
            RHS(I,1)=RHS(I,1)-AN(I)*AINTW(INPT)*DTM*THCK*
     1        ((GCPAR/CLPAR+TWO*HIST)*PHASE-TWO*HIST)
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTO
            SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM,I+NSTVTT,INPT)=SVARS(NSTVTO*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points







































C     ==================================================================
C     ******************************************************************
C     Constructing Plane-strain element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2D4 Isotropic Isothermal Elasticity with
C         isotropic energy degradation
C     ------------------------------------------------------------------
      ELSEIF (JTYPE .EQ. 21.D0) THEN
        STEPITER=USRVAR(JELEM-N_ELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        EMOD = PROPS(1)
        ENU = PROPS(2)
        THCK = PROPS(3)
        PARK = PROPS(4)
C
        EBULK3=EMOD/(ONE-TWO*ENU)
        EG2=EMOD/(ONE+ENU)
        EG=EG2/TWO
        EG3=THREE*EG
        ELAM=(EBULK3-EG2)/THREE
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = FOUR
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initial variables
          DO I=1,NSTVTT
            SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUN(AN,dNdxi,XI)
C     Shape functions
          IY=ZERO
          DO I = 1,NNODE
            IX=IY+1
            IY=IX+1
            VNI(1,IX)=AN(I)
            VNI(1,IY)=ZERO
            VNI(2,IX)=ZERO
            VNI(2,IY)=AN(I)
          END DO
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          IY=0
          DO INODE=1,NNODE
            IX=IY+1
            IY=IX+1
            BB(1,IX)= dNdx(INODE,1)
            BB(1,IY)= ZERO
            BB(2,IX)= ZERO
            BB(2,IY)= dNdx(INODE,2)
            BB(3,IX)= dNdx(INODE,2)
            BB(3,IY)= dNdx(INODE,1)
          END DO
C
C     ==================================================================
C     Calculating materials elastic stiffness matrix (plane strain)
C     ==================================================================
          DO I=1,3
            DO J=1,3
              CMAT(I,J)=ZERO
            END DO
          END DO
          CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
          CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
          CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
          CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
C     ==================================================================
C     Nodal displacements
C     ==================================================================
          DO J=1,2
            ULOC(J)=ZERO
          END DO
          DO J=1,2
            DO I=1,NDOFEL
              ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
            END DO
          END DO
          DO J=1,2
            SDV(J)=ULOC(J)
          END DO
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            PHASE=USRVAR(JELEM-N_ELEM,15,INPT)
          ELSE
            PHASE=USRVAR(JELEM-N_ELEM,14,INPT)
          ENDIF
C
          SDV(14)=PHASE
C     ==================================================================
C     Calculating strain
C     ==================================================================
          DO J=1,3
            EPS(J)=ZERO
          END DO
          DO I=1,3
            DO J=1,NDOFEL
              EPS(I)=EPS(I)+BB(I,J)*U(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+2)=EPS(J)
          END DO
C
C     ==================================================================
C     Calculating stresses
C     ==================================================================
          DO K1=1,3
            STRESS(K1)=ZERO
          END DO
          DO K1=1,3
            DO K2=1,3
              STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
            END DO
          END DO
          DO J=1,3
            SDV(J+5)=STRESS(J)*((ONE-PHASE)**TWO+PARK)
          END DO
          DO J=1,3
            SDV(J+8)=STRESS(J)
          END DO
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
          ENG=ZERO
          DO I=1,3
            ENG=ENG+STRESS(I)*EPS(I)*HALF
          END DO
          SDV(12)=ENG*((ONE-PHASE)**TWO+PARK)
          SDV(13)=ENG
          ENERGY(2)=ENG
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
          DO K=1,NDOFEL
            DO L=1,NDOFEL
              DO I=1,3
                DO J=1,3
                  AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1              BB(J,L)*DTM*THCK*((ONE-PHASE)**TWO+PARK)
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO K1=1,NDOFEL
            DO K4=1,3
              RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*
     1          STRESS(K4)*DTM*THCK*((ONE-PHASE)**TWO+PARK)
            END DO
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTT
            SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM-N_ELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points
C
C     ==================================================================
C     ******************************************************************
C     Constructing Plane-strain element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2D4 Isotropic Isothermal Elasticity with
C         anisotropic energy degradation by
C         spectral decomposition
C     ------------------------------------------------------------------
      ELSEIF (JTYPE .EQ. 22.D0) THEN
        STEPITER=USRVAR(JELEM-N_ELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        EMOD = PROPS(1)
        ENU = PROPS(2)
        THCK = PROPS(3)
        PARK = PROPS(4)
C
        ELAMEL=EMOD*ENU/((ONE+ENU)*(ONE-TWO*ENU))
        ELAMEG=EMOD/(TWO*(ONE+ENU))
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = FOUR
        DO I=1,INNODE
          AINTW(I) = ONE
        END DO
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initial variables
          DO I=1,NSTVTT
            SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUN(AN,dNdxi,XI)
C     Shape functions
          IY=ZERO
          DO I = 1,NNODE
            IX=IY+1
            IY=IX+1
            VNI(1,IX)=AN(I)
            VNI(1,IY)=ZERO
            VNI(2,IX)=ZERO
            VNI(2,IY)=AN(I)
          END DO
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          IY=0
          DO INODE=1,NNODE
            IX=IY+1
            IY=IX+1
            BB(1,IX)= dNdx(INODE,1)
            BB(1,IY)= ZERO
            BB(2,IX)= ZERO
            BB(2,IY)= dNdx(INODE,2)
            BB(3,IX)= dNdx(INODE,2)
            BB(3,IY)= dNdx(INODE,1)
          END DO
C     ==================================================================
C     Nodal displacements
C     ==================================================================
          DO J=1,2
            ULOC(J)=ZERO
          END DO
          DO J=1,2
            DO I=1,NDOFEL
              ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
            END DO
          END DO
          DO J=1,2
            SDV(J)=ULOC(J)
          END DO
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            PHASE=USRVAR(JELEM-N_ELEM,15,INPT)
          ELSE
            PHASE=USRVAR(JELEM-N_ELEM,14,INPT)
          ENDIF
C
          SDV(14)=PHASE
C     ==================================================================
C     Calculating strain
C     ==================================================================
          DO J=1,3
            EPS(J)=ZERO
          END DO
          DO I=1,3
            DO J=1,NDOFEL
              EPS(I)=EPS(I)+BB(I,J)*U(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+2)=EPS(J)
          END DO
C     ==================================================================
C     Calculating eigenvalue decomposition
C     ==================================================================
C
C    Only updating the stiffness matrix in the fist iteration step
C
          IF (STEPITER.LE.ONE) THEN
            DO K1=1,3
              EPSC(K1)=EPS(K1)
            END DO
          ELSE
            DO K1=1,3
              EPSC(K1)=USRVAR(JELEM-N_ELEM,2+K1,INPT)
            END DO
          ENDIF
C
          EPSMAX=MAXVAL(ABS(EPSC))
          IF (EPSMAX.GT.1e-8) THEN
            DE=EPSMAX*1e-4
          ELSE
            EPSC(1)=1e-8
            DE=1e-12
          ENDIF
C
          CALL EIG(EPSC,EIGV)
          SDV(9)=EIGV(1)
          SDV(10)=EIGV(2)
          SDV(11)=EIGV(1)+EIGV(2)
          ALPHA=ZERO
          IF ((EIGV(1)+EIGV(2)).GT.ZERO) THEN
            ALPHA=ONE
          ENDIF
          DO K1=1,2
            ALPHAI(K1)=ZERO
            IF (EIGV(K1).GT.ZERO) THEN
              ALPHAI(K1)=ONE
            ENDIF
          END DO
C
C     Derivatives of eigenvelus respect to original strain values
          CALL DIFFEIG(EPSC,DE,DEIGDEPS)
C
C     Second derivatives of eigenvelus respect to original strain values
          CALL DIFF2EIG(EPSC,DE,DEIG2DEPS)
C
C     ==================================================================
C     Calculating materials stiffness matrix
C     ==================================================================
          DO I=1,2
            DO J=1,2
              CLMAT(I,J)=ZERO
            END DO
          END DO
C
        GAMMA=ENU/(ONE-TWO*ENU)
        CLMAT(1,1)=((ONE-ALPHAI(1)*PHASE)**TWO+PARK)+GAMMA*
     1    ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,2)=((ONE-ALPHAI(2)*PHASE)**TWO+PARK)+GAMMA*
     1    ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,2)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,1)=CLMAT(1,2)
        DO I=1,2
          DO J=1,2
            CLMAT(I,J)=CLMAT(I,J)*EMOD/(ONE+ENU)
          END DO
        END DO
C
C     ==================================================================
C     Calculating degrated stiffness matrix
C     ==================================================================
          DO K1=1,3
            DO K2=1,3
              CMAT(K1,K2)=ZERO
            END DO
          END DO
C
          DO I=1,3
            DO L=1,3
              DO K=1,2
                DO J=1,2
            CMAT(I,L)=CMAT(I,L)+DEIGDEPS(J,I)*CLMAT(J,K)*DEIGDEPS(K,L)
                END DO
              END DO
            END DO
          END DO
          DO I=1,3
            DO L=1,3
              DO K=1,2
                DO J=1,2
            CMAT(I,L)=CMAT(I,L)+EIGV(K)*CLMAT(K,J)*DEIG2DEPS(I,L,J)
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Calculating stresses
C     ==================================================================
          DO K1=1,3
            STRESS(K1)=ZERO
          END DO
          DO I=1,3
            DO J=1,3
              STRESS(I)=STRESS(I)+CMAT(I,J)*EPS(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+5)=STRESS(J)
          END DO
C
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
          CALL EIG(EPS,EIGV)
          SDV(9)=EIGV(1)
          SDV(10)=EIGV(2)
          SDV(11)=EIGV(1)+EIGV(2)
          ENG=ZERO
          DO K2=1,3
            ENG=ENG+STRESS(K2)*EPS(K2)*HALF
          END DO
          SDV(12)=ENG
          SDV(13)=((ELAMEL*(ALPHA*(EIGV(1)+EIGV(2)))**TWO)/TWO+ELAMEG*
     1      ((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*ALPHAI(2))**TWO))
          ENERGY(2)=ENG
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
          DO K=1,NDOFEL
            DO L=1,NDOFEL
              DO I=1,3
                DO J=1,3
                  AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*
     1              CMAT(I,J)*BB(J,L)*DTM*THCK
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO K1=1,NDOFEL
            DO K4=1,3
              RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM*
     1          THCK
            END DO
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTT
            SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM-N_ELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points
C
C     ==================================================================
C     ******************************************************************
C     Constructing Plane-strain element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2D4 Isotropic Isothermal Elasticity with
C         anisotropic energy degradation by
C         volumetric-deviatoric decomposition
C     ------------------------------------------------------------------
      ELSEIF (JTYPE .EQ. 23.D0) THEN
        STEPITER=USRVAR(JELEM-N_ELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        EMOD = PROPS(1)
        ENU = PROPS(2)
        THCK = PROPS(3)
        PARK = PROPS(4)
C
        EBULK3=EMOD/(ONE-TWO*ENU)
        EG2=EMOD/(ONE+ENU)
        EG=EG2/TWO
        EG3=THREE*EG
        ELAM=(EBULK3-EG2)/THREE
        ebulk = EMOD/(ONE-TWO*ENU)/three
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = FOUR
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initial variables
          DO I=1,NSTVTT
            SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUN(AN,dNdxi,XI)
C     Shape functions
          IY=ZERO
          DO I = 1,NNODE
            IX=IY+1
            IY=IX+1
            VNI(1,IX)=AN(I)
            VNI(1,IY)=ZERO
            VNI(2,IX)=ZERO
            VNI(2,IY)=AN(I)
          END DO
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          IY=0
          DO INODE=1,NNODE
            IX=IY+1
            IY=IX+1
            BB(1,IX)= dNdx(INODE,1)
            BB(1,IY)= ZERO
            BB(2,IX)= ZERO
            BB(2,IY)= dNdx(INODE,2)
            BB(3,IX)= dNdx(INODE,2)
            BB(3,IY)= dNdx(INODE,1)
          END DO
C
C     ==================================================================
C     Nodal displacements
C     ==================================================================
          DO J=1,2
            ULOC(J)=ZERO
          END DO
          DO J=1,2
            DO I=1,NDOFEL
              ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
            END DO
          END DO
          DO J=1,2
            SDV(J)=ULOC(J)
          END DO
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            PHASE=USRVAR(JELEM-N_ELEM,15,INPT)
          ELSE
            PHASE=USRVAR(JELEM-N_ELEM,14,INPT)
          ENDIF
C
          SDV(14)=PHASE
C     ==================================================================
C     Calculating strain
C     ==================================================================
          DO J=1,3
            EPS(J)=ZERO
          END DO
          DO I=1,3
            DO J=1,NDOFEL
              EPS(I)=EPS(I)+BB(I,J)*U(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+2)=EPS(J)
          END DO
C     ==================================================================
C     Calculating volumetric-deviatoric decomposition
C     ==================================================================
      ! volumetric strain
          eps0 = zero
          do k1 = 1, 2
            eps0 = eps0 + eps(k1)
          enddo
          eps0 = eps0 / three
          ALPHA=ZERO
          IF (eps0.GT.ZERO) THEN
            ALPHA=ONE
          ENDIF
          ! deviatoric strain
          DO K1=1, 3
            EPSC(K1) = eps(k1)
          END DO
          do k1 = 1, 2
            EPSC(K1) = EPSC(K1) - eps0
          enddo
          sdv(9) = epsc(1)
          sdv(10) = epsc(2)
          sdv(11) = epsc(3)
C     ==================================================================
C     Calculating degrated materials elastic stiffness matrix
C     ==================================================================
          DO I=1,3
            DO J=1,3
              CMAT(I,J)=ZERO
            END DO
          END DO
C
          DO I=1,2
            DO J=1,2
              CMAT(I,J) = ebulk*alpha - eg*two/three
            END DO
            cmat(i,i) = ebulk*alpha + eg*four/three
          END DO
          cmat(3,3) = EG
          do k1 = 1, 3
            do k2 = 1, 3
              cmat(k1, K2) = cmat(k1, K2) * ((ONE-PHASE)**TWO+PARK)
            enddo
          enddo
          do k1 = 1, 2
            do k2 = 1, 2
              cmat(k1, K2) = cmat(k1, K2) + ebulk*(one-alpha)
            enddo
          enddo
C
C     ==================================================================
C     Calculating degraded stresses
C     ==================================================================
          DO K1=1,3
            STRESS(K1)=ZERO
          END DO
          DO K1=1,3
            DO K2=1,3
              STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
            END DO
          END DO
          DO J=1,3
            SDV(J+5)=STRESS(J)
          END DO
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
          ENG=ZERO
          DO I=1,3
            ENG=ENG+STRESS(I)*EPS(I)*HALF
          END DO
          SDV(12)=ENG
          ENERGY(2)=ENG
          ! phi_0^+
          eng0 = zero
          eng0 = half*ebulk*(alpha*eps0*three)**two
          do k1 = 1, 2
            eng0 = eng0 + eg*epsc(k1)*epsc(k1)
          enddo
          eng0 = eng0 + eg*epsc(3)*epsc(3)/two
          SDV(13)=eng0
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
          DO K=1,NDOFEL
            DO L=1,NDOFEL
              DO I=1,3
                DO J=1,3
                  AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1              BB(J,L)*DTM*THCK
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO K1=1,NDOFEL
            DO K4=1,3
              RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*
     1          STRESS(K4)*DTM*THCK
            END DO
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTT
            SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM-N_ELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points





































C     ==================================================================
C     ******************************************************************
C     Constructing Plane-strain element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2D3 Isotropic Isothermal Elasticity with
C         isotropic energy degradation
C     ------------------------------------------------------------------
      ELSEIF (JTYPE .EQ. 51.D0) THEN
        STEPITER=USRVAR(JELEM-N_ELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        EMOD = PROPS(1)
        ENU = PROPS(2)
        THCK = PROPS(3)
        PARK = PROPS(4)
C
        EBULK3=EMOD/(ONE-TWO*ENU)
        EG2=EMOD/(ONE+ENU)
        EG=EG2/TWO
        EG3=THREE*EG
        ELAM=(EBULK3-EG2)/THREE
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = ONE
        AINTW(1) = HALF
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initial variables
          DO I=1,NSTVTT
            SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUNT(AN,dNdxi,XI)
C     Shape functions
          IY=ZERO
          DO I = 1,NNODE
            IX=IY+1
            IY=IX+1
            VNI(1,IX)=AN(I)
            VNI(1,IY)=ZERO
            VNI(2,IX)=ZERO
            VNI(2,IY)=AN(I)
          END DO
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          IY=0
          DO INODE=1,NNODE
            IX=IY+1
            IY=IX+1
            BB(1,IX)= dNdx(INODE,1)
            BB(1,IY)= ZERO
            BB(2,IX)= ZERO
            BB(2,IY)= dNdx(INODE,2)
            BB(3,IX)= dNdx(INODE,2)
            BB(3,IY)= dNdx(INODE,1)
          END DO
C
C     ==================================================================
C     Calculating materials elastic stiffness matrix (plane strain)
C     ==================================================================
          DO I=1,3
            DO J=1,3
              CMAT(I,J)=ZERO
            END DO
          END DO
          CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
          CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
          CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
          CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
C     ==================================================================
C     Nodal displacements
C     ==================================================================
          DO J=1,2
            ULOC(J)=ZERO
          END DO
          DO J=1,2
            DO I=1,NDOFEL
              ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
            END DO
          END DO
          DO J=1,2
            SDV(J)=ULOC(J)
          END DO
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            PHASE=USRVAR(JELEM-N_ELEM,15,INPT)
          ELSE
            PHASE=USRVAR(JELEM-N_ELEM,14,INPT)
          ENDIF
C
          SDV(14)=PHASE
C     ==================================================================
C     Calculating strain
C     ==================================================================
          DO J=1,3
            EPS(J)=ZERO
          END DO
          DO I=1,3
            DO J=1,NDOFEL
              EPS(I)=EPS(I)+BB(I,J)*U(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+2)=EPS(J)
          END DO
C
C     ==================================================================
C     Calculating stresses
C     ==================================================================
          DO K1=1,3
            STRESS(K1)=ZERO
          END DO
          DO K1=1,3
            DO K2=1,3
              STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
            END DO
          END DO
          DO J=1,3
            SDV(J+5)=STRESS(J)*((ONE-PHASE)**TWO+PARK)
          END DO
          DO J=1,3
            SDV(J+8)=STRESS(J)
          END DO
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
          ENG=ZERO
          DO I=1,3
            ENG=ENG+STRESS(I)*EPS(I)*HALF
          END DO
          SDV(12)=ENG*((ONE-PHASE)**TWO+PARK)
          SDV(13)=ENG
          ENERGY(2)=ENG
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
          DO K=1,NDOFEL
            DO L=1,NDOFEL
              DO I=1,3
                DO J=1,3
                  AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1              BB(J,L)*DTM*THCK*((ONE-PHASE)**TWO+PARK)
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO K1=1,NDOFEL
            DO K4=1,3
              RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*
     1          STRESS(K4)*DTM*THCK*((ONE-PHASE)**TWO+PARK)
            END DO
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTT
            SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM-N_ELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points
C
C     ==================================================================
C     ******************************************************************
C     Constructing Plane-strain element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2D3 Isotropic Isothermal Elasticity with
C         anisotropic energy degradation by
C         spectral decomposition
C     ------------------------------------------------------------------
      ELSEIF (JTYPE .EQ. 52.D0) THEN
        STEPITER=USRVAR(JELEM-N_ELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        EMOD = PROPS(1)
        ENU = PROPS(2)
        THCK = PROPS(3)
        PARK = PROPS(4)
C
        ELAMEL=EMOD*ENU/((ONE+ENU)*(ONE-TWO*ENU))
        ELAMEG=EMOD/(TWO*(ONE+ENU))
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = ONE
        AINTW(1) = HALF
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initial variables
          DO I=1,NSTVTT
            SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUNT(AN,dNdxi,XI)
C     Shape functions
          IY=ZERO
          DO I = 1,NNODE
            IX=IY+1
            IY=IX+1
            VNI(1,IX)=AN(I)
            VNI(1,IY)=ZERO
            VNI(2,IX)=ZERO
            VNI(2,IY)=AN(I)
          END DO
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          IY=0
          DO INODE=1,NNODE
            IX=IY+1
            IY=IX+1
            BB(1,IX)= dNdx(INODE,1)
            BB(1,IY)= ZERO
            BB(2,IX)= ZERO
            BB(2,IY)= dNdx(INODE,2)
            BB(3,IX)= dNdx(INODE,2)
            BB(3,IY)= dNdx(INODE,1)
          END DO
C     ==================================================================
C     Nodal displacements
C     ==================================================================
          DO J=1,2
            ULOC(J)=ZERO
          END DO
          DO J=1,2
            DO I=1,NDOFEL
              ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
            END DO
          END DO
          DO J=1,2
            SDV(J)=ULOC(J)
          END DO
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            PHASE=USRVAR(JELEM-N_ELEM,15,INPT)
          ELSE
            PHASE=USRVAR(JELEM-N_ELEM,14,INPT)
          ENDIF
C
          SDV(14)=PHASE
C     ==================================================================
C     Calculating strain
C     ==================================================================
          DO J=1,3
            EPS(J)=ZERO
          END DO
          DO I=1,3
            DO J=1,NDOFEL
              EPS(I)=EPS(I)+BB(I,J)*U(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+2)=EPS(J)
          END DO
C     ==================================================================
C     Calculating eigenvalue decomposition
C     ==================================================================
C
C    Only updating the stiffness matrix in the fist iteration step
C
          IF (STEPITER.LE.ONE) THEN
            DO K1=1,3
              EPSC(K1)=EPS(K1)
            END DO
          ELSE
            DO K1=1,3
              EPSC(K1)=USRVAR(JELEM-N_ELEM,2+K1,INPT)
            END DO
          ENDIF
C
          EPSMAX=MAXVAL(ABS(EPSC))
          IF (EPSMAX.GT.1e-8) THEN
            DE=EPSMAX*1e-4
          ELSE
            EPSC(1)=1e-8
            DE=1e-12
          ENDIF
C
          CALL EIG(EPSC,EIGV)
          SDV(9)=EIGV(1)
          SDV(10)=EIGV(2)
          SDV(11)=EIGV(1)+EIGV(2)
          ALPHA=ZERO
          IF ((EIGV(1)+EIGV(2)).GT.ZERO) THEN
            ALPHA=ONE
          ENDIF
          DO K1=1,2
            ALPHAI(K1)=ZERO
            IF (EIGV(K1).GT.ZERO) THEN
              ALPHAI(K1)=ONE
            ENDIF
          END DO
C
C     Derivatives of eigenvelus respect to original strain values
          CALL DIFFEIG(EPSC,DE,DEIGDEPS)
C
C     Second derivatives of eigenvelus respect to original strain values
          CALL DIFF2EIG(EPSC,DE,DEIG2DEPS)
C
C     ==================================================================
C     Calculating materials stiffness matrix
C     ==================================================================
          DO I=1,2
            DO J=1,2
              CLMAT(I,J)=ZERO
            END DO
          END DO
C
        GAMMA=ENU/(ONE-TWO*ENU)
        CLMAT(1,1)=((ONE-ALPHAI(1)*PHASE)**TWO+PARK)+GAMMA*
     1    ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,2)=((ONE-ALPHAI(2)*PHASE)**TWO+PARK)+GAMMA*
     1    ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,2)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,1)=CLMAT(1,2)
        DO I=1,2
          DO J=1,2
            CLMAT(I,J)=CLMAT(I,J)*EMOD/(ONE+ENU)
          END DO
        END DO
C
C     ==================================================================
C     Calculating degrated stiffness matrix
C     ==================================================================
          DO K1=1,3
            DO K2=1,3
              CMAT(K1,K2)=ZERO
            END DO
          END DO
C
          DO I=1,3
            DO L=1,3
              DO K=1,2
                DO J=1,2
            CMAT(I,L)=CMAT(I,L)+DEIGDEPS(J,I)*CLMAT(J,K)*DEIGDEPS(K,L)
                END DO
              END DO
            END DO
          END DO
          DO I=1,3
            DO L=1,3
              DO K=1,2
                DO J=1,2
            CMAT(I,L)=CMAT(I,L)+EIGV(K)*CLMAT(K,J)*DEIG2DEPS(I,L,J)
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Calculating stresses
C     ==================================================================
          DO K1=1,3
            STRESS(K1)=ZERO
          END DO
          DO I=1,3
            DO J=1,3
              STRESS(I)=STRESS(I)+CMAT(I,J)*EPS(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+5)=STRESS(J)
          END DO
C
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
          CALL EIG(EPS,EIGV)
          SDV(9)=EIGV(1)
          SDV(10)=EIGV(2)
          SDV(11)=EIGV(1)+EIGV(2)
          ENG=ZERO
          DO K2=1,3
            ENG=ENG+STRESS(K2)*EPS(K2)*HALF
          END DO
          SDV(12)=ENG
          SDV(13)=((ELAMEL*(ALPHA*(EIGV(1)+EIGV(2)))**TWO)/TWO+ELAMEG*
     1      ((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*ALPHAI(2))**TWO))
          ENERGY(2)=ENG
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
          DO K=1,NDOFEL
            DO L=1,NDOFEL
              DO I=1,3
                DO J=1,3
                  AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*
     1              CMAT(I,J)*BB(J,L)*DTM*THCK
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO K1=1,NDOFEL
            DO K4=1,3
              RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM*
     1          THCK
            END DO
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTT
            SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM-N_ELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points
C
C     ==================================================================
C     ******************************************************************
C     Constructing Plane-strain element
C     ******************************************************************
C     ==================================================================
C     ------------------------------------------------------------------
C     C2D3 Isotropic Isothermal Elasticity with
C         anisotropic energy degradation by
C         volumetric-deviatoric decomposition
C     ------------------------------------------------------------------
      ELSEIF (JTYPE .EQ. 53.D0) THEN
        STEPITER=USRVAR(JELEM-N_ELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
        EMOD = PROPS(1)
        ENU = PROPS(2)
        THCK = PROPS(3)
        PARK = PROPS(4)
C
        EBULK3=EMOD/(ONE-TWO*ENU)
        EG2=EMOD/(ONE+ENU)
        EG=EG2/TWO
        EG3=THREE*EG
        ELAM=(EBULK3-EG2)/THREE
        ebulk = EMOD/(ONE-TWO*ENU)/three
C     ==================================================================
C     Initial preparations
C     ==================================================================
        DO K1 = 1, NDOFEL
          DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
          END DO
          DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
          END DO
        END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = ONE
        AINTW(1) = HALF
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        DO INPT=1,INNODE
C ======================================================================
C     Initial variables
          DO I=1,NSTVTT
            SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C
C     Local coordinates of the integration point
          XI(1) = XII(INPT,1)
          XI(2) = XII(INPT,2)
C     Shape functions and local derivatives
          CALL SHAPEFUNT(AN,dNdxi,XI)
C     Shape functions
          IY=ZERO
          DO I = 1,NNODE
            IX=IY+1
            IY=IX+1
            VNI(1,IX)=AN(I)
            VNI(1,IY)=ZERO
            VNI(2,IX)=ZERO
            VNI(2,IY)=AN(I)
          END DO
C     Jacobian
          DO I = 1,2
            DO J = 1,2
              VJACOB(I,J) = ZERO
              DO K = 1,NNODE
                VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
              END DO
            END DO
          END DO
C         ! calculate determinant of Jacobian matrix
          DTM = ZERO
          DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
C
          IF (DTM.LT.ZERO) THEN
            WRITE(7,*) 'Negative Jacobian',DTM
            CALL XIT
          ENDIF
C
C     Inverse of Jacobian
          VJABOBINV(1,1)=VJACOB(2,2)/DTM
          VJABOBINV(1,2)=-VJACOB(1,2)/DTM
          VJABOBINV(2,1)=-VJACOB(2,1)/DTM
          VJABOBINV(2,2)=VJACOB(1,1)/DTM
C
C     Derivatives of shape functions respect to global ccordinates
          DO K = 1,NNODE
            DO I = 1,2
              dNdx(K,I) = ZERO
              DO J = 1,2
                dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
              END DO
            END DO
          END DO
C
C     Calculating B matrix (B=LN)
          IY=0
          DO INODE=1,NNODE
            IX=IY+1
            IY=IX+1
            BB(1,IX)= dNdx(INODE,1)
            BB(1,IY)= ZERO
            BB(2,IX)= ZERO
            BB(2,IY)= dNdx(INODE,2)
            BB(3,IX)= dNdx(INODE,2)
            BB(3,IY)= dNdx(INODE,1)
          END DO
C
C     ==================================================================
C     Nodal displacements
C     ==================================================================
          DO J=1,2
            ULOC(J)=ZERO
          END DO
          DO J=1,2
            DO I=1,NDOFEL
              ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
            END DO
          END DO
          DO J=1,2
            SDV(J)=ULOC(J)
          END DO
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
          IF (STEPITER.EQ.ZERO) THEN
            PHASE=USRVAR(JELEM-N_ELEM,15,INPT)
          ELSE
            PHASE=USRVAR(JELEM-N_ELEM,14,INPT)
          ENDIF
C
          SDV(14)=PHASE
C     ==================================================================
C     Calculating strain
C     ==================================================================
          DO J=1,3
            EPS(J)=ZERO
          END DO
          DO I=1,3
            DO J=1,NDOFEL
              EPS(I)=EPS(I)+BB(I,J)*U(J)
            END DO
          END DO
          DO J=1,3
            SDV(J+2)=EPS(J)
          END DO
C     ==================================================================
C     Calculating volumetric-deviatoric decomposition
C     ==================================================================
      ! volumetric strain
          eps0 = zero
          do k1 = 1, 2
            eps0 = eps0 + eps(k1)
          enddo
          eps0 = eps0 / three
          ALPHA=ZERO
          IF (eps0.GT.ZERO) THEN
            ALPHA=ONE
          ENDIF
          ! deviatoric strain
          DO K1=1, 3
            EPSC(K1) = eps(k1)
          END DO
          do k1 = 1, 2
            EPSC(K1) = EPSC(K1) - eps0
          enddo
          sdv(9) = epsc(1)
          sdv(10) = epsc(2)
          sdv(11) = epsc(3)
C     ==================================================================
C     Calculating degrated materials elastic stiffness matrix
C     ==================================================================
          DO I=1,3
            DO J=1,3
              CMAT(I,J)=ZERO
            END DO
          END DO
C
          DO I=1,2
            DO J=1,2
              CMAT(I,J) = ebulk*alpha - eg*two/three
            END DO
            cmat(i,i) = ebulk*alpha + eg*four/three
          END DO
          cmat(3,3) = EG
          do k1 = 1, 3
            do k2 = 1, 3
              cmat(k1, K2) = cmat(k1, K2) * ((ONE-PHASE)**TWO+PARK)
            enddo
          enddo
          do k1 = 1, 2
            do k2 = 1, 2
              cmat(k1, K2) = cmat(k1, K2) + ebulk*(one-alpha)
            enddo
          enddo
C
C     ==================================================================
C     Calculating degraded stresses
C     ==================================================================
          DO K1=1,3
            STRESS(K1)=ZERO
          END DO
          DO K1=1,3
            DO K2=1,3
              STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
            END DO
          END DO
          DO J=1,3
            SDV(J+5)=STRESS(J)
          END DO
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
          ENG=ZERO
          DO I=1,3
            ENG=ENG+STRESS(I)*EPS(I)*HALF
          END DO
          SDV(12)=ENG
          ENERGY(2)=ENG
          ! phi_0^+
          eng0 = zero
          eng0 = half*ebulk*(alpha*eps0*three)**two
          do k1 = 1, 2
            eng0 = eng0 + eg*epsc(k1)*epsc(k1)
          enddo
          eng0 = eng0 + eg*epsc(3)*epsc(3)/two
          SDV(13)=eng0
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
          DO K=1,NDOFEL
            DO L=1,NDOFEL
              DO I=1,3
                DO J=1,3
                  AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1              BB(J,L)*DTM*THCK
                END DO
              END DO
            END DO
          END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
          DO K1=1,NDOFEL
            DO K4=1,3
              RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*
     1          STRESS(K4)*DTM*THCK
            END DO
          END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
          DO I=1,NSTVTT
            SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
            USRVAR(JELEM-N_ELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
          END DO
C ======================================================================
        END DO ! end of loop on integration points






      ENDIF

      RETURN
      END

C
C ======================================================================
C Shape functions for triangular elements
C ======================================================================
C
      SUBROUTINE SHAPEFUNT(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(3,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

C     Values of shape functions as a function of local coord.
      AN(1) = XI(1)
      AN(2) = XI(2)
      AN(3) = ONE-XI(1)-XI(2)
C
C     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
          dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  ONE
      dNdxi(1,2) =  ZERO
      dNdxi(2,1) =  ZERO
      dNdxi(2,2) =  ONE
      dNdxi(3,1) =  MONE
      dNdxi(3,2) =  MONE
C
      RETURN
      END
C
C ======================================================================
C Shape functions for brick elements C3D8
C ======================================================================
C
      SUBROUTINE SHAPEFUN(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      REAL*8 AN(4),dNdxi(4,2)
      REAL*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

C     Values of shape functions as a function of local coord.
      AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))

C     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END

C
C ======================================================================
C Eigenstrains from Voigt notation in 2D
C ======================================================================
C
      SUBROUTINE EIG(AMAT,EIGV)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.0)
      INTEGER I, J, K, NDIM
      REAL*8 AMAT(3),EIGV(2)
C
      EIGV(1)=ZERO
      EIGV(2)=ZERO
      EIGV(1)=(AMAT(1)+AMAT(2)-SQRT((AMAT(1)+AMAT(2))**TWO-FOUR*
     1  (AMAT(1)*AMAT(2)-AMAT(3)**TWO/FOUR)))/TWO;
      EIGV(2)=(AMAT(1)+AMAT(2)+SQRT((AMAT(1)+AMAT(2))**TWO-FOUR*
     1  (AMAT(1)*AMAT(2)-AMAT(3)**TWO/FOUR)))/TWO;
C
      RETURN
      END
C
C ======================================================================
C Derivatives of the eigenvalues respect to the original matrix
C ======================================================================
C
      SUBROUTINE DIFFEIG(EPS,DE,DEIGDEPS)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.0)
      INTEGER I, J, K, NDIM
      REAL*8 EPS(3),DEIGDEPS(2,3),EIGINI(2),EPSNEW(3),EIGNEW(2)
C
      CALL EIG(EPS,EIGINI)
C
      DO J=1,3
        DO I=1,3
          EPSNEW(I)=EPS(I)
        END DO
        EPSNEW(J)=EPSNEW(J)+DE
        CALL EIG(EPSNEW,EIGNEW)
        DO K=1,2
          DEIGDEPS(K,J)=(EIGNEW(K)-EIGINI(K))/DE
        END DO
      END DO
C
      RETURN
      END
C
C ======================================================================
C Second derivatives of the eigenvalues respect to the original matrix
C ======================================================================
C
      SUBROUTINE DIFF2EIG(EPS,DE,DEIG2DEPS)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.0)
      INTEGER I, J, K, L, NDIM
      REAL*8 EPS(3),DEIG2DEPS(3,3,2),EIGINI(2),EPSUP(3),EPSDW(3),
     1 EIGUP(2),EIGDW(2),EPS21(3),EPS12(3),EIG21(2),EIG12(2)
C
      CALL EIG(EPS,EIGINI)
C
C Diagonal elements
      DO J=1,3
        DO I=1,3
          EPSUP(I)=EPS(I)
          EPSDW(I)=EPS(I)
        END DO
        EPSUP(J)=EPSUP(J)+DE
        EPSDW(J)=EPSDW(J)-DE
        CALL EIG(EPSUP,EIGUP)
        CALL EIG(EPSDW,EIGDW)
        DO K=1,2
          DEIG2DEPS(J,J,K)=(EIGUP(K)-TWO*EIGINI(K)+EIGDW(K))/DE**TWO
        END DO
      END DO
C
C Of-diagonal elements
      DO J=1,3
      DO L=J+1,3
        DO I=1,3
          EPSUP(I)=EPS(I)
          EPSDW(I)=EPS(I)
          EPS12(I)=EPS(I)
          EPS21(I)=EPS(I)
        END DO
        EPSUP(J)=EPSUP(J)+DE
        EPSUP(L)=EPSUP(L)+DE
        EPSDW(J)=EPSDW(J)-DE
        EPSDW(L)=EPSDW(L)-DE
        EPS12(J)=EPS12(J)+DE
        EPS12(L)=EPS12(L)-DE
        EPS21(J)=EPS21(J)-DE
        EPS21(L)=EPS21(L)+DE
        CALL EIG(EPSUP,EIGUP)
        CALL EIG(EPSDW,EIGDW)
        CALL EIG(EPS12,EIG12)
        CALL EIG(EPS21,EIG21)
        DO K=1,2
          DEIG2DEPS(J,L,K)=(EIGUP(K)-EIG21(K)-EIG12(K)+EIGDW(K))/
     1      DE**TWO/FOUR
          DEIG2DEPS(L,J,K)=DEIG2DEPS(J,L,K)
        END DO
      END DO
      END DO
C
      RETURN
      END
C
C Subroutine UMAT  :
C Dummy material
C
C ==============================================================
C !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
C ==============================================================
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF=0.5,
     1 N_ELEM=4434,NSTV=18)
      DATA NEWTON,TOLER/40,1.D-6/
C
      COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)
C
C -----------------------------------------------------------
C          Material properties
C -----------------------------------------------------------
C          PROPS(1) - Young's modulus
C          PROPS(2) - Poisson ratio
C -----------------------------------------------------------
C
C   Elastic properties
C
      EMOD=PROPS(1)
      ENU=PROPS(2)
      EG=EMOD/(TWO*(ONE+ENU))
      EG2=EG*TWO
      ELAM=EG2*ENU/(ONE-TWO*ENU)
C
C   Stiffness tensor
C
      DO K1=1, NTENS
        DO K2=1, NTENS
          DDSDDE(K2, K1)=0.0
        END DO
      END DO
C
      DO K1=1, NDI
        DO K2=1, NDI
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
C
      DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
      END DO
C
C     Calculate Stresses
C
      DO K1=1, NTENS
        DO K2=1, NTENS
          STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
      END DO
C
      NELEMAN=NOEL-TWO*N_ELEM
      IF (NPT.EQ.3) THEN
        NPT=4
      ELSEIF (NPT.EQ.4) THEN
        NPT=3
      ENDIF
C
      DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
      END DO
C
      RETURN
      END
