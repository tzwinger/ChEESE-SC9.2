!------------------------------------------------------------------------------
! IcebergFrontLoad
!
! Turns the tabulated TOTAL contact force of an iceberg (Newtons) into surface
! tractions on the calving front, applied in a boundary condition as
!   Normal Force = Equals <normal traction variable>
!   Force 3      = Equals <vertical tangential traction variable>
! ("External Pressure" is left free -- it carries the point-cloud pressure from
!  FrontFieldFromCloud; the two normal contributions simply add on the BC.)
!
! Sign: the raw contact is compressive (a push). With the default OUTWARD
! normals of a single-body run, "Normal Force" loads along +normal (outward),
! so an inward push needs a NEGATIVE value. The normal output therefore carries
! "Normal Force Sign" (default -1.0) so that a push comes out negative and adds
! correctly to the (positive, compressive) External Pressure. Use +1.0 if you
! instead feed a compressive-positive keyword such as External Pressure.
!
! Inputs
!   * Force File   : columns  t[s]  Fn[N]  Ft[N]   (normal, vertical-tangential;
!                    text headers / comment lines are skipped automatically).
!   * Motion File  : same times in column 1; the vertical position z_att(t) of
!                    the iceberg attachment line in column "Motion Z Column".
!   * Constants    : Start Coordinate(3), Stop Coordinate(3) -- the constant
!                    (x,y) endpoints of the attachment line (z ignored).
!
! Method (per timestep, at the current simulation Time t):
!   1. Linearly interpolate Fn(t), Ft(t) from the force table and z_att(t) from
!      the motion table. If t lies outside EITHER table -> both tractions set to
!      zero and an INFO message (Level 3) is written; nothing else is done.
!   2. The contact patch is the flagged front boundary ("Iceberg Attachment")
!      whose nodes (a) project onto the Start->Stop segment with parameter in
!      [0,1] and (b) have |z - z_att| <= 0.5*"Attachment Area Bandwidth".
!   3. Unilateral contact: an iceberg can only PUSH the front, not pull it.
!      Compression is taken to be positive (Force Sign)*Fn. When that is
!      negative the iceberg has separated -> both tractions are set to zero
!      (no contact => no normal load and no friction). Disable with
!      "Unilateral Normal Contact = Logical False".
!   4. The patch area A is obtained by integrating the nodal in-band indicator
!      over the boundary faces (A = int phi dA). The tractions are then
!         normal_traction = Fn / A ,   tangential_traction = Ft / A
!      on the band nodes (zero elsewhere), the normal one also times
!      "Normal Force Sign". Because A = int phi dA is computed with the SAME
!      shape functions Elmer uses, the applied load integrates back to the
!      tabulated totals: |int NormalForce dA| = Fn, int Force3 dA = Ft.
!   5. Element-wise contact mask: the solver CREATES an elemental variable
!      (default "Iceberg Contact", override "Contact Mask Name") with one DOF
!      per BOUNDARY FACE -- bulk elements are excluded (Perm=0). +1 on faces
!      carrying a band node while the load is finite (in contact), -1 on the
!      other flagged faces / on the no-load returns. It is built in code, NOT
!      via "Exported Variable = -elem", because a bulk-active solver's -elem
!      export is permuted over bulk elements and would miss the faces. Do NOT
!      declare it in the SIF; it appears automatically in the output.
!
! The patch (and thus A) is recomputed every step, so it correctly follows the
! moving z_att(t) and a deforming front mesh.
!
! Compile:  elmerf90 IcebergFrontLoad.F90 -o IcebergFrontLoad
!
! SIF usage:
!
!   Constants
!     Start Coordinate(3) = Real 509577.5 8621975.0 0.0
!     Stop Coordinate(3)  = Real 509232.5 8622540.0 0.0
!   End
!
!   Solver 6
!     Exec Solver = "Before Timestep"
!     Equation = "IcebergLoad"
!     Procedure = "IcebergFrontLoad" "IcebergFrontLoad"
!     Variable = -nooutput "dummy_ibl"
!     Variable DOFs = 1
!     Exported Variable 1 = "IcebergNormalTraction"
!     Exported Variable 1 DOFs = 1
!     Exported Variable 2 = "IcebergForce3"
!     Exported Variable 2 DOFs = 1
!     ! NB: the "Iceberg Contact" element mask is created by the solver over the
!     !     boundary faces -- do NOT declare it as an Exported Variable here.
!
!     Force File  = File "Fc.dat"
!     Motion File = File "Mvt_Frame_body_TopLeftCorner.dat"
!     Motion Z Column = Integer 3
!     Attachment Area Bandwidth = Real 20.0   ! [m] vertical extent of the band
!   ! Force Sign = Real 1.0                   ! sign s.t. compression = positive s*Fn
!   ! Normal Force Sign = Real -1.0           ! push -> negative for "Normal Force"
!   ! Unilateral Normal Contact = Logical True! push-only; tensile phase -> zero
!   ! Vertical Offset = Real 0.0              ! added to z_att (datum shift)
!   ! Contact Mask Name = String "Iceberg Contact"  ! name of the -elem mask var
!   End
!
!   Boundary Condition 1                      ! the front
!     Target Boundaries(1) = 1
!     Flow Force BC = Logical True
!     Front Interpolation = Logical True      ! FrontFieldFromCloud (point cloud)
!     Iceberg Attachment  = Logical True      ! this solver
!     External Pressure = Equals FrontPressure          ! point-cloud pressure
!     Normal Force      = Equals IcebergNormalTraction  ! iceberg normal (signed)
!     Force 3           = Equals IcebergForce3          ! iceberg vertical
!   End
!------------------------------------------------------------------------------
SUBROUTINE IcebergFrontLoad( Model, Solver, dt, TransientSimulation )
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp)  :: dt
  LOGICAL        :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(Mesh_t),     POINTER :: Mesh
  TYPE(Variable_t), POINTER :: VarN, VarT, TimeVar, MaskVar
  TYPE(Element_t),  POINTER :: Element
  TYPE(ValueList_t),POINTER :: BC, SolverParams
  TYPE(Nodes_t), SAVE       :: Nodes
  TYPE(GaussIntegrationPoints_t) :: IP
  REAL(KIND=dp), POINTER    :: Carr(:,:)
  LOGICAL, ALLOCATABLE      :: InBand(:)
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), dummyB(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: ForceFile, MotionFile, nName, tName
  CHARACTER(LEN=MAX_NAME_LEN) :: Caller = 'IcebergFrontLoad'
  LOGICAL :: Found, stat, InRangeF, InRangeM, HaveMask, elemActive, anyForce
  INTEGER :: t, i, gp, n, nmax, node, nb, ei, nactive, e, nbe, ntot
  REAL(KIND=dp) :: tq, Fn, Ft, zatt, area, tracN, tracT
  REAL(KIND=dp) :: sdx, sdy, slen2, detJ, wdA, phigp, sproj, xn, yn, zn

  ! --- persistent tables + settings ---------------------------------------
  REAL(KIND=dp), ALLOCATABLE, SAVE :: fcT(:), fcN(:), fcF(:), mvT(:), mvZ(:)
  INTEGER, SAVE :: nfc = 0, nmv = 0, mzcol = 3
  REAL(KIND=dp), SAVE :: bandw = 0.0_dp, fsign = 1.0_dp, voff = 0.0_dp, ncsign = -1.0_dp
  REAL(KIND=dp), SAVE :: sx1, sy1, sx2, sy2
  LOGICAL, SAVE :: Initialized = .FALSE., Unilateral = .TRUE.
  ! element-wise contact mask, created over the boundary faces by this solver
  CHARACTER(LEN=MAX_NAME_LEN), SAVE :: mName = ''
  INTEGER,       POINTER, SAVE :: MaskPerm(:) => NULL()
  REAL(KIND=dp), POINTER, SAVE :: MaskVal(:)  => NULL()
!------------------------------------------------------------------------------

  Mesh => Solver % Mesh
  SolverParams => GetSolverParams(Solver)

  ! ------------------------------------------------- one-time setup / reads
  IF ( .NOT. Initialized ) THEN
    ForceFile = GetString( SolverParams, 'Force File', Found )
    IF ( .NOT. Found ) CALL Fatal(Caller,'"Force File" missing')
    MotionFile = GetString( SolverParams, 'Motion File', Found )
    IF ( .NOT. Found ) CALL Fatal(Caller,'"Motion File" missing')
    mzcol = GetInteger( SolverParams, 'Motion Z Column', Found )
    IF ( .NOT. Found ) mzcol = 3
    bandw = GetConstReal( SolverParams, 'Attachment Area Bandwidth', Found )
    IF ( .NOT. Found ) CALL Fatal(Caller,'"Attachment Area Bandwidth" missing')
    fsign = GetConstReal( SolverParams, 'Force Sign', Found )
    IF ( .NOT. Found ) fsign = 1.0_dp
    ncsign = GetConstReal( SolverParams, 'Normal Force Sign', Found )
    IF ( .NOT. Found ) ncsign = -1.0_dp
    voff = GetConstReal( SolverParams, 'Vertical Offset', Found )
    IF ( .NOT. Found ) voff = 0.0_dp
    Unilateral = GetLogical( SolverParams, 'Unilateral Normal Contact', Found )
    IF ( .NOT. Found ) Unilateral = .TRUE.

    CALL ReadTable( ForceFile, 2, 3, fcT, fcN, fcF )   ! t, Fn, Ft
    CALL ReadTable( MotionFile, mzcol, mzcol, mvT, mvZ, dummyB )
    nfc = SIZE(fcT) ; nmv = SIZE(mvT)

    Carr => ListGetConstRealArray( Model % Constants, 'Start Coordinate', Found )
    IF ( .NOT. Found ) CALL Fatal(Caller,'"Start Coordinate" missing in Constants')
    CALL GetXY( Carr, sx1, sy1 )
    Carr => ListGetConstRealArray( Model % Constants, 'Stop Coordinate', Found )
    IF ( .NOT. Found ) CALL Fatal(Caller,'"Stop Coordinate" missing in Constants')
    CALL GetXY( Carr, sx2, sy2 )

    ! --- create the element-wise contact mask over the BOUNDARY FACES --------
    ! One DOF per boundary element; bulk elements get Perm=0 (excluded). Built
    ! here rather than via "Exported Variable = -elem" because a bulk-active
    ! solver's -elem export would not cover the boundary faces.
    mName = GetString( SolverParams, 'Contact Mask Name', Found )
    IF ( .NOT. Found ) mName = 'Iceberg Contact'
    IF ( .NOT. ASSOCIATED( VariableGet(Mesh % Variables, mName) ) ) THEN
      nbe  = Mesh % NumberOfBoundaryElements
      ntot = Mesh % NumberOfBulkElements + nbe
      ALLOCATE( MaskPerm(ntot) ) ; MaskPerm = 0
      DO e = 1, nbe
        MaskPerm( Mesh % NumberOfBulkElements + e ) = e
      END DO
      ALLOCATE( MaskVal(nbe) ) ; MaskVal = 0.0_dp
      CALL VariableAdd( Mesh % Variables, Mesh, Solver, mName, 1, &
           MaskVal, MaskPerm, TYPE = Variable_on_elements )
      CALL Info(Caller,'created elemental contact mask "'//TRIM(mName)// &
           '" on the boundary faces', Level=4)
    END IF

    Initialized = .TRUE.
    WRITE(Message,'(A,I0,A,I0,A)') 'read ',nfc,' force rows, ',nmv,' motion rows'
    CALL Info(Caller, Message, Level=4)
  END IF

  ! -------------------------------------------------------- output variables
  nName = GetString( SolverParams, 'Exported Variable 1', Found )
  IF ( .NOT. Found ) nName = 'IcebergNormalTraction'
  tName = GetString( SolverParams, 'Exported Variable 2', Found )
  IF ( .NOT. Found ) tName = 'IcebergForce3'
  VarN => VariableGet( Mesh % Variables, nName )
  VarT => VariableGet( Mesh % Variables, tName )
  IF ( .NOT. ASSOCIATED(VarN) .OR. .NOT. ASSOCIATED(VarT) ) &
       CALL Fatal(Caller,'exported traction variables not found')
  VarN % Values = 0.0_dp
  VarT % Values = 0.0_dp

  ! element-wise contact mask:  +1 active (face sees a finite force), -1 passive.
  ! Reset to passive on the flagged FACES only, so early returns leave it passive
  ! and bulk is never touched.
  MaskVar => VariableGet( Mesh % Variables, mName )
  HaveMask = ASSOCIATED(MaskVar)
  IF ( HaveMask ) THEN
    DO t = 1, Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE
      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE
      IF ( .NOT. GetLogical(BC,'Iceberg Attachment',Found) ) CYCLE
      ei = Element % ElementIndex
      IF ( ei >= 1 .AND. ei <= SIZE(MaskVar % Perm) ) THEN
        IF ( MaskVar % Perm(ei) > 0 ) MaskVar % Values( MaskVar % Perm(ei) ) = -1.0_dp
      END IF
    END DO
  END IF

  ! ---------------------------------------------- interpolate force + motion
  TimeVar => VariableGet( Mesh % Variables, 'Time' )
  IF ( ASSOCIATED(TimeVar) ) THEN ; tq = TimeVar % Values(1) ; ELSE ; tq = 0.0_dp ; END IF

  CALL Interp1( fcT, fcN, nfc, tq, Fn,   InRangeF )
  CALL Interp1( fcT, fcF, nfc, tq, Ft,   InRangeF )
  CALL Interp1( mvT, mvZ, nmv, tq, zatt, InRangeM )

  IF ( .NOT. ( InRangeF .AND. InRangeM ) ) THEN
    WRITE(Message,'(A,ES12.4,A)') 'time=',tq, &
         ' s outside tabulated force/motion range - zero forces applied'
    CALL Info(Caller, Message, Level=3)
    RETURN                                       ! variables already zeroed
  END IF

  Fn = fsign*Fn ; Ft = fsign*Ft
  IF ( Unilateral .AND. Fn < 0.0_dp ) THEN       ! tensile -> iceberg not in contact
    WRITE(Message,'(A,ES12.4,A)') 'time=',tq, &
         ' s: tensile normal - iceberg separated, zero contact force'
    CALL Info(Caller, Message, Level=4)
    RETURN                                       ! traction variables already zeroed
  END IF
  zatt = zatt + voff

  ! ------------------------------- pass 1: mark band nodes + integrate area
  ALLOCATE( InBand(Mesh % NumberOfNodes) ) ; InBand = .FALSE.
  nmax = Mesh % MaxElementNodes
  ALLOCATE( Basis(nmax), dBasisdx(nmax,3) )
  sdx = sx2-sx1 ; sdy = sy2-sy1 ; slen2 = sdx*sdx + sdy*sdy
  area = 0.0_dp

  DO t = 1, Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)
    IF ( .NOT. ActiveBoundaryElement() ) CYCLE
    BC => GetBC()
    IF ( .NOT. ASSOCIATED(BC) ) CYCLE
    IF ( .NOT. GetLogical(BC,'Iceberg Attachment',Found) ) CYCLE
    n = GetElementNOFNodes()

    DO i = 1, n                                  ! band membership (nodal)
      node = Element % NodeIndexes(i)
      xn = Mesh % Nodes % x(node)
      yn = Mesh % Nodes % y(node)
      zn = Mesh % Nodes % z(node)
      IF ( slen2 > 0.0_dp ) THEN
        sproj = ( (xn-sx1)*sdx + (yn-sy1)*sdy ) / slen2
      ELSE
        sproj = 0.0_dp
      END IF
      InBand(node) = ( sproj >= 0.0_dp .AND. sproj <= 1.0_dp .AND. &
                       ABS(zn - zatt) <= 0.5_dp*bandw )
    END DO

    CALL GetElementNodes( Nodes, Element )       ! integrate phi over the face
    IP = GaussPoints( Element )
    DO gp = 1, IP % n
      stat = ElementInfo( Element, Nodes, IP%u(gp), IP%v(gp), IP%w(gp), &
                          detJ, Basis, dBasisdx )
      wdA = IP % s(gp) * detJ
      phigp = 0.0_dp
      DO i = 1, n
        IF ( InBand(Element % NodeIndexes(i)) ) phigp = phigp + Basis(i)
      END DO
      area = area + wdA*phigp
    END DO
  END DO

  IF ( area <= 0.0_dp ) THEN
    WRITE(Message,'(A,ES12.4,A,ES12.4,A)') 'no contact area at z_att=',zatt, &
         ' (bandwidth ',bandw,') - zero forces applied'
    CALL Info(Caller, Message, Level=3)
    DEALLOCATE( InBand, Basis, dBasisdx )
    RETURN
  END IF

  tracN = ncsign * Fn/area ; tracT = Ft/area
  anyForce = ( ABS(Fn) > 0.0_dp .OR. ABS(Ft) > 0.0_dp )

  ! ---------------- pass 2: write tractions on band nodes + set element mask
  nb = 0 ; nactive = 0
  DO t = 1, Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)
    IF ( .NOT. ActiveBoundaryElement() ) CYCLE
    BC => GetBC()
    IF ( .NOT. ASSOCIATED(BC) ) CYCLE
    IF ( .NOT. GetLogical(BC,'Iceberg Attachment',Found) ) CYCLE
    n = GetElementNOFNodes()
    elemActive = .FALSE.
    DO i = 1, n
      node = Element % NodeIndexes(i)
      IF ( .NOT. InBand(node) ) CYCLE
      IF ( VarN % Perm(node) > 0 ) VarN % Values(VarN % Perm(node)) = tracN
      IF ( VarT % Perm(node) > 0 ) VarT % Values(VarT % Perm(node)) = tracT
      nb = nb + 1
      elemActive = .TRUE.                        ! this face carries a band node
    END DO
    ! active if the face is in the loaded band AND the load is finite
    IF ( HaveMask .AND. elemActive .AND. anyForce ) THEN
      ei = Element % ElementIndex
      IF ( ei >= 1 .AND. ei <= SIZE(MaskVar % Perm) ) THEN
        IF ( MaskVar % Perm(ei) > 0 ) THEN
          MaskVar % Values( MaskVar % Perm(ei) ) = 1.0_dp
          nactive = nactive + 1
        END IF
      END IF
    END IF
  END DO

  WRITE(Message,'(A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4)') 'Fn=',Fn, &
       ' Ft=',Ft,' area=',area,' tN=',tracN,' tT=',tracT
  CALL Info(Caller, Message, Level=4)
  WRITE(Message,'(A,I0,A,I0,A,ES12.4)') 'iceberg load on ',nb, &
       ' band nodes / ',nactive,' active elements at z_att=',zatt
  CALL Info(Caller, Message, Level=4)

  DEALLOCATE( InBand, Basis, dBasisdx )

CONTAINS

  !--- x,y from a Constants coordinate array, robust to (3,1) or (1,3) -----
  SUBROUTINE GetXY( arr, x, y )
    REAL(KIND=dp), INTENT(IN)  :: arr(:,:)
    REAL(KIND=dp), INTENT(OUT) :: x, y
    x = arr(1,1)
    IF ( SIZE(arr,1) >= 2 ) THEN
      y = arr(2,1)
    ELSE
      y = arr(1,2)
    END IF
  END SUBROUTINE GetXY

  !--- read numeric table: col1 -> T, colA -> A, colB -> B -----------------
  SUBROUTINE ReadTable( fname, vcolA, vcolB, T, A, B )
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER,          INTENT(IN) :: vcolA, vcolB
    REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: T(:), A(:), B(:)
    CHARACTER(LEN=1024) :: line
    REAL(KIND=dp) :: row(MAX(vcolA,vcolB))
    INTEGER :: u, ios, k, nrows, ncol
    ncol = MAX(vcolA,vcolB)
    OPEN(NEWUNIT=u, FILE=TRIM(fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF ( ios /= 0 ) CALL Fatal('IcebergFrontLoad','cannot open "'//TRIM(fname)//'"')
    nrows = 0
    DO
      READ(u,'(A)',IOSTAT=ios) line ; IF ( ios /= 0 ) EXIT
      IF ( IsDataLine(line) ) nrows = nrows + 1
    END DO
    IF ( nrows < 1 ) CALL Fatal('IcebergFrontLoad','no data in "'//TRIM(fname)//'"')
    REWIND(u)
    ALLOCATE( T(nrows), A(nrows), B(nrows) )
    k = 0
    DO
      READ(u,'(A)',IOSTAT=ios) line ; IF ( ios /= 0 ) EXIT
      IF ( .NOT. IsDataLine(line) ) CYCLE
      k = k + 1
      READ(line,*,IOSTAT=ios) row(1:ncol)
      IF ( ios /= 0 ) CALL Fatal('IcebergFrontLoad', &
           'bad numeric row in "'//TRIM(fname)//'"')
      T(k) = row(1) ; A(k) = row(vcolA) ; B(k) = row(vcolB)
    END DO
    CLOSE(u)
  END SUBROUTINE ReadTable

  !--- true if the (trimmed) line begins a numeric record ------------------
  FUNCTION IsDataLine( line ) RESULT( ok )
    CHARACTER(LEN=*), INTENT(IN) :: line
    LOGICAL :: ok
    CHARACTER(LEN=1) :: c
    CHARACTER(LEN=:), ALLOCATABLE :: s
    s = ADJUSTL(line)
    ok = .FALSE.
    IF ( LEN_TRIM(s) == 0 ) RETURN
    c = s(1:1)
    ok = ( (c >= '0' .AND. c <= '9') .OR. c == '+' .OR. c == '-' .OR. c == '.' )
  END FUNCTION IsDataLine

  !--- linear interpolation of V(T) at tq; ascending T; out-of-range -> F --
  SUBROUTINE Interp1( T, V, n, tq, val, inrange )
    REAL(KIND=dp), INTENT(IN)  :: T(:), V(:), tq
    INTEGER,       INTENT(IN)  :: n
    REAL(KIND=dp), INTENT(OUT) :: val
    LOGICAL,       INTENT(OUT) :: inrange
    INTEGER :: j
    val = 0.0_dp ; inrange = .FALSE.
    IF ( n < 1 ) RETURN
    IF ( tq < T(1) .OR. tq > T(n) ) RETURN
    inrange = .TRUE.
    DO j = 1, n-1
      IF ( tq >= T(j) .AND. tq <= T(j+1) ) THEN
        IF ( T(j+1) > T(j) ) THEN
          val = V(j) + (V(j+1)-V(j)) * (tq-T(j)) / (T(j+1)-T(j))
        ELSE
          val = V(j)
        END IF
        RETURN
      END IF
    END DO
    val = V(n)                                   ! tq == T(n)
  END SUBROUTINE Interp1

END SUBROUTINE IcebergFrontLoad
!------------------------------------------------------------------------------
