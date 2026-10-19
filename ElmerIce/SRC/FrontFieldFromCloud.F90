!------------------------------------------------------------------------------
! FrontFieldFromCloud
!
! Reads scattered 3D data clouds (columns:  X  Y  Z  ... value ...) and
! interpolates the value onto the nodes of a boundary of a 3D mesh (e.g. the
! Bowdoin calving front). Spatial interpolation is k-nearest-neighbour inverse-
! distance weighting (Shepard) in true 3D -- no external nn/csa libraries.
!
! Time handling: the data files are named <Prefix><itime><Suffix>, one per
! supported integer data-time. Three modes, selected by the SIF keywords:
!
!   (a) TIME INTERPOLATION  -- set "Data Time End" (with "Data Time Start").
!       At first call the solver probes which files actually exist in
!       [Start,End] (stride "Data Time Step", default 1) and records their
!       times. The list may be NON-UNIFORM: missing files (gaps) are simply
!       absent from the list. For an arbitrary simulation time t the two
!       bracketing AVAILABLE times t_lo <= t <= t_hi are located and blended:
!           value = (1-w)*field(t_lo) + w*field(t_hi),  w=(t-t_lo)/(t_hi-t_lo).
!       Across a gap the blend simply spans the wider interval (e.g. 180..184
!       if 182 is missing). Both bracketing clouds are held resident (2 slots)
!       and reused across steps; a new file is read only when a bracket
!       endpoint changes. The spatial kNN is re-evaluated every step against
!       the CURRENT node positions, so this is correct on a moving front mesh.
!       "Data Time Digits" zero-pads the index (e.g. 4 -> Pressure_t_0178.csv).
!
!   (b) SNAP (legacy)       -- "Data File Prefix" but no "Data Time Step":
!       one file per step, itime = NINT(t * multiplier).
!
!   (c) STATIC              -- "Data File" (single file, no time dependence).
!
! Fast kNN: each cloud is indexed with a uniform cubic bucket grid (counting-
! sort CSR); a ring-expanding query stops once the k-th neighbour is provably
! nearest -- O(1) average instead of the O(ndata) brute-force scan.
!
! Hydrostatic correction (optional): the total pressure in the cloud carries a
! strong hydrostatic gradient P_hyd(z) = rho*g*(SeaLevel - z) for z < SeaLevel
! (0 above), which, combined with a varying input point density, seeds
! interpolation artefacts. With "Subtract Hydrostatic Pressure = Logical True"
! the analytic P_hyd(z) (using the cloud's z, column 3) is removed from each
! datum BEFORE the interpolant is built, so only the smooth DYNAMIC pressure is
! interpolated. Set "Restore Hydrostatic Pressure = Logical True" to add the
! analytic P_hyd(z) back at each mesh node afterwards (artefact-free TOTAL);
! leave it off to output the dynamic part alone. Keywords: "Water Density"
! (default 1025.0), "Gravity" (default 9.81), "Sea Level" (default 0.0).
!
! File format: whitespace- OR comma-separated; blank lines and lines that do
! not start with a number (e.g. a "X, Y, Z, Press_tot" header, or #/! comments)
! are skipped automatically.
!
! Compile:   elmerf90 FrontFieldFromCloud.F90 -o FrontFieldFromCloud
!
! SIF usage (time interpolation):
!
!   Solver 5
!     Exec Solver = "Before Timestep"
!     Equation = "FrontField"
!     Procedure = "FrontFieldFromCloud" "FrontFieldFromCloud"
!     Variable = -nooutput "dummy_ffc"
!     Variable DOFs = 1
!     Exported Variable 1 = "FrontPressure"
!     Exported Variable 1 DOFs = 1
!
!     Data File Prefix = File "loads/forces_Bowdoin_v1/Pressure_t_"
!     Data File Suffix = File ".csv"
!     Value Column    = Integer 4       ! Press_tot (coords are cols 1..3)
!     Data Time Start = Integer 2       ! first candidate time (clamp below)
!     Data Time End   = Integer 470     ! last candidate time  (clamp above) -> enables blending
!     Data Time Step  = Integer 2       ! probe stride over candidates (default 1)
!   ! Data Time Digits = Integer 0      ! zero-pad width; 0 = no padding (default)
!   ! Data Time Multiplier = Real 1.0   ! sim-time -> data-time scaling (default 1)
!     Number Of Neighbours = Integer 8
!     IDW Power = Real 2.0
!     Snap Tolerance = Real 1.0e-3
!     Subtract Hydrostatic Pressure = Logical True   ! interpolate dynamic part only
!   ! Restore Hydrostatic Pressure  = Logical True   ! add P_hyd(z) back -> total
!   ! Water Density = Real 1025.0       ! rho [kg/m3]  (default 1025)
!   ! Gravity       = Real 9.81         ! g   [m/s2]   (default 9.81)
!   ! Sea Level     = Real 0.0          ! z of sea level (default 0)
!   ! Points Per Cell = Real 1.0
!   ! Grid Cell Size  = Real 20.0
!   End
!
!   Boundary Condition 1
!     Target Boundaries(1) = 1
!     Front Interpolation = Logical True
!   End
!------------------------------------------------------------------------------

!==============================================================================
MODULE FrontCloudUtils
!==============================================================================
  USE DefUtils
  IMPLICIT NONE

  ! one scattered cloud + its bucket grid + cache key
  TYPE CloudData_t
    CHARACTER(LEN=MAX_NAME_LEN) :: fname = ''
    INTEGER :: ndata = 0, knn = 6
    INTEGER :: gnx = 1, gny = 1, gnz = 1, ncells = 1
    REAL(KIND=dp) :: gh = 1.0_dp, gx0 = 0.0_dp, gy0 = 0.0_dp, gz0 = 0.0_dp
    REAL(KIND=dp), ALLOCATABLE :: cx(:), cy(:), cz(:), cv(:)
    INTEGER,       ALLOCATABLE :: cellstart(:), pidx(:)
  END TYPE CloudData_t

CONTAINS

  !--- read "X Y Z ... value ..." (space/comma sep); skip header/comment lines
  SUBROUTINE ReadCloud( cd, fname, valcol )
    TYPE(CloudData_t), INTENT(INOUT) :: cd
    CHARACTER(LEN=*),  INTENT(IN)    :: fname
    INTEGER,           INTENT(IN)    :: valcol
    CHARACTER(LEN=1024) :: line
    REAL(KIND=dp) :: row(MAX(3,valcol))
    INTEGER :: u, ios, k, ncol
    ncol = MAX(3, valcol)
    IF ( ALLOCATED(cd%cx) ) DEALLOCATE( cd%cx, cd%cy, cd%cz, cd%cv )

    OPEN(NEWUNIT=u, FILE=TRIM(fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF ( ios /= 0 ) CALL Fatal('FrontFieldFromCloud', &
         'cannot open data file "'//TRIM(fname)//'"')

    cd%ndata = 0                                ! pass 1: count data lines
    DO
      READ(u,'(A)',IOSTAT=ios) line
      IF ( ios /= 0 ) EXIT
      IF ( IsDataLine(line) ) cd%ndata = cd%ndata + 1
    END DO
    IF ( cd%ndata < 1 ) CALL Fatal('FrontFieldFromCloud', &
         'no numeric data lines in "'//TRIM(fname)//'"')

    REWIND(u)
    ALLOCATE( cd%cx(cd%ndata), cd%cy(cd%ndata), cd%cz(cd%ndata), cd%cv(cd%ndata) )
    k = 0                                       ! pass 2: parse
    DO
      READ(u,'(A)',IOSTAT=ios) line
      IF ( ios /= 0 ) EXIT
      IF ( .NOT. IsDataLine(line) ) CYCLE
      k = k + 1
      READ(line,*,IOSTAT=ios) row(1:ncol)
      IF ( ios /= 0 ) THEN
        WRITE(Message,'(A,I0,A)') 'could not read ',ncol,' columns from a data line'
        CALL Fatal('FrontFieldFromCloud', Message)
      END IF
      cd%cx(k)=row(1) ; cd%cy(k)=row(2) ; cd%cz(k)=row(3) ; cd%cv(k)=row(valcol)
    END DO
    CLOSE(u)
    WRITE(Message,'(A,A)') 'Read Datafile: ', TRIM(fname)
    CALL Info('FrontFieldFromCloud', Message, Level=5)
  END SUBROUTINE ReadCloud

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

  !--- linear cell index (0-based) of cloud point p ------------------------
  FUNCTION CellOf( cd, p ) RESULT( c )
    TYPE(CloudData_t), INTENT(IN) :: cd
    INTEGER,           INTENT(IN) :: p
    INTEGER :: c, ix, iy, iz
    ix = MIN( cd%gnx-1, MAX( 0, INT( (cd%cx(p)-cd%gx0)/cd%gh ) ) )
    iy = MIN( cd%gny-1, MAX( 0, INT( (cd%cy(p)-cd%gy0)/cd%gh ) ) )
    iz = MIN( cd%gnz-1, MAX( 0, INT( (cd%cz(p)-cd%gz0)/cd%gh ) ) )
    c  = ix + cd%gnx*( iy + cd%gny*iz )
  END FUNCTION CellOf

  !--- build uniform cubic bucket grid over the cloud (counting-sort CSR) --
  SUBROUTINE BuildGrid( cd, ppc, cellsize )
    TYPE(CloudData_t), INTENT(INOUT) :: cd
    REAL(KIND=dp),     INTENT(IN)    :: ppc, cellsize
    INTEGER :: p, c
    INTEGER, ALLOCATABLE :: cnt(:), cur(:)
    REAL(KIND=dp) :: xmin,xmax,ymin,ymax,zmin,zmax, dx,dy,dz, vol, spanmax

    xmin = MINVAL(cd%cx) ; xmax = MAXVAL(cd%cx)
    ymin = MINVAL(cd%cy) ; ymax = MAXVAL(cd%cy)
    zmin = MINVAL(cd%cz) ; zmax = MAXVAL(cd%cz)
    dx = xmax-xmin ; dy = ymax-ymin ; dz = zmax-zmin
    spanmax = MAX( dx, dy, dz, TINY(1.0_dp) )

    IF ( cellsize > 0.0_dp ) THEN
      cd%gh = cellsize
    ELSE
      vol = MAX(dx,1.0e-8_dp*spanmax) * MAX(dy,1.0e-8_dp*spanmax) &
            * MAX(dz,1.0e-8_dp*spanmax)
      cd%gh = ( vol * ppc / REAL(cd%ndata,dp) ) ** (1.0_dp/3.0_dp)
    END IF
    cd%gh = MAX( cd%gh, 1.0e-8_dp*spanmax )

    cd%gx0 = xmin ; cd%gy0 = ymin ; cd%gz0 = zmin
    cd%gnx = MAX(1,INT(dx/cd%gh)+1)
    cd%gny = MAX(1,INT(dy/cd%gh)+1)
    cd%gnz = MAX(1,INT(dz/cd%gh)+1)
    DO WHILE ( REAL(cd%gnx,dp)*REAL(cd%gny,dp)*REAL(cd%gnz,dp) &
               > 8.0_dp*REAL(cd%ndata,dp)+1000.0_dp )
      cd%gh = cd%gh * 1.259921_dp               ! *2^(1/3)
      cd%gnx = MAX(1,INT(dx/cd%gh)+1)
      cd%gny = MAX(1,INT(dy/cd%gh)+1)
      cd%gnz = MAX(1,INT(dz/cd%gh)+1)
    END DO
    cd%ncells = cd%gnx*cd%gny*cd%gnz

    IF ( ALLOCATED(cd%cellstart) ) DEALLOCATE(cd%cellstart)
    IF ( ALLOCATED(cd%pidx) )      DEALLOCATE(cd%pidx)
    ALLOCATE( cd%cellstart(0:cd%ncells), cd%pidx(cd%ndata) )
    ALLOCATE( cnt(0:cd%ncells-1) ) ; cnt = 0

    DO p = 1, cd%ndata                          ! histogram
      c = CellOf(cd,p) ; cnt(c) = cnt(c) + 1
    END DO
    cd%cellstart(0) = 0                         ! prefix sum -> CSR offsets
    DO c = 1, cd%ncells
      cd%cellstart(c) = cd%cellstart(c-1) + cnt(c-1)
    END DO
    ALLOCATE( cur(0:cd%ncells-1) ) ; cur = cd%cellstart(0:cd%ncells-1)
    DO p = 1, cd%ndata                          ! scatter point ids
      c = CellOf(cd,p) ; cur(c) = cur(c) + 1 ; cd%pidx(cur(c)) = p
    END DO
    DEALLOCATE( cnt, cur )
  END SUBROUTINE BuildGrid

  !--- exact kNN IDW at (xp,yp,zp) via ring-expanding bucket search ---------
  SUBROUTINE Interp( cd, xp, yp, zp, power, snaptol, res )
    TYPE(CloudData_t), INTENT(IN)  :: cd
    REAL(KIND=dp),     INTENT(IN)  :: xp, yp, zp, power, snaptol
    REAL(KIND=dp),     INTENT(OUT) :: res
    REAL(KIND=dp) :: bestd(cd%knn), d, dd, w, wsum, vsum
    INTEGER :: besti(cd%knn), ix, iy, iz, r, rmax, di, dj, dk, jx, jy, jz
    INTEGER :: c, pp, p, k, pos

    bestd = HUGE(1.0_dp) ; besti = 0
    ix = MIN( cd%gnx-1, MAX( 0, INT( (xp-cd%gx0)/cd%gh ) ) )
    iy = MIN( cd%gny-1, MAX( 0, INT( (yp-cd%gy0)/cd%gh ) ) )
    iz = MIN( cd%gnz-1, MAX( 0, INT( (zp-cd%gz0)/cd%gh ) ) )
    rmax = MAX( ix, cd%gnx-1-ix, iy, cd%gny-1-iy, iz, cd%gnz-1-iz )

    r = 0
    DO
      DO dk = -r, r                             ! scan the shell at radius r
        jz = iz + dk ; IF ( jz < 0 .OR. jz >= cd%gnz ) CYCLE
        DO dj = -r, r
          jy = iy + dj ; IF ( jy < 0 .OR. jy >= cd%gny ) CYCLE
          DO di = -r, r
            IF ( MAX(ABS(di),ABS(dj),ABS(dk)) /= r ) CYCLE   ! shell surface only
            jx = ix + di ; IF ( jx < 0 .OR. jx >= cd%gnx ) CYCLE
            c = jx + cd%gnx*( jy + cd%gny*jz )
            DO pp = cd%cellstart(c)+1, cd%cellstart(c+1)
              p = cd%pidx(pp)
              d = (xp-cd%cx(p))**2 + (yp-cd%cy(p))**2 + (zp-cd%cz(p))**2
              IF ( d < bestd(cd%knn) ) THEN
                pos = cd%knn
                DO WHILE ( pos > 1 )
                  IF ( d < bestd(pos-1) ) THEN
                    bestd(pos)=bestd(pos-1) ; besti(pos)=besti(pos-1) ; pos=pos-1
                  ELSE
                    EXIT
                  END IF
                END DO
                bestd(pos) = d ; besti(pos) = p
              END IF
            END DO
          END DO
        END DO
      END DO
      ! done once the k-th neighbour is closer than any unscanned shell can be
      IF ( besti(cd%knn) /= 0 ) THEN
        IF ( SQRT(bestd(cd%knn)) <= REAL(r,dp)*cd%gh ) EXIT
      END IF
      IF ( r >= rmax ) EXIT
      r = r + 1
    END DO

    wsum = 0.0_dp ; vsum = 0.0_dp
    DO k = 1, cd%knn
      IF ( besti(k) == 0 ) CYCLE
      dd = SQRT(bestd(k))
      IF ( dd < snaptol ) THEN                   ! node sits on a data point
        res = cd%cv(besti(k)) ; RETURN
      END IF
      w = 1.0_dp / dd**power
      wsum = wsum + w
      vsum = vsum + w*cd%cv(besti(k))
    END DO
    res = vsum / wsum
  END SUBROUTINE Interp

END MODULE FrontCloudUtils


!==============================================================================
SUBROUTINE FrontFieldFromCloud( Model, Solver, dt, TransientSimulation )
!==============================================================================
  USE DefUtils
  USE FrontCloudUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp)  :: dt
  LOGICAL        :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(Mesh_t),     POINTER :: Mesh
  TYPE(Variable_t), POINTER :: Var, TimeVar
  TYPE(Element_t),  POINTER :: Element
  TYPE(ValueList_t),POINTER :: BC, SolverParams
  LOGICAL, ALLOCATABLE      :: Done(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, Prefix, Suffix, NumStr
  CHARACTER(LEN=MAX_NAME_LEN) :: StaticFile, FloName, FhiName
  CHARACTER(LEN=MAX_NAME_LEN) :: Caller = 'FrontFieldFromCloud'
  LOGICAL :: Found, GotPrefix, HaveEnd, Temporal
  INTEGER :: t, i, node, n, nset, ia, ib
  INTEGER :: itime, tend, ilo, ihi, kk
  REAL(KIND=dp) :: xq, yq, zq, val, vlo, vhi, rtime, dtime, tmult, wblend

  ! --- two resident cloud slots + settings (persist across calls) ---------
  TYPE(CloudData_t), SAVE :: slot(2)
  INTEGER, ALLOCATABLE, SAVE :: avail(:)         ! sorted list of available data times
  INTEGER, SAVE :: knn_req = 6, valcol = 4, digits = 0, navail = 0
  REAL(KIND=dp), SAVE :: power = 2.0_dp, snaptol = 1.0e-9_dp
  REAL(KIND=dp), SAVE :: ppc = 1.0_dp, cellsize = 0.0_dp
  REAL(KIND=dp), SAVE :: rhow = 1025.0_dp, grav = 9.81_dp, sealev = 0.0_dp
  LOGICAL, SAVE :: HaveParams = .FALSE., HaveAvail = .FALSE.
  LOGICAL, SAVE :: SubHydro = .FALSE., AddHydro = .FALSE.
!------------------------------------------------------------------------------

  Mesh => Solver % Mesh
  SolverParams => GetSolverParams(Solver)

  CALL INFO(Caller,"******************************************",Level=3)
  CALL INFO(Caller,"Start",Level=3)
  CALL INFO(Caller,"******************************************",Level=3)

  ! ------------------------------------------------- read static params once
  IF ( .NOT. HaveParams ) THEN
    knn_req = GetInteger( SolverParams, 'Number Of Neighbours', Found )
    IF ( .NOT. Found ) knn_req = 6
    power = GetConstReal( SolverParams, 'IDW Power', Found )
    IF ( .NOT. Found ) power = 2.0_dp
    snaptol = GetConstReal( SolverParams, 'Snap Tolerance', Found )
    IF ( .NOT. Found ) snaptol = 1.0e-9_dp
    valcol = GetInteger( SolverParams, 'Value Column', Found )
    IF ( .NOT. Found ) valcol = 4
    ppc = GetConstReal( SolverParams, 'Points Per Cell', Found )
    IF ( .NOT. Found ) ppc = 1.0_dp
    cellsize = GetConstReal( SolverParams, 'Grid Cell Size', Found )
    IF ( .NOT. Found ) cellsize = 0.0_dp
    digits = GetInteger( SolverParams, 'Data Time Digits', Found )
    IF ( .NOT. Found ) digits = 0
    SubHydro = GetLogical( SolverParams, 'Subtract Hydrostatic Pressure', Found )
    IF ( .NOT. Found ) SubHydro = .FALSE.
    AddHydro = GetLogical( SolverParams, 'Restore Hydrostatic Pressure', Found )
    IF ( .NOT. Found ) AddHydro = .FALSE.
    rhow = GetConstReal( SolverParams, 'Water Density', Found )
    IF ( .NOT. Found ) rhow = 1025.0_dp
    grav = GetConstReal( SolverParams, 'Gravity', Found )
    IF ( .NOT. Found ) grav = 9.81_dp
    sealev = GetConstReal( SolverParams, 'Sea Level', Found )
    IF ( .NOT. Found ) sealev = 0.0_dp
    IF ( SubHydro ) THEN
      WRITE(Message,'(A,F0.2,A,F0.3,A,F0.2)') 'hydrostatic correction ON: rho=', &
           rhow,' g=',grav,' sea level=',sealev
      CALL Info(Caller, Message, Level=4)
    END IF
    HaveParams = .TRUE.
  END IF

  ! --------------------------------------------------- current (data) time
  Prefix = GetString( SolverParams, 'Data File Prefix', GotPrefix )
  Suffix = GetString( SolverParams, 'Data File Suffix', Found )
  IF ( .NOT. Found ) Suffix = ''
  tmult = GetConstReal( SolverParams, 'Data Time Multiplier', Found )
  IF ( .NOT. Found ) tmult = 1.0_dp
  rtime = GetConstReal( SolverParams, 'Data Time', Found )
  IF ( .NOT. Found ) THEN
    TimeVar => VariableGet( Mesh % Variables, 'Time' )
    IF ( ASSOCIATED(TimeVar) ) THEN
      rtime = TimeVar % Values(1)
    ELSE
      rtime = 0.0_dp
    END IF
  END IF
  dtime = rtime * tmult

  tend = GetInteger( SolverParams, 'Data Time End', HaveEnd )
  Temporal = GotPrefix .AND. HaveEnd            ! a Data Time End enables blending

  ! ------------------------------------------- resolve bracketing file names
  wblend = 0.0_dp
  IF ( .NOT. GotPrefix ) THEN                   ! (c) static single file
    StaticFile = GetString( SolverParams, 'Data File', Found )
    IF ( .NOT. Found ) CALL Fatal(Caller, &
         'give "Data File Prefix" (+time) or a static "Data File"')
    FloName = StaticFile ; FhiName = StaticFile

  ELSE IF ( Temporal ) THEN                     ! (a) time interpolation
    IF ( .NOT. HaveAvail ) CALL BuildAvail()    ! probe which files exist (once)
    ! bracket dtime within the (possibly NON-UNIFORM) list of available times
    IF ( dtime <= REAL(avail(1),dp) ) THEN
      ilo = avail(1) ; ihi = avail(1)
    ELSE IF ( dtime >= REAL(avail(navail),dp) ) THEN
      ilo = avail(navail) ; ihi = avail(navail)
    ELSE
      kk = 1
      DO i = 1, navail-1                         ! avail() is ascending
        IF ( REAL(avail(i),dp) <= dtime .AND. dtime < REAL(avail(i+1),dp) ) THEN
          kk = i ; EXIT
        END IF
      END DO
      ilo = avail(kk) ; ihi = avail(kk+1)
      wblend = (dtime - REAL(ilo,dp)) / REAL(ihi-ilo,dp)
    END IF
    IF ( wblend <= 0.0_dp ) ihi = ilo           ! exact hit / clamp -> one file
    FloName = MakeName(ilo) ; FhiName = MakeName(ihi)

  ELSE                                          ! (b) snap to nearest (legacy)
    itime = NINT( dtime )
    FloName = MakeName(itime)
    FhiName = FloName
  END IF

  ! --------------------------------- make both brackets resident in slot 1/2
  ia = SlotHolding( FloName )
  IF ( TRIM(FloName) == TRIM(FhiName) ) THEN
    IF ( ia < 0 ) THEN
      ia = 1 ; CALL LoadSlot( 1, FloName )
    END IF
    ib = ia
  ELSE
    ib = SlotHolding( FhiName )
    IF ( ia < 0 ) THEN                          ! load lo into the slot != ib
      IF ( ib == 1 ) THEN ; ia = 2 ; ELSE ; ia = 1 ; END IF
      CALL LoadSlot( ia, FloName )
    END IF
    IF ( ib < 0 ) THEN                          ! load hi into the other slot
      IF ( ia == 1 ) THEN ; ib = 2 ; ELSE ; ib = 1 ; END IF
      CALL LoadSlot( ib, FhiName )
    END IF
  END IF

  WRITE(Message,'(A,ES12.4,A,A,A,A,A,F8.5)') 'data time=',dtime, &
       '  lo="',TRIM(slot(ia)%fname),'" hi="',TRIM(slot(ib)%fname),'" w=',wblend
  CALL Info(Caller, Message, Level=4)

  ! ------------------------------------------------------------- output field
  VarName = GetString( SolverParams, 'Exported Variable 1', Found )
  IF ( .NOT. Found ) VarName = 'FrontField'
  Var => VariableGet( Mesh % Variables, VarName )
  IF ( .NOT. ASSOCIATED(Var) ) &
       CALL Fatal(Caller,'Exported variable "'//TRIM(VarName)//'" not found')

  ! ------------------------------- loop the front nodes, blend lo/hi in time
  ALLOCATE( Done(Mesh % NumberOfNodes) )
  Done = .FALSE.
  nset = 0

  DO t = 1, Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)
    IF ( .NOT. ActiveBoundaryElement() ) CYCLE
    BC => GetBC()
    IF ( .NOT. ASSOCIATED(BC) ) CYCLE
    IF ( .NOT. GetLogical(BC,'Front Interpolation',Found) ) CYCLE

    n = GetElementNOFNodes()
    DO i = 1, n
      node = Element % NodeIndexes(i)
      IF ( Done(node) ) CYCLE
      IF ( Var % Perm(node) <= 0 ) CYCLE
      xq = Mesh % Nodes % x(node)
      yq = Mesh % Nodes % y(node)
      zq = Mesh % Nodes % z(node)
      CALL Interp( slot(ia), xq, yq, zq, power, snaptol, vlo )
      IF ( ib /= ia ) THEN
        CALL Interp( slot(ib), xq, yq, zq, power, snaptol, vhi )
        val = (1.0_dp - wblend)*vlo + wblend*vhi
      ELSE
        val = vlo
      END IF
      IF ( AddHydro ) val = val + Phydro( zq )     ! rebuild total at node z
      Var % Values( Var % Perm(node) ) = val
      Done(node) = .TRUE.
      nset = nset + 1
    END DO
  END DO

  DEALLOCATE( Done )
  WRITE(Message,'(A,I0,A)') 'Set field on ',nset,' front nodes'
  CALL Info(Caller, Message, Level=4)

CONTAINS

  !--- hydrostatic seawater pressure at height z (0 above sea level) -------
  FUNCTION Phydro( zz ) RESULT( ph )
    REAL(KIND=dp), INTENT(IN) :: zz
    REAL(KIND=dp) :: ph
    IF ( zz < sealev ) THEN
      ph = rhow * grav * ( sealev - zz )
    ELSE
      ph = 0.0_dp
    END IF
  END FUNCTION Phydro

  !--- index of the slot currently holding "name", or -1 -------------------
  FUNCTION SlotHolding( name ) RESULT( idx )
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: idx, s
    idx = -1
    DO s = 1, 2
      IF ( slot(s)%ndata > 0 .AND. TRIM(slot(s)%fname) == TRIM(name) ) THEN
        idx = s ; RETURN
      END IF
    END DO
  END FUNCTION SlotHolding

  !--- (re)load a data file into slot s and (re)build its grid -------------
  SUBROUTINE LoadSlot( s, name )
    INTEGER,          INTENT(IN) :: s
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: j
    CALL ReadCloud( slot(s), name, valcol )
    IF ( SubHydro ) THEN                          ! strip hydrostatic -> dynamic
      DO j = 1, slot(s)%ndata
        slot(s)%cv(j) = slot(s)%cv(j) - Phydro( slot(s)%cz(j) )
      END DO
    END IF
    slot(s)%knn = MIN( knn_req, slot(s)%ndata )
    CALL BuildGrid( slot(s), ppc, cellsize )
    slot(s)%fname = name
    WRITE(Message,'(A,I0,A,A,A,I0,A)') 'slot ',s,' <- "',TRIM(name), &
         '"  (',slot(s)%ndata,' pts)'
    CALL Info(Caller, Message, Level=5)
  END SUBROUTINE LoadSlot

  !--- build "<Prefix><it><Suffix>", zero-padded to "digits" if digits>0 ---
  FUNCTION MakeName( it ) RESULT( name )
    INTEGER, INTENT(IN) :: it
    CHARACTER(LEN=MAX_NAME_LEN) :: name
    CHARACTER(LEN=32) :: numstr, fmt
    IF ( digits > 0 ) THEN
      WRITE( fmt, '(A,I0,A,I0,A)' ) '(I', digits, '.', digits, ')'
      WRITE( numstr, fmt ) it
    ELSE
      WRITE( numstr, '(I0)' ) it
    END IF
    name = TRIM(Prefix) // TRIM(ADJUSTL(numstr)) // TRIM(Suffix)
  END FUNCTION MakeName

  !--- probe [Start,End] once and record the times whose files exist -------
  SUBROUTINE BuildAvail()
    INTEGER :: it, j, tstart, tstep
    LOGICAL :: ex, f
    tstart = GetInteger( SolverParams, 'Data Time Start', f )
    IF ( .NOT. f ) CALL Fatal(Caller,'time interpolation needs "Data Time Start"')
    tstep = GetInteger( SolverParams, 'Data Time Step', f )
    IF ( .NOT. f .OR. tstep < 1 ) tstep = 1

    navail = 0                                   ! pass 1: count existing files
    DO it = tstart, tend, tstep
      INQUIRE( FILE=TRIM(MakeName(it)), EXIST=ex )
      IF ( ex ) navail = navail + 1
    END DO
    IF ( navail < 1 ) CALL Fatal(Caller, &
         'no data files found in the "Data Time" range (check path/digits/suffix)')

    IF ( ALLOCATED(avail) ) DEALLOCATE(avail)
    ALLOCATE( avail(navail) )
    j = 0                                        ! pass 2: record (ascending)
    DO it = tstart, tend, tstep
      INQUIRE( FILE=TRIM(MakeName(it)), EXIST=ex )
      IF ( ex ) THEN ; j = j + 1 ; avail(j) = it ; END IF
    END DO
    HaveAvail = .TRUE.
    WRITE(Message,'(A,I0,A,I0,A,I0,A,I0,A)') 'found ',navail, &
         ' data files in [',tstart,',',tend,'] (stride ',tstep,')'
    CALL Info(Caller, Message, Level=4)
  END SUBROUTINE BuildAvail

END SUBROUTINE FrontFieldFromCloud
!------------------------------------------------------------------------------
