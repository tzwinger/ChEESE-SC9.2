SUBROUTINE ReadIsisData( Model,Solver,dt,TransientSimulation )
  USE DefUtils  
  USE ElementDescription
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt

  TYPE(Variable_t), POINTER :: VarFtang, VarFnorm, VarZcoord, VarTime
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="ReadIsisData"
  REAL(KIND=dp), POINTER :: Ftang(:), Fnorm(:), Zcoord(:), Time(:)
  INTEGER, POINTER :: FtangPerm(:), FnormPerm(:), ZcoordPerm(:)
  INTEGER, SAVE :: datafieldentries, I, switch=0
  REAL(KIND=dp) :: CurrentTime
  REAL(KIND=dp), SAVE, POINTER :: datafield(:,:) 
  LOGICAL, SAVE :: FirstTime=.TRUE.

  
  IF (FirstTime) THEN
    FirstTime = .FALSE.
    CALL read_fc(datafield,datafieldentries)
    WRITE (Message,*) "Read files with", datafieldentries, "entries"
    CALL INFO(SolverName, "****************************************************", Level=1)
    CALL INFO(SolverName, Message, Level=1)
    CALL INFO(SolverName, "****************************************************", Level=1)
    DO I=1,datafieldentries
      PRINT *, I,":",datafield(1:4,I)
    END DO
    CALL INFO(SolverName, "****************************************************", Level=1)
  END IF

  VarZcoord => Solver % Variable
  IF (.NOT.ASSOCIATED(VarZcoord)) THEN
    CALL FATAL(SolverName,'Solver Variable (vertical cordinate) not found')
  END IF
  Zcoord => VarZcoord % Values
  ZcoordPerm => VarZcoord % Perm
  
  VarFtang => VariableGet( Solver % Mesh % Variables,'Vertical Tangential Force')
  IF (.NOT.ASSOCIATED(VarFtang)) THEN
    CALL FATAL(SolverName,'Variable "Vertical Tangential Force" not found')
  END IF
  FtangPerm => VarFtang % Perm
  Ftang => VarFtang % Values
  
  VarFnorm => VariableGet( Solver % Mesh % Variables,'Normal Force')
  IF (.NOT.ASSOCIATED(VarFnorm)) THEN
    CALL FATAL(SolverName,'Variable "Normal Force" not found')
  END IF
  FnormPerm => VarFnorm % Perm
  Fnorm => VarFnorm % Values

  VarTime => VariableGet( Solver % Mesh % Variables,"Time")
  IF (.NOT.ASSOCIATED(VarFnorm)) THEN
    CALL FATAL(SolverName,'Variable "Time" not found')
  END IF
  Time => VarTime % Values
  CurrentTime = Time(1)

  ! time, pos, normal, tang
  !PRINT *, CurrentTime, datafield(1,1), datafield(1,2), datafieldentries
  I=2
  DO
    !PRINT *, I, CurrentTime, datafield(1,i)
    IF (((CurrentTime < datafield(1,i)) .AND. (CurrentTime >  datafield(1,i-1))) .AND. (I <= datafieldentries)) EXIT
    I = I +1
  END DO

  !WRITE(Message,'(a,es6.3,a,d,es6.3,es6.3)') "time", Currenttime, "time intervall", I, datafield(1,I), datafield(1,I-1)
  WRITE(Message,*) "time", Currenttime, "time intervall", I, datafield(1,I), datafield(1,I-1)
  CALL INFO(SolverName, Message, Level=3)
  !WRITE(Message,'(a,es6.3,a,es6.3,es6.3)') "position", datafield(2,I), "forces", datafield(3,I), datafield(4,I)
  WRITE(Message,*) "position", datafield(2,I), "forces", datafield(3,I), datafield(4,I)
  CALL INFO(SolverName, Message, Level=3)

  Fnorm(1) =  datafield(4,I)
  Ftang(1) =  datafield(3,I)
  Zcoord(1) = datafield(2,I)
  
CONTAINS
  
  SUBROUTINE read_fc(tzf,n_lines)

    USE DefUtils
    IMPLICIT NONE

    REAL(KIND=dp), POINTER, INTENT(out):: tzf(:,:)
    REAL(KIND=dp), ALLOCATABLE:: time_forces(:,:), time_zc(:,:)
    INTEGER :: io,io2, stat, i, n_lines, n_lines_z
    CHARACTER(len=512) :: msg
    REAL(KIND=dp) :: a

    ! We want to read both Fc.dat (header lines= 2) and Mvt_FrameBody_TopLeftCorner.dat (header lines = 8)

    OPEN(unit = 67, file = "F_wall_iceb.dat", status = "old", action = "read", iostat=stat, iomsg=msg)
    n_lines=0

    DO 
      READ(67,*, iostat=io)
      IF (io/=0) EXIT
      n_lines = n_lines + 1
    END DO
    !PRINT *, n_lines
    CLOSE(67)

    n_lines = n_lines-2   ! header lines


    OPEN(unit = 67, file = "Mvt_Frame_body_TopLeftCorner.dat", status = "old", action = "read", iostat=stat, iomsg=msg)
    n_lines_z=0
    DO 
      READ(67,*, iostat=io)
      IF ((io/=0) .OR. (n_lines_z == n_lines)) EXIT
      n_lines_z = n_lines_z + 1
    END DO
    !PRINT *, n_lines_z
    CLOSE(67)

    n_lines_z = n_lines_z-8   ! header lines

    ALLOCATE(time_forces(2, n_lines))
    !ALLOCATE(time_zc(2,n_lines_z))
    !ALLOCATE(tzf(4,n_lines_z))
    ALLOCATE(time_zc(2,n_lines))
    ALLOCATE(tzf(4,n_lines))

    OPEN(newunit = io, file = "F_wall_iceb.dat", status = "old", action = "read", iostat=stat, iomsg=msg)
    OPEN(unit = 66, file = "Mvt_Frame_body_TopLeftCorner.dat", status = "old", action = "read", iostat=stat, iomsg=msg)

    DO i = 1,2
      READ(io, *)
    END DO

    DO i = 1,8
      READ(66, *)
    END DO


    DO i = 1, n_lines
      READ(io,*) a, time_forces(1,i), time_forces(2,i)  ! read time, Fcx (normal efforts), Fcz (tangential efforts)
    END DO

    DO i = 1, n_lines! n_lines_z
      READ(66,*) time_zc(1,i), a, time_zc(2,i)  ! read altitude of the contact,  header +1
    END DO

    CLOSE(io)
    CLOSE(66)


    DO i = 1, n_lines !n_lines_z
      tzf(1,i) = time_zc(1,i) 
      tzf(2,i) = time_zc(2,i)
    END DO

    DO i = 1, n_lines 
      tzf(3,i) = time_forces(1,i) 
      tzf(4,i) = time_forces(2,i)
    END DO


  END SUBROUTINE read_fc
END SUBROUTINE ReadIsisData
