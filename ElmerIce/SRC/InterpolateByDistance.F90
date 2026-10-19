!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  Authors: Thomas Zwinger
! *  Email:   Thomas.Zwinger@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 17 May 2002
! *  Added limiters: 09 Jan 2008
! *
! ****************************************************************************/

!-----------------------------------------------------------------------------
!>  Initialisation routine to set initial free surface height 
!>  when the surface normal is not aligned to Z axis (rotated FS)
!> \ingroup Solvers
!-----------------------------------------------------------------------------
SUBROUTINE InterpolateByDistance( Model,Solver,dt,TransientSimulation )
  USE DefUtils  
  USE ElementDescription
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt

  TYPE(ValueList_t), POINTER :: SolverParams, Constants
  TYPE(Variable_t), POINTER :: DistVar, ContactVar
  !TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t), POINTER :: CurrentElement
  
  INTEGER :: DIM, i, j, n, t
  INTEGER, POINTER ::  DistPerm(:)

  REAL(KIND=dp), POINTER :: Dist(:), Contact(:),PArray(:,:) => NULL()
  REAL(KIND=dp) :: NodeHolder(3), RefCoord(3), mindist, maxdist

  LOGICAL ::  Found, FirstTime=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN)  :: VariableName, SolverName

 

  SAVE DIM, FirstTime, maxdist, mindist, Contact
  
  SolverParams => GetSolverParams()
  Constants => Model % Constants
  
  IF (FirstTime) THEN
    DIM = CoordinateSystemDimension()
    SolverName = "InterpolateByDistance"
    ! Read in reference point coordinates
    PArray => ListGetConstRealArray( Constants,'Reference Coordinate', Found)
    IF ( .NOT.Found ) THEN
      WRITE(Message, '(a,a,a)') "Can't find >Reference Coordinate<. Using (0,0,0)"  
      CALL WARN(SolverName, Message)
    ELSE
      DO i=1,DIM
        RefCoord(i) = PArray(i,1)
      END DO
      WRITE(Message, *) "Found coordinates of reference point:", RefCoord(1:DIM)
      CALL INFO(SolverName,Message,Level=1)
    END IF
    mindist = ListGetConstReal( Constants,"Minimum Distance", UnfoundFatal=.TRUE.)
    maxdist = ListGetConstReal( Constants,"Maximum Distance", UnfoundFatal=.TRUE.)
    WRITE(Message,*) "Min/max distance:", mindist, maxdist
    CALL INFO(SolverName,Message,Level=1)
    ! For each point of an active element compute the distance
    !ALLOCATE(Contact(Model % Mesh % NumberOfNodes))
    !CALL VariableAdd(Solver % Mesh % Variables, Solver % Mesh, Solver, "In Contact",1,Contact,Distperm)
    ContactVar => VariableGet(Solver % Mesh % Variables, "In Contact")
    IF (.NOT.ASSOCIATED(ContactVar)) THEN
      CALL FATAL(SolverName,"Contact Variable not associated")
    END IF
    DistVar => Solver % Variable
    IF (.NOT.ASSOCIATED(DistVar)) THEN
      CALL FATAL(SolverName,"Variable not associated")
    END IF
    Dist => DistVar % Values
    DistPerm => DistVar % Perm
    Contact => ContactVar % Values
    Contact = 0.0_dp
    IF (.NOT.ASSOCIATED(DistPerm)) &
         CALL FATAL(SolverName,"Permutation not associated")
    !CALL DefaultInitialize()
    DO t=1,Solver % NumberOfActiveElements
      CurrentElement => GetActiveElement(t)
      n = GetElementNOFNodes()
      IF(GetElementFamily(CurrentElement) == 1) CYCLE
      DO i=1,N
        j =  CurrentElement % NodeIndexes(i)
        !PRINT *, j, DistPerm(j), DIM
        Dist(DistPerm(j)) = &
             ComputeDistance(Solver,Model,RefCoord,j,DIM)        
        IF ( (Dist(DistPerm(j)) >= mindist) .AND. (Dist(DistPerm(j)) <= maxdist) ) THEN
          Contact(DistPerm(j)) = 1.0_dp
        !  PRINT *, "HUH"
        END IF
      END DO
    END DO
    FirstTime = .FALSE.
  END IF

  
  
CONTAINS

  FUNCTION ComputeDistance(Solver,Model,RefCoord,NodeIndex,DIM) RESULT(Distance)

    IMPLICIT NONE
    TYPE(Model_t), TARGET  :: Model
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp) :: RefCoord(3), Distance
    INTEGER :: DIM
    INTEGER :: NodeIndex
    INTEGER :: i
    REAL(KIND=dp) :: aux, Coord(3), RelCoord(3)

    
    Coord(1) = Solver % Mesh % Nodes % x(NodeIndex)
    Coord(2) = Solver % Mesh % Nodes % y(NodeIndex)
    IF (DIM == 3) &
         Coord(3) = Solver % Mesh % Nodes % z(NodeIndex)
    aux = 0.0_dp
    !PRINT *, "Coord=", Coord, "RefCoord=", RefCoord
    !PRINT *, "RelCoord=", RelCoord
    RelCoord(1:DIM) = Coord(1:DIM) - RefCoord(1:DIM)
    aux = SUM(RelCoord(1:DIM-1)*RelCoord(1:DIM-1))
    !DO i = 1,DIM
    !  aux = aux + (Coord(i) - RefCoord(i) )**2.0_dp
    !END DO
    Distance = SQRT(aux)
  END FUNCTION ComputeDistance
  
END SUBROUTINE InterpolateByDistance
