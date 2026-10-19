This folder contains the material related to Simulation Case 9.2 on the modeling of Bowdoin glacier mechanical response to an iceberg capsize.

Two solvers are used: 
	- ISIS-CFD computing the water flow field around a rigid (idealized) iceberg. The output from this simulation are 
			(1) the magnitude of the contact force applied by the iceberg on the glacier front (Fc.dat) with the position of the iceberg edge that is in contact with the glacier Mvt_FrameBody_TopLeftCorner.dat. 
			(2) the pressure field on the glacier front (series of .csv files).
	- Elmer/Ice computing the deformation of Bowdoin glacier while taking the output from ISIS-CFD as input. The force and pressure fields are applied as boundary conditions. 
