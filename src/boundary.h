#ifndef AMRSOLVER_BOUNDARY_H
#define AMRSOLVER_BOUNDARY_H

enum boundary_flags {

	NONE,
	DIRICHLET,
	NEUMANN,
	MPI_BOUNDARY
};

typedef boundary_flags BC;

class BoundVal {

	double east_bc_value;
	double west_bc_value; 
	double north_bc_value; 
	double south_bc_value; 
	double top_bc_value; 
	double bottom_bc_value;

};

#endif
