#ifndef AMRSOLVER_BOUNDARY_H
#define AMRSOLVER_BOUNDARY_H

enum boundary_flags {

	NONE,
	DIRICHLET,
	NEUMANN,
	MPI_BOUNDARY
};

typedef boundary_flags BC;


#endif
