#ifndef AMRSOLVER_BOUNDARY_H
#define AMRSOLVER_BOUNDARY_H

namespace myOctree {

enum node_boundary_flags {

	NONE,
	BOUNDARY,
	MPI_BOUNDARY
};

typedef node_boundary_flags NodeBc;

enum field_boundary_flags {

	none,
	dirichlet,
	neumann,
	mpi_boundary,

};

typedef field_boundary_flags FieldBc;

}
#endif
