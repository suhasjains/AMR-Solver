#ifndef AMRSOLVER_BOUNDARY_H
#define AMRSOLVER_BOUNDARY_H
#include <string>

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




void set_FieldBc_FieldBcVal(int, std::string, FieldBc**, double** ); 
void set_VecFieldBc_VecFieldBcVal(int, std::string, FieldBc**, FieldBc**, FieldBc**, double**, double**, double** ); 

}
#endif
