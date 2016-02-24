#include <iostream>
#include <fstream>
#include "octree.h"
#include <stdlib.h>

using myOctree::Octree;
using myOctree::nodes;

namespace std {

string BC_to_string(BC bcc) {

	string bc;

	if(bcc==NONE) 
		bc = "NONE";

	if(bcc==DIRICHLET) 
		bc = "DIRICHLET";

	if(bcc==NEUMANN) 
		bc = "NEUMANN";

	if(bcc==MPI_BOUNDARY) 
		bc = "MPI_BOUNDARY";


	return bc;
}


void write_output_file() {

	int count = 0;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	int level;
	string eastbc, westbc, northbc, southbc, topbc, bottombc;
        BC east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc;



	ofstream file ("output.pfs");
	if(file.fail()) {
                cerr << "Error opening input file\n";
                exit(1);
        }
		
	file << "blocks" <<endl;
	
	for (std::list<Octree*>::iterator it = nodes.begin(), end = nodes.end(); it != end; ++it)     {
	count++;

	xmin = (*it)->x_min;	
	ymin = (*it)->y_min;	
	zmin = (*it)->z_min;	
	xmax = (*it)->x_max;	
	ymax = (*it)->y_max;	
	zmax = (*it)->z_max;	
	level = (*it)->get_level();
	east_bc = (*it)->east_bc;
	west_bc = (*it)->west_bc;
	north_bc = (*it)->north_bc;
	south_bc = (*it)->south_bc;
	top_bc = (*it)->top_bc;
	bottom_bc = (*it)->bottom_bc;
	


	eastbc = BC_to_string(east_bc);
	westbc = BC_to_string(west_bc);
	northbc = BC_to_string(north_bc);
	southbc = BC_to_string(south_bc);
	topbc = BC_to_string(top_bc);
	bottombc = BC_to_string(bottom_bc);

	file << count << " " << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << " " << level<< " " ;			 
	file << eastbc << " " << westbc << " " << northbc << " " << southbc << " " << topbc << " " << bottombc << endl;
		
	}

	file.close();
	
}

}
