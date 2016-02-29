#include <iostream>
#include <fstream>
#include "octree.h"
#include <stdlib.h>

using myOctree::Octree;
using myOctree::nodes;

namespace std {

string NodeBc_to_string(myOctree::NodeBc bcc) {

	string bc;

	if(bcc==myOctree::NONE) 
		bc = "N";

	if(bcc==myOctree::BOUNDARY) 
		bc = "B";

	if(bcc==myOctree::MPI_BOUNDARY) 
		bc = "MB";


	return bc;
}


void write_output_file() {

	int count = 0;
	int blocknumber;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	int level;
	string eastbc, westbc, northbc, southbc, topbc, bottombc;
        myOctree::NodeBc east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc;



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
	
	blocknumber = (*it)->number;

	eastbc = NodeBc_to_string(east_bc);
	westbc = NodeBc_to_string(west_bc);
	northbc = NodeBc_to_string(north_bc);
	southbc = NodeBc_to_string(south_bc);
	topbc = NodeBc_to_string(top_bc);
	bottombc = NodeBc_to_string(bottom_bc);

	file << blocknumber << " " << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << " " << level<< " " ;			 
	file << eastbc << " " << westbc << " " << northbc << " " << southbc << " " << topbc << " " << bottombc << endl;
		
	}

	file.close();
	
}

}
