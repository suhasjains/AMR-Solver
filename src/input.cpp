#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "boundary.h"
#include "octreegrid.h"

namespace std {

BC string_to_BC(string bc) {

	BC bcc;

	if(bc=="NONE")
		bcc = NONE;

	if(bc=="DIRICHLET")
		bcc = DIRICHLET;

	if(bc=="NEUMANN")
		bcc = NEUMANN;

	if(bc=="MPI_BOUNDARY")
		bcc = MPI_BOUNDARY;
		
	return bcc;
}

void read_input_file() {


	string line;
	int blocknumber, level;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	string eastbc, westbc, northbc, southbc, topbc, bottombc;
	BC east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc;

	ifstream file ("../input/input.pfs"); 
	if(file.fail()) {
		cerr << "Error opening input file\n";
		exit(1);
	}

	if(file.peek() == ifstream::traits_type::eof()) {
		cerr << "Nothing in the input file\n";
		exit(1);
	}

	else {

		getline(file,line);
		if(line=="blocks") {
	
			while(file) {
				file >> blocknumber >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
				file >> level;
				file >> eastbc >> westbc >> northbc >> southbc >> topbc >> bottombc;
							
				east_bc = string_to_BC(eastbc);		
				west_bc = string_to_BC(westbc);		
				north_bc = string_to_BC(northbc);		
				south_bc = string_to_BC(southbc);		
				top_bc = string_to_BC(topbc);		
				bottom_bc = string_to_BC(bottombc);		
				
				if(file.eof() ) break;	
				
				myOctree::create_node(xmin, xmax, ymin, ymax, zmin, zmax, level, east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc);

			}

		}

	}
	file.close();
}

}
