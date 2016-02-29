#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "boundary.h"
#include "octreegrid.h"
#include "adapt.h"

namespace std {

myOctree::NodeBc string_to_NodeBc(string bc) {

	myOctree::NodeBc bcc;

	if(bc=="N")
		bcc = myOctree::NONE;

	if(bc=="B")
		bcc = myOctree::BOUNDARY;

	if(bc=="MPI_BOUNDARY")
		bcc = myOctree::MPI_BOUNDARY;
		
	return bcc;
}

myOctree::FieldBc string_to_FieldBc(string bc) {

	myOctree::FieldBc bcc;

	if(bc=="N")
		bcc = myOctree::none;

	if(bc=="D")
		bcc = myOctree::dirichlet;
	
	if(bc=="NE")
		bcc = myOctree::neumann;

	if(bc=="MB")
		bcc = myOctree::mpi_boundary;
		
	return bcc;
}

int read_blocks(ifstream& file) {
	
	string line, str;

	int blocknumber, level;	
	double xmin, xmax, ymin, ymax, zmin, zmax;
	string eastbc, westbc, northbc, southbc, topbc, bottombc;
	myOctree::NodeBc east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc;

	//skips line if line is empty
	getline(file,line);
	while(line.empty()) {
		getline(file, line);
	}


	if(line=="blocks") {


		file >> str;
		if(str != "{") {
			cerr << "Expected an opening bracket" <<endl;				
			exit(1);	
		}
		
		while(file) {
			if(file >> str && str == "}") break;
			else	blocknumber = atoi(str.c_str());

			//file >> blocknumber;
			file >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
			file >> level;
			file >> eastbc >> westbc >> northbc >> southbc >> topbc >> bottombc;
						
			east_bc = string_to_NodeBc(eastbc);		
			west_bc = string_to_NodeBc(westbc);		
			north_bc = string_to_NodeBc(northbc);		
			south_bc = string_to_NodeBc(southbc);		
			top_bc = string_to_NodeBc(topbc);		
			bottom_bc = string_to_NodeBc(bottombc);		

			//cerr << blocknumber << xmin << xmax << ymin << ymax << zmin << zmax << endl;

			myOctree::create_node(blocknumber, xmin, xmax, ymin, ymax, zmin, zmax, level, east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc);

		}

	}

	return blocknumber;
}

void read_scalar_fields(ifstream& file, int number) {


	int blocknumber;	
	string line, str;
	string eastbc, westbc, northbc, southbc, topbc, bottombc;
	myOctree::FieldBc east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc;
	double eastbcval, westbcval, northbcval, southbcval, topbcval, bottombcval;

	//Skips if a line is empty		
	getline(file, line);
	while(line.empty()) {
		getline(file, line);
	}

	if(line=="scalar_fields") {

		file >> str;
                if(str != "{") {
                        cerr << "Expected an opening bracket" <<endl;
                        exit(1);
                }
		
		cerr << "\nUser defined scalar fields" << endl;

		while(file) { 
	
		
			if(file >> str && str == "}") break;
                        else    myOctree::scalar_fields.push_back(str);
			cerr << str << endl;

			for(int i = 0; i < number ; ++i) {
		
				file >> blocknumber; 
				file >> eastbc >> westbc >> northbc >> southbc >> topbc >> bottombc;
				file >> eastbcval >> westbcval >> northbcval >> southbcval >> topbcval >> bottombcval;
				
				east_bc = string_to_FieldBc(eastbc);		
				west_bc = string_to_FieldBc(westbc);		
				north_bc = string_to_FieldBc(northbc);		
				south_bc = string_to_FieldBc(southbc);		
				top_bc = string_to_FieldBc(topbc);		
				bottom_bc = string_to_FieldBc(bottombc);		

				cerr << blocknumber << eastbc << westbc << northbc << southbc << topbc << bottombc << endl;

				/*add a function here which sets boundary condition and boundary values*/
			}
		}	
	}
}


void read_vector_fields(ifstream& file) {

	string line, str;
	
	//Skips if a line is empty		
	getline(file, line);
	while(line.empty()) {
		getline(file, line);
	}

	if(line=="vector_fields") {

		file >> str;
	        if(str != "{") {
	                cerr << "Expected an opening bracket" <<endl;
	                exit(1);
	        }
		
		cerr << "\nUser defined vector Fields" << endl;

		while(file) { 
	
		
			if(file >> str && str == "}") break;
	                else    myOctree::vector_fields.push_back(str);
			cerr << str << endl;	
		}	
	}

}

void read_max_level(ifstream& file) {

	string line, str;
	
	//Skips if a line is empty		
	getline(file, line);
	while(line.empty()) {
		getline(file, line);
	}

	if(line=="max_level") {

		file >> myOctree::max_level;
	
	}

}
		

void read_input_file() {

	int blocknumber;
	string line, str;

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

	
		blocknumber = read_blocks(file);

		read_scalar_fields(file, blocknumber);

		read_vector_fields(file);

		read_max_level(file);

	}

	file.close();
}

}
