#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "boundary.h"
#include "octreegrid.h"
#include "adapt.h"

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


	string line, str;
	char c;
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
							
				east_bc = string_to_BC(eastbc);		
				west_bc = string_to_BC(westbc);		
				north_bc = string_to_BC(northbc);		
				south_bc = string_to_BC(southbc);		
				top_bc = string_to_BC(topbc);		
				bottom_bc = string_to_BC(bottombc);		

				//cerr << blocknumber << xmin << xmax << ymin << ymax << zmin << zmax << endl;

				myOctree::create_node(xmin, xmax, ymin, ymax, zmin, zmax, level, east_bc, west_bc, north_bc, south_bc, top_bc, bottom_bc);

			}

		}
		
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
			
			cerr << "\nUser defined scalar Fields" << endl;

			while(file) { 
		
			
				if(file >> str && str == "}") break;
        	                else    myOctree::scalar_fields.push_back(str);
				cerr << str << endl;	
			}	
		}
		
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

		//Skips if a line is empty		
		getline(file, line);
		while(line.empty()) {
			getline(file, line);
		}
	
		if(line=="max_level") {

			file >> myOctree::max_level;
		
		}
	
	}
	file.close();
}

}
