#include <iostream>
#include "boundary.h"
#include "octree.h"
#include "direction.h"

namespace myOctree {


void set_FieldBc_FieldBcVal(int number, std::string name, FieldBc **bc ) {

	


	for (int i=0; i<3; i++) {
	      	for (int j=0; j<2; j++) {
			std::cerr << bc[i][j] ;
		}
	}

	//std::cerr << name ;
	std::cerr << "\n" ;

	for (std::list<Octree*>::iterator it = nodes.begin(), end = nodes.end(); it != end; ++it) {	

	
		if((*it)->number==number) {

			//std::cerr << " Yes found it " <<  (*it)->number << std::endl;	


			for(int i = 0; i<scalar_fields.size() ; i++) {	
			
				if( (*it)->get_block_data()->scalarfields[i]->name == name ) {
			
					//std::cerr << "Yes found it " << (*it)->get_block_data()->scalarfields[i]->name << std::endl;
					for (int i=0; i<3; i++) {
					      	for (int j=0; j<2; j++) {
							 //(*it)->get_block_data()->scalarfields[i]->bc[i][j] = bc[i][j];

							 //std::cerr << (*it)->get_block_data()->scalarfields[i]->bc[i][j];
						}
					}
				}
			}		
		}
	}	
	std::cerr << "\n" ;

}

}

