#include <iostream>
#include "boundary.h"
#include "octree.h"
#include "direction.h"

namespace myOctree {


void set_FieldBc_FieldBcVal(int number, std::string name, FieldBc **bc ) {

	


//	for (int i=0; i<3; i++) {
//	      	for (int j=0; j<2; j++) {
//			std::cerr << bc[i][j] ;
//		}
//	}

	//std::cerr << name ;
	//std::cerr << "\n" ;

	for (std::list<Octree*>::iterator it = nodes.begin(), end = nodes.end(); it != end; ++it) {	

	
		if((*it)->number==number) {

			std::cerr << "Block " <<  (*it)->number << std::endl;	


			for(int l = 0; l<scalar_fields.size() ; l++) {	
			
				if( (*it)->get_block_data()->scalarfields[l]->name == name ) {
			
					//std::cerr << "Setting BC to " << (*it)->get_block_data()->scalarfields[l]->name << std::endl;
					for (int i=0; i<3; i++) {
					      	for (int j=0; j<2; j++) {
							 (*it)->get_block_data()->scalarfields[l]->bc[i][j] = bc[i][j];

							 //std::cerr << (*it)->get_block_data()->scalarfields[l]->bc[i][j];
						}
					}
				}
			}		
		}
	}	
	//std::cerr << "\n\n" ;

}

void set_VecFieldBc_VecFieldBcVal(int number, std::string name, FieldBc **xbc, FieldBc **ybc, FieldBc **zbc ) {


//	for (int i=0; i<3; i++) {
//	      	for (int j=0; j<2; j++) {
//			std::cerr << xbc[i][j] ;
//		}
//	}
//	for (int i=0; i<3; i++) {
//	      	for (int j=0; j<2; j++) {
//			std::cerr << ybc[i][j] ;
//		}
//	}
//	for (int i=0; i<3; i++) {
//	      	for (int j=0; j<2; j++) {
//			std::cerr << zbc[i][j] ;
//		}
//	}

	//std::cerr << name ;
	//std::cerr << "\n" ;

	for (std::list<Octree*>::iterator it = nodes.begin(), end = nodes.end(); it != end; ++it) {	

	
		if((*it)->number==number) {

			std::cerr << "Block " <<  (*it)->number << std::endl;	


			for(int l = 0; l<vector_fields.size() ; l++) {	
			
				if( (*it)->get_block_data()->vectorfields[l]->name == name ) {
			
					//std::cerr << "Setting BC to " << (*it)->get_block_data()->vectorfields[l]->name << std::endl;
					for (int i=0; i<3; i++) {
					      	for (int j=0; j<2; j++) {
							 (*it)->get_block_data()->vectorfields[l]->xbc[i][j] = xbc[i][j];
							 (*it)->get_block_data()->vectorfields[l]->ybc[i][j] = ybc[i][j];
							 (*it)->get_block_data()->vectorfields[l]->zbc[i][j] = zbc[i][j];

							 //std::cerr << (*it)->get_block_data()->scalarfields[l]->bc[i][j];
						}
					}
				}
			}		
		}
	}	
	//std::cerr << "\n\n" ;

}


}

