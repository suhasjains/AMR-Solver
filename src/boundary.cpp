#include <iostream>
#include "boundary.h"
#include "octree.h"
#include "direction.h"

namespace myOctree {

/*!Sets boundary conditions to blocks and boundary condition values to ghost cells in the blocks of the given scalar field.*/
void set_FieldBc_FieldBcVal(int number, std::string name, FieldBc **bc, double **bcval ) {

	for (std::list<Octree*>::iterator it = nodes.begin(), end = nodes.end(); it != end; ++it) {	

	
		if((*it)->number==number) {

			std::cerr << "Block " <<  (*it)->number << std::endl;	


			for(int l = 0; l<scalar_fields.size() ; l++) {	

				Field *f = (*it)->get_block_data()->scalarfields[l];			

				if( f->name == name ) {
		
					//setting BC	
					//std::cerr << "Setting BC to " << (*it)->get_block_data()->scalarfields[l]->name << std::endl;
					for (int i=0; i<3; i++) {
					      	for (int j=0; j<2; j++) {
							 f->bc[i][j] = bc[i][j];

							 //std::cerr << (*it)->get_block_data()->scalarfields[l]->bc[i][j];
						}
					}
					

					//Setting BC values
					for(int i=0; i<f->Nx; i++) {
						for(int j=0; j<f->Ny; j++) {
							for(int k=0; k<f->Nz; k++) {
								if(i<pad)	f->val[i][j][k] = bcval[XDIR][LEFT]; 			
								if(i>=(nx_block + pad))	f->val[i][j][k] = bcval[XDIR][RIGHT]; 			
								if(j<pad)	f->val[i][j][k] = bcval[YDIR][LEFT]; 			
								if(j>=(ny_block + pad))	f->val[i][j][k] = bcval[YDIR][RIGHT]; 			
								if(k<pad)	f->val[i][j][k] = bcval[ZDIR][LEFT]; 			
								if(k>=(nz_block + pad))	f->val[i][j][k] = bcval[ZDIR][RIGHT]; 			
									
							}
						}
					}
				}
			}		
		}
	}	
}

/*!Sets boundary conditions to blocks and boundary condition values to ghost cells in the blocks of the given vector field.*/  
void set_VecFieldBc_VecFieldBcVal(int number, std::string name, FieldBc **xbc, FieldBc **ybc, FieldBc **zbc, double** xbcval, double** ybcval, double** zbcval ) {


	for (std::list<Octree*>::iterator it = nodes.begin(), end = nodes.end(); it != end; ++it) {	

	
		if((*it)->number==number) {

			std::cerr << "Block " <<  (*it)->number << std::endl;	


			for(int l = 0; l<vector_fields.size() ; l++) {	
		
				VecField *f = (*it)->get_block_data()->vectorfields[l];
	
				if( f->name == name ) {
			
					//std::cerr << "Setting BC to " << (*it)->get_block_data()->vectorfields[l]->name << std::endl;
					for (int i=0; i<3; i++) {
					      	for (int j=0; j<2; j++) {
							 f->xbc[i][j] = xbc[i][j];
							 f->ybc[i][j] = ybc[i][j];
							 f->zbc[i][j] = zbc[i][j];

							 //std::cerr << (*it)->get_block_data()->scalarfields[l]->bc[i][j];
						}
					}
					
					//Setting BC values
                                        for(int i=0; i<f->Nx; i++) {
                                                for(int j=0; j<f->Ny; j++) {
                                                        for(int k=0; k<f->Nz; k++) {
                                                                if(i<pad) {
								       f->x[i][j][k] = xbcval[XDIR][LEFT];
								       f->y[i][j][k] = ybcval[XDIR][LEFT];
								       f->z[i][j][k] = zbcval[XDIR][LEFT];
								}
								if(i>=(f->Nx + pad)) {
								    f->x[i][j][k] = xbcval[XDIR][RIGHT];
								    f->y[i][j][k] = ybcval[XDIR][RIGHT];
								    f->z[i][j][k] = zbcval[XDIR][RIGHT];
                                                                }
								if(j<pad) {
								       f->x[i][j][k] = xbcval[YDIR][LEFT];
								       f->y[i][j][k] = ybcval[YDIR][LEFT];
								       f->z[i][j][k] = zbcval[YDIR][LEFT];
                                                                }
								if(j>=(f->Ny + pad)) {
									f->x[i][j][k] = xbcval[YDIR][RIGHT];
									f->y[i][j][k] = ybcval[YDIR][RIGHT];
									f->z[i][j][k] = zbcval[YDIR][RIGHT];
                                                                }
								if(k<pad) {
								     	f->x[i][j][k] = xbcval[ZDIR][LEFT];
								     	f->y[i][j][k] = ybcval[ZDIR][LEFT];
								     	f->z[i][j][k] = zbcval[ZDIR][LEFT];
                                                                }
								if(k>=(f->Nz + pad)) {
									f->x[i][j][k] = xbcval[ZDIR][RIGHT];
									f->y[i][j][k] = ybcval[ZDIR][RIGHT];
									f->z[i][j][k] = zbcval[ZDIR][RIGHT];
								}
                                                        }
                                                }
                                        }
				}
			}		
		}
	}	
}


}

