#include "octreegrid.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include "multigrid.h"
#include "direction.h"

using myOctree::Field;
using myOctree::Block;
using myOctree::VecField;
using myOctree::Octree;
using myOctree::level_nodes;
using myOctree::multilevel_nodes;
using myOctree::pad;
using myOctree::nx_block;
using myOctree::ny_block;
using myOctree::nz_block;
using myOctree::scalar_fields;
using myOctree::vector_fields;

using myOctree::XDIR;
using myOctree::YDIR;
using myOctree::ZDIR;
using myOctree::RIGHT;
using myOctree::LEFT;


using myOctree::NONE;





namespace amrsolver {

/*!Exchanges ghost values of the given field at the given level*/	
void exchange_ghost_val(int level, std::string name) {

	myOctree::create_lists_of_level_nodes();

	if(level_nodes[level].empty()) {
		std::cerr << "Error!: No blocks to exchange ghost values in level " << level << std::endl; 
		exit(1);
	}

	for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {	

		for(int l = 0; l<scalar_fields.size() ; l++) {
			
			Field *f = (*it)->get_block_data()->scalarfields[l];

                        if( f->name == name ) {
			
				//exchanging ghost cells
				for (int m=0; m<3; m++) {
                        		for (int n=0; n<2; n++) {

						if((*it)->bc[m][n] == NONE && (*it)->neighbour[m][n] != NULL) {
				
										
							if(n==RIGHT && m==XDIR) {
						
								for(int i=nx_block+pad; i<f->Nx; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->val[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i-nx_block][j][k]; 
										}
									}
								}
							}
							if(n==LEFT && m==XDIR) {
						
								for(int i=0; i<pad; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->val[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i+nx_block][j][k]; 
										}
									}
								}	
							}
							if(n==RIGHT && m==YDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=ny_block+pad; j<f->Ny; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->val[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j-ny_block][k]; 
										}
									}
								}
							}
							if(n==LEFT && m==YDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=0; j<pad; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->val[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j+ny_block][k]; 
										}
									}
								}	
							}
							if(n==RIGHT && m==ZDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=nz_block+pad; k<f->Nz; k++) {
								
											f->val[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j][k-nz_block]; 
										}
									}
								}
							}
							if(n==LEFT && m==ZDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=0; k<pad; k++) {
								
											f->val[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j][k+nz_block]; 
										}
									}
								}	
							}
						}
                                   	}
                          	}
			}
		}
		
		for(int l = 0; l<vector_fields.size() ; l++) {
			
			VecField *f = (*it)->get_block_data()->vectorfields[l];

                        if( f->name == name ) {
			
				//exchanging ghost cells
				for (int m=0; m<3; m++) {
                        		for (int n=0; n<2; n++) {

						if((*it)->bc[m][n] == NONE && (*it)->neighbour[m][n] != NULL) {
										
						
							if(n==RIGHT && m==XDIR) {
						
								for(int i=nx_block+pad; i<f->Nx; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->x[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->x[i-nx_block][j][k]; 
											f->y[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->y[i-nx_block][j][k]; 
											f->z[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->z[i-nx_block][j][k]; 
										}
									}
								}
							}
							if(n==LEFT && m==XDIR) {
						
								for(int i=0; i<pad; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->x[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->x[i+nx_block][j][k]; 
											f->y[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->y[i+nx_block][j][k]; 
											f->z[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->z[i+nx_block][j][k]; 
										}
									}
								}	
							}
							if(n==RIGHT && m==YDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=ny_block+pad; j<f->Ny; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->x[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->x[i][j-ny_block][k]; 
											f->y[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->y[i][j-ny_block][k]; 
											f->z[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->z[i][j-ny_block][k]; 
										}
									}
								}
							}
							if(n==LEFT && m==YDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=0; j<pad; j++) {
                                                        			for(int k=pad; k<nz_block+pad; k++) {
								
											f->x[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->x[i][j+ny_block][k]; 
											f->y[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->y[i][j+ny_block][k]; 
											f->z[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->z[i][j+ny_block][k]; 
										}
									}
								}	
							}
							if(n==RIGHT && m==ZDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=nz_block+pad; k<f->Nz; k++) {
								
											f->x[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->x[i][j][k-nz_block]; 
											f->y[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->y[i][j][k-nz_block]; 
											f->z[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->z[i][j][k-nz_block]; 
										}
									}
								}
							}
							if(n==LEFT && m==ZDIR) {
						
								for(int i=pad; i<nx_block+pad; i++) {
                                                			for(int j=pad; j<ny_block+pad; j++) {
                                                        			for(int k=0; k<pad; k++) {
								
											f->x[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->x[i][j][k+nz_block]; 
											f->y[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->y[i][j][k+nz_block]; 
											f->z[i][j][k] = \
							 (*it)->neighbour[m][n]->get_block_data()->vectorfields[l]->z[i][j][k+nz_block]; 
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
}

/*!If this node is a left child, add 0, else if its a right child, add 1.
 
 */
void prolongate_multilevel_ghost_exchange_for_child(Octree* node, Octree* par_neigh, int l, int i, int j, int k, int iv, int jv, int kv, int addx, int addy, int addz) {
	//indices
	int i0, i1, j0, j1, k0, k1;	
	

	//grid and parent's neighbour's grid
	VecField* mesh = node->get_block_data()->mesh;
	VecField* meshp = par_neigh->get_block_data()->mesh; 

	//field and parent's neighbour's field
       	Field *f = node->get_block_data()->scalarfields[l];
       	Field *fp = par_neigh->get_block_data()->scalarfields[l];

	if((((iv-pad)%2==0)&&((iv-pad)!=0))||((iv-pad)==(f->Nx-1-2*pad))) {
		if((((jv-pad)%2==0)&&((jv-pad)!=0))||((jv-pad)==(f->Ny-1-2*pad))) {
			if((((kv-pad)%2==0)&&((kv-pad)!=0))||((kv-pad)==(f->Nz-1-2*pad))) {
				i0=addx+(iv+pad)/2 - 1;
				i1=addx+(iv+pad)/2;
				j0=addy+(jv+pad)/2 - 1;
				j1=addy+(jv+pad)/2;
				k0=addz+(kv+pad)/2 - 1;
				k1=addz+(kv+pad)/2;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
			else if((((kv-pad)%2!=0)&&((kv-pad)!=(f->Nz-1-2*pad)))||((kv-pad)==0)) {
				i0=addx+(iv+pad)/2 - 1;
				i1=addx+(iv+pad)/2;
				j0=addy+(jv+pad)/2 - 1;
				j1=addy+(jv+pad)/2;
				k0=addz+(kv+pad)/2 ;
				k1=addz+(kv+pad)/2 + 1;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
		}
		else if((((jv-pad)%2!=0)&&((jv-pad)!=(f->Ny-1-2*pad)))||((jv-pad)==0)) {
			if((((kv-pad)%2==0)&&((kv-pad)!=0))||((kv-pad)==(f->Nz-1-2*pad))) {
				i0=addx+(iv+pad)/2 - 1;
				i1=addx+(iv+pad)/2;
				j0=addy+(jv+pad)/2;
				j1=addy+(jv+pad)/2 + 1;
				k0=addz+(kv+pad)/2 - 1;
				k1=addz+(kv+pad)/2;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
			else if((((kv-pad)%2!=0)&&((kv-pad)!=(f->Nz-1-2*pad)))||((kv-pad)==0)) {
				i0=addx+(iv+pad)/2 - 1;
				i1=addx+(iv+pad)/2;
				j0=addy+(jv+pad)/2;
				j1=addy+(jv+pad)/2 + 1;
				k0=addz+(kv+pad)/2; 
				k1=addz+(kv+pad)/2 + 1;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
		}
	}
	else if((((iv-pad)%2!=0)&&((iv-pad)!=(f->Nx-1-2*pad)))||((iv-pad)==0)) {
		if((((jv-pad)%2==0)&&((jv-pad)!=0))||((jv-pad)==(f->Ny-1-2*pad))) {
			if((((kv-pad)%2==0)&&((kv-pad)!=0))||((kv-pad)==(f->Nz-1-2*pad))) {
				i0=addx+(iv+pad)/2; 
				i1=addx+(iv+pad)/2 + 1;
				j0=addy+(jv+pad)/2 - 1;
				j1=addy+(jv+pad)/2;
				k0=addz+(kv+pad)/2 - 1;
				k1=addz+(kv+pad)/2;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
			else if((((kv-pad)%2!=0)&&((kv-pad)!=(f->Nz-1-2*pad)))||((kv-pad)==0)) {
				i0=addx+(iv+pad)/2; 
				i1=addx+(iv+pad)/2 + 1;
				j0=addy+(jv+pad)/2 - 1;
				j1=addy+(jv+pad)/2;
				k0=addz+(kv+pad)/2; 
				k1=addz+(kv+pad)/2 + 1;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
		}
		else if((((jv-pad)%2!=0)&&((jv-pad)!=(f->Ny-1-2*pad)))||((jv-pad)==0)) {
			if((((kv-pad)%2==0)&&((kv-pad)!=0))||((kv-pad)==(f->Nz-1-2*pad))) {
				i0=addx+(iv+pad)/2; 
				i1=addx+(iv+pad)/2 + 1;
				j0=addy+(jv+pad)/2; 
				j1=addy+(jv+pad)/2 + 1;
				k0=addz+(kv+pad)/2 - 1;
				k1=addz+(kv+pad)/2;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
			else if((((kv-pad)%2!=0)&&((kv-pad)!=(f->Nz-1-2*pad)))||((kv-pad)==0)) {
				i0=addx+(iv+pad)/2; 
				i1=addx+(iv+pad)/2 + 1;
				j0=addy+(jv+pad)/2; 
				j1=addy+(jv+pad)/2 + 1;
				k0=addz+(kv+pad)/2; 
				k1=addz+(kv+pad)/2 + 1;
				f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
									meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
									meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
									fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
									fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
									fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
									fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
			}
		}
	}
}



/*!address, parent's neighbour's address, indices, indices of the virtual neighbour, relative location of virtual neighbour wrt siblings. */ 
void prolongate_for_multilevel_ghost_exchange_at(Octree* node, Octree* par_neigh, int l, int i, int j, int k, int iv, int jv, int kv, int loc_x, int loc_y, int loc_z) {

	Block* b =node->get_block_data(); 

	if(loc_x==1) {
    		if(loc_y==1) {
                  	if(loc_z==1) {
                       		prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,b->iNx/2,b->iNy/2,b->iNz/2);
      	           	}
         		else if(loc_z==0) {
                         	prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,b->iNx/2,b->iNy/2,0);
                 	}
          	}
                else if(loc_y==0) {
                	if(loc_z==1) {
                       		prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,b->iNx/2,0,b->iNz/2);
                   	}
                 	else if(loc_z==0) {
                        	prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,b->iNx/2,0,0);
                    	}
      		}
        }
	if(loc_x==0) {
            	if(loc_y==1) {
               		if(loc_z==1) {
                    		prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,0,b->iNy/2,b->iNz/2);
                   	}
          		else if(loc_z==0) {
                        	prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,0,b->iNy/2,0);
                  	}
            	}
    		else if(loc_y==0) {
         		if(loc_z==1) {
          			prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,0,0,b->iNz/2);
            		}
           		else if(loc_z==0) {
                  		prolongate_multilevel_ghost_exchange_for_child(node, par_neigh,l,i,j,k,iv, jv, kv,0,0,0);
           		}
    		}
        }
}

/*!Exchanges ghost values of the given field at the given multilevel*/	
void exchange_multilevel_ghost_val(int level, std::string name) {

	//relative locations wrt siblings
	int *loc_x, *loc_y, *loc_z;
	loc_x = new int;
	loc_y = new int;
	loc_z = new int;

	myOctree::create_lists_of_multilevel_nodes();

	if(multilevel_nodes[level].empty()) {
		std::cerr << "Error!: No blocks to exchange ghost values in multilevel " << level << std::endl; 
		exit(1);
	}

	for (std::list<Octree*>::iterator it = multilevel_nodes[level].begin(), end = multilevel_nodes[level].end(); it != end; ++it) {	

  		//get location relative to siblings
          	(*it)->get_relative_location(loc_x,loc_y,loc_z);

		for(int l = 0; l<scalar_fields.size() ; l++) {
			
			Field *f = (*it)->get_block_data()->scalarfields[l];

                        if( f->name == name ) {

				//traverse through all cells
				for (int i=0; i<f->Nx; i++) {
					for (int j=0; j<f->Ny; j++) {
						for (int k=0; k<f->Nz; k++) {
			
							//traverse through directions
							for (int m=0; m<3; m++) {
                        					for (int n=0; n<2; n++) {

									if((*it)->bc[m][n] == NONE && (*it)->neighbour[m][n] != NULL) {
										
										//exchange with neighbours
										if((*it)->neighbour[m][n]->get_level() == level) {
								
											if(n==RIGHT && m==XDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::EAST_GHOST) {
												f->val[i][j][k] = \
							 					(*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i-nx_block][j][k]; 
											}
											if(n==LEFT && m==XDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::WEST_GHOST) {
												f->val[i][j][k] = \
							 					(*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i+nx_block][j][k]; 
											}
											if(n==RIGHT && m==YDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::NORTH_GHOST) {
												f->val[i][j][k] = \
							 					(*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j-ny_block][k]; 
											}
											if(n==LEFT && m==YDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::SOUTH_GHOST) {
												f->val[i][j][k] = \
							 					(*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j+ny_block][k]; 
											}
											if(n==RIGHT && m==ZDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::TOP_GHOST) {
												f->val[i][j][k] = \
							 					(*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j][k-nz_block]; 
											}
											if(n==LEFT && m==ZDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::BOTTOM_GHOST) {
												f->val[i][j][k] = \
							 					(*it)->neighbour[m][n]->get_block_data()->scalarfields[l]->val[i][j][k+nz_block]; 
											}
										}

										//exchange with neighbour's child
										else {

											if(n==RIGHT && m==XDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::EAST_GHOST) {
												f->val[i][j][k] = restrict_at((*it)->neighbour[m][n], l, i-nx_block, j, k);
											}
											if(n==LEFT && m==XDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::WEST_GHOST) {
												f->val[i][j][k] = restrict_at((*it)->neighbour[m][n], l, i+nx_block, j, k);
											}
											if(n==RIGHT && m==YDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::NORTH_GHOST) {
												f->val[i][j][k] = restrict_at((*it)->neighbour[m][n], l, i, j-ny_block, k);
											}
											if(n==LEFT && m==YDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::SOUTH_GHOST) {
												f->val[i][j][k] = restrict_at((*it)->neighbour[m][n], l, i, j+ny_block, k);
											}
											if(n==RIGHT && m==ZDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::TOP_GHOST) {
												f->val[i][j][k] = restrict_at((*it)->neighbour[m][n], l, i, j, k-nz_block);
											}
											if(n==LEFT && m==ZDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::BOTTOM_GHOST) {
												f->val[i][j][k] = restrict_at((*it)->neighbour[m][n], l, i, j, k+nz_block);
											}

										}
									}

									//exchange with neighbour's parent
									else if((*it)->bc[m][n] == NONE && (*it)->neighbour[m][n] == NULL) {
										
										if(n==RIGHT && m==XDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::EAST_GHOST) {
											prolongate_for_multilevel_ghost_exchange_at((*it),(*it)->get_parent()->neighbour[m][n], l, \
												       i, j, k, i-nx_block, j, k, 0, *loc_y, *loc_z);
										}
										if(n==LEFT && m==XDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::WEST_GHOST) {
											prolongate_for_multilevel_ghost_exchange_at((*it),(*it)->get_parent()->neighbour[m][n], l, \
												       i, j, k, i+nx_block, j, k, 1, *loc_y, *loc_z);
										}
										if(n==RIGHT && m==YDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::NORTH_GHOST) {
											prolongate_for_multilevel_ghost_exchange_at((*it),(*it)->get_parent()->neighbour[m][n], l, \
												       i, j, k, i, j-ny_block, k, *loc_x, 0, *loc_z);
										}
										if(n==LEFT && m==YDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::SOUTH_GHOST) {
											prolongate_for_multilevel_ghost_exchange_at((*it),(*it)->get_parent()->neighbour[m][n], l, \
												       i, j, k, i, j+ny_block, k, *loc_x, 1, *loc_z);
										}
										if(n==RIGHT && m==ZDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::TOP_GHOST) {
											prolongate_for_multilevel_ghost_exchange_at((*it),(*it)->get_parent()->neighbour[m][n], l, \
												       i, j, k, i, j, k-nz_block, *loc_x, *loc_y, 0);
										}
										if(n==LEFT && m==ZDIR && (*it)->get_block_data()->flag[i][j][k]==myOctree::BOTTOM_GHOST) {
											prolongate_for_multilevel_ghost_exchange_at((*it),(*it)->get_parent()->neighbour[m][n], l, \
												       i, j, k, i, j, k+nz_block, *loc_x, *loc_y, 1);
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
}


}
