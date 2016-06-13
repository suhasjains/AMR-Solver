#include <stdlib.h>
#include <string>
#include "octreegrid.h"
#include <iostream>
//#include <cmath>
#include "multigrid.h"

using myOctree::Field;
using myOctree::VecField;
using myOctree::Octree;
using myOctree::pad;
using myOctree::level_nodes;


namespace amrsolver {

/*!Bilinear interpolation of the given field's ghost cells from level-1 to level.*/
void prolongate_ghost_val_to(int level, std::string name) {

	//relative locations wrt siblings
	int *loc_x, *loc_y, *loc_z;
	loc_x = new int;
	loc_y = new int;
	loc_z = new int;

	myOctree::create_lists_of_level_nodes();
	
	if(myOctree::level_nodes[level].empty()) {
                std::cerr << "No blocks in level " << level << std::endl;
                exit(1);
        }

	//loop over nodes
	for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {

        	//get location relative to siblings
            	(*it)->get_relative_location(loc_x,loc_y,loc_z);

            	//std::cerr << *loc_x << " " << *loc_y << " " << *loc_z << std::endl;

            	//loop over scalar fields
            	for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {

       			Field *f = (*it)->get_block_data()->scalarfields[l];

               		if( f->name == name ) {

  				for(int i=0; i<f->Nx; i++) {
         				for(int j=0; j<f->Ny; j++) {
                         			for(int k=0; k<f->Nz; k++) {
									
										//replace relative location with neighbour == none
							if(((*it)->get_block_data()->flag[i][j][k]==EAST_GHOST && *loc_x==0)||
							   ((*it)->get_block_data()->flag[i][j][k]==WEST_GHOST && *loc_x==1)||
							   ((*it)->get_block_data()->flag[i][j][k]==NORTH_GHOST && *loc_y==0)||
							   ((*it)->get_block_data()->flag[i][j][k]==SOUTH_GHOST && *loc_y==1)||
							   ((*it)->get_block_data()->flag[i][j][k]==TOP_GHOST && *loc_z==0)||
							   ((*it)->get_block_data()->flag[i][j][k]==BOTTOM_GHOST && *loc_z==1)) {
								if((i-pad)%2==0) {
									if((j-pad)%2==0) {
										if((k-pad)%2==0) {
											(i+pad)/2 - 1
											(i+pad)/2
											(j+pad)/2 - 1
											(j+pad)/2
											(k+pad)/2 - 1
											(k+pad)/2
										}
										else if((k-pad)%2!=0) {
											(i+pad)/2 - 1
											(i+pad)/2
											(j+pad)/2 - 1
											(j+pad)/2
											(k+pad)/2 
											(k+pad)/2 + 1
										}
									}
									else if((j-pad)%2!=0) {
										if((k-pad)%2==0) {
											(i+pad)/2 - 1
											(i+pad)/2
											(j+pad)/2 
											(j+pad)/2 + 1
											(k+pad)/2 - 1
											(k+pad)/2
										}
										else if((k-pad)%2!=0) {
											(i+pad)/2 - 1
											(i+pad)/2
											(j+pad)/2 
											(j+pad)/2 + 1
											(k+pad)/2 
											(k+pad)/2 + 1
										}
									}
								}
								else if((i-pad)%2!=0) {
									if((j-pad)%2==0) {
										if((k-pad)%2==0) {
											(i+pad)/2 
											(i+pad)/2 + 1
											(j+pad)/2 - 1
											(j+pad)/2
											(k+pad)/2 - 1
											(k+pad)/2
										}
										else if((k-pad)%2!=0) {
											(i+pad)/2 
											(i+pad)/2 + 1
											(j+pad)/2 - 1
											(j+pad)/2
											(k+pad)/2 
											(k+pad)/2 + 1
										}
									}
									else if((j-pad)%2!=0) {
										if((k-pad)%2==0) {
											(i+pad)/2 
											(i+pad)/2 + 1
											(j+pad)/2 
											(j+pad)/2 + 1
											(k+pad)/2 - 1
											(k+pad)/2
										}
										else if((k-pad)%2!=0) {
											(i+pad)/2 
											(i+pad)/2 + 1
											(j+pad)/2 
											(j+pad)/2 + 1
											(k+pad)/2 
											(k+pad)/2 + 1
										}
									}
								}
							}
							if(((*it)->get_block_data()->flag[i][j][k]==WEST_GHOST)&&()) {
								x transverse;
								y along;
								z along;
							}		
							if((*it)->get_block_data()->flag[i][j][k]==NORTH_GHOST) {
								x along;
								y transverse;
								z along;
							}		
							if((*it)->get_block_data()->flag[i][j][k]==SOUTH_GHOST) {
								x along;
								y transverse;
								z along;
							}		
							if((*it)->get_block_data()->flag[i][j][k]==TOP_GHOST) {
								x along;
								y along;
								z transverse;
							}		
							if((*it)->get_block_data()->flag[i][j][k]==BOTTOM_GHOST) {
								x along;
								y along;
								z transverse;
							}		

                           			}
          				}
                 		}
       			}
         	}
        }

}

}
