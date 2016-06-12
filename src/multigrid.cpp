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

							if((*it)->get_block_data()->flag[i][j][k]==EAST_GHOST) {
								if(*loc_x == 0) {
									left child;
								}
								else if(*loc_x == 1) {
									right child;
								}
								if(*loc_y == 0) {
									left child;
								}
								else if(*loc_y == 1) {
									right child;
								}
								if(*loc_z == 0) {
									left child;
								}
								else if(*loc_z == 1) {
									right child;
								}
								x transverse;
								y along;
								z along;
							}		
							if((*it)->get_block_data()->flag[i][j][k]==WEST_GHOST) {
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
