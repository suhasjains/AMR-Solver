#include <stdlib.h>
#include <string>
#include "octreegrid.h"
#include "ghost.h"
#include <iostream>
#include <cmath>

using myOctree::Field;
using myOctree::VecField;
using myOctree::Octree;

using myOctree::nx_block;
using myOctree::ny_block;
using myOctree::nz_block;
using myOctree::pad;

using myOctree::level_nodes;

namespace amrsolver {

//
void jacobi(int level, std::string name) {

	myOctree::create_lists_of_level_nodes();

	if(myOctree::level_nodes[level].empty()) {
        	std::cerr << "No blocks in level " << level << std::endl;
	 	exit(1);
	}

	double global_res = 1.0;
	double res, err;

	double dx, dy, dz;
	double ue, uw, un, us, ut, ub, up;

	double force = 0.0;

	//temporary field
	Field temp(nx_block+2*pad,ny_block+2*pad,nz_block+2*pad, "temp");

	//iteration loop
	while(global_res > pow(10,-7)) {

		global_res = 0.0;

		//loop over all the nodes
		for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {
	
			for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {
	
				Field *f = (*it)->get_block_data()->scalarfields[l];
	
				if( f->name == name ) {
					
					
					res = 0.0;

					temp.set_field(0.0);	

					dx = (*it)->get_block_data()->dx;
					dy = (*it)->get_block_data()->dy;
					dz = (*it)->get_block_data()->dz;


					//jocobi work - add here
					for(int i=pad; i<nx_block+pad; i++) {
						for(int j=pad; j<ny_block+pad; j++) {
							for(int k=pad; k<nz_block+pad; k++) {
									
								ue = f->val[i+1][j][k];		
								uw = f->val[i-1][j][k];		
								un = f->val[i][j+1][k];		
								us = f->val[i][j-1][k];		
								ut = f->val[i][j][k+1];		
								ub = f->val[i][j][k-1];		
								up = f->val[i][j][k];		
								

								if(i==pad)		temp.val[i][j][k] += (2.0*uw - up)/(dx*dx); 	
								else			temp.val[i][j][k] += uw/(dx*dx);	
								if(i==nx_block+pad-1)	temp.val[i][j][k] += (2.0*ue - up)/(dx*dx); 	
								else			temp.val[i][j][k] += ue/(dx*dx);	
								if(j==pad)		temp.val[i][j][k] += (2.0*us - up)/(dy*dy); 	
								else			temp.val[i][j][k] += us/(dy*dy);	
								if(j==ny_block+pad-1)	temp.val[i][j][k] += (2.0*un - up)/(dy*dy); 	
								else			temp.val[i][j][k] += un/(dy*dy);	
								if(k==pad)		temp.val[i][j][k] += (2.0*ub - up)/(dz*dz); 	
								else			temp.val[i][j][k] += ub/(dz*dz);	
								if(k==nz_block+pad-1)	temp.val[i][j][k] += (2.0*ut - up)/(dz*dz); 	
								else			temp.val[i][j][k] += ut/(dz*dz);	

								temp.val[i][j][k] -= force;
								temp.val[i][j][k] /= 2.0*(1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz));


								/*magnitude to be considered*/
								err = temp.val[i][j][k] - f->val[i][j][k];
								if(err>res) 	res = err;	


							}
						}
					}	

					//Transferring values from temp to field
					for(int i=pad; i<nx_block+pad; i++) {
						for(int j=pad; j<ny_block+pad; j++) {
							for(int k=pad; k<nz_block+pad; k++) {
				
								f->val[i][j][k] = temp.val[i][j][k];	
					
							}
						}
					}	


					if(res>global_res)	global_res = res;


				}	
			}
		}

		//exchange ghost values
		myOctree::exchange_ghost_val(level, name);

		std::cerr << "Hi" << global_res << std::endl;

	}	
}

void jacobi_for_field(Field* f) {



}

}	
