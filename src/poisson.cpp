#include <stdlib.h>
#include <string>
#include "octreegrid.h"
#include "poisson.h"
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

	long int loop = 0;

	//iteration loop
	while(loop <250 ) {
	//while(global_res > pow(10,-7)) {

		global_res = 0.0;
		
		//exchange ghost values
		myOctree::exchange_ghost_val(level, name);

		//loop over all the nodes
		for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {
	
			for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {
	
				Field *f = (*it)->get_block_data()->scalarfields[l];
							//	std::cerr << f->name << std::endl;
	
				if( f->name == name ) {
					
					global_res = jacobi_for_field(*it, f, global_res);	

				}	
			}
		}

		std::cout << "Residual " << global_res << std::endl;

		loop++;

	}	
}

/*!Solves poisson equation using multi-block red-black gauss-seidel procedure at the given level for the given scalar field.*/
void gauss_seidel(int level, std::string name) {

	myOctree::create_lists_of_level_nodes();

	if(myOctree::level_nodes[level].empty()) {
        	std::cerr << "No blocks in level " << level << std::endl;
	 	exit(1);
	}

	std::cerr << "Solving" << std::endl;

	double global_res = 1.0;

	long int loop = 0;

	//iteration loop
	//while(loop <250 ) {
	while(global_res > pow(10,-7)) {

		global_res = 0.0;
		
		//exchange ghost values
		myOctree::exchange_ghost_val(level, name);

		//loop over all the nodes
		for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {
	
			for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {
	
				Field *f = (*it)->get_block_data()->scalarfields[l];
				if( f->name == name ) {
					
					global_res = gauss_seidel_red(*it, f, global_res);	

				}	
			}
		}

		//exchange ghost values
		myOctree::exchange_ghost_val(level, name);
		
		//loop over all the nodes
		for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {
	
			for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {
	
				Field *f = (*it)->get_block_data()->scalarfields[l];
				if( f->name == name ) {
					
					global_res = gauss_seidel_black(*it, f, global_res);	

				}	
			}
		}

		//std::cout << "Residual " << global_res << std::endl; 

		loop++;

	}	
		std::cout << "Residual " << global_res << std::endl; 
		std::cout << "No of loops " << loop << std::endl;
}

double jacobi_for_field(Octree* node, Field* f, double global_res) {

	double res, err;
	double dx, dy, dz;
	double ue, uw, un, us, ut, ub, up;
	double force = 1.0;
	res = 0.0;
	
	Field temp(nx_block+2*pad,ny_block+2*pad,nz_block+2*pad, "temp");
	//temp.set_field(0.0);	
	
	dx = node->get_block_data()->dx;
	dy = node->get_block_data()->dy;
	dz = node->get_block_data()->dz;
	
	
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
				
	
				if(i==pad)		temp.val[i][j][k] += (2.0*uw - up)/pow(dx,2); 	
				else			temp.val[i][j][k] += uw/pow(dx,2);	
				if(i==nx_block+pad-1)	temp.val[i][j][k] += (2.0*ue - up)/pow(dx,2); 	
				else			temp.val[i][j][k] += ue/pow(dx,2);	
				if(j==pad)		temp.val[i][j][k] += (2.0*us - up)/pow(dy,2); 	
				else			temp.val[i][j][k] += us/pow(dy,2);	
				if(j==ny_block+pad-1)	temp.val[i][j][k] += (2.0*un - up)/pow(dy,2); 	
				else			temp.val[i][j][k] += un/pow(dy,2);	
				if(k==pad)		temp.val[i][j][k] += (2.0*ub - up)/pow(dz,2); 	
				else			temp.val[i][j][k] += ub/pow(dz,2);	
				if(k==nz_block+pad-1)	temp.val[i][j][k] += (2.0*ut - up)/pow(dz,2); 	
				else			temp.val[i][j][k] += ut/pow(dz,2);	
	
				temp.val[i][j][k] -= force;
				temp.val[i][j][k] /= 2.0;
				temp.val[i][j][k] /= (1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
	
			//	std::cerr << temp.val[i][j][k] << std::endl;
	
				/*magnitude to be considered*/
				err = std::abs(temp.val[i][j][k] - f->val[i][j][k]);
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

	return global_res;

}

/*!Solving grid points marked red in gauss-seidel algorithm.*/
double gauss_seidel_red(Octree* node, Field* f, double global_res) {

	double res, err;
	double dx, dy, dz;
	double ue, uw, un, us, ut, ub, up;
	double prev_val, curr_val;
	double force = 100.0;
	double relax = 1.5;
	res = 0.0;
	
	Field temp(nx_block+2*pad,ny_block+2*pad,nz_block+2*pad, "temp");
	//temp.set_field(0.0);	
	
	dx = node->get_block_data()->dx;
	dy = node->get_block_data()->dy;
	dz = node->get_block_data()->dz;
	
	
	//jocobi work - add here
	for(int i=pad; i<nx_block+pad; i++) {
		for(int j=pad; j<ny_block+pad; j++) {
			for(int k=pad; k<nz_block+pad; k++) {
			

				if((i+j+k)%2==0) {

					ue = f->val[i+1][j][k];		
					uw = f->val[i-1][j][k];		
					un = f->val[i][j+1][k];		
					us = f->val[i][j-1][k];		
					ut = f->val[i][j][k+1];		
					ub = f->val[i][j][k-1];		
					up = f->val[i][j][k];		
					
	
					if(i==pad)		temp.val[i][j][k] += (2.0*uw - up)/pow(dx,2); 	
					else			temp.val[i][j][k] += uw/pow(dx,2);	
					if(i==nx_block+pad-1)	temp.val[i][j][k] += (2.0*ue - up)/pow(dx,2); 	
					else			temp.val[i][j][k] += ue/pow(dx,2);	
					if(j==pad)		temp.val[i][j][k] += (2.0*us - up)/pow(dy,2); 	
					else			temp.val[i][j][k] += us/pow(dy,2);	
					if(j==ny_block+pad-1)	temp.val[i][j][k] += (2.0*un - up)/pow(dy,2); 	
					else			temp.val[i][j][k] += un/pow(dy,2);	
					if(k==pad)		temp.val[i][j][k] += (2.0*ub - up)/pow(dz,2); 	
					else			temp.val[i][j][k] += ub/pow(dz,2);	
					if(k==nz_block+pad-1)	temp.val[i][j][k] += (2.0*ut - up)/pow(dz,2); 	
					else			temp.val[i][j][k] += ut/pow(dz,2);	
	
					temp.val[i][j][k] -= force;
					temp.val[i][j][k] /= 2.0;
					temp.val[i][j][k] /= (1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
	
					prev_val = f->val[i][j][k];

					f->val[i][j][k] = (1.0 - relax) * f->val[i][j][k] + relax * temp.val[i][j][k];
					
					curr_val = f->val[i][j][k];
	
					/*magnitude to be considered*/
					err = std::abs(curr_val - prev_val);
					if(err>res) 	res = err;	
				
				}	
	
			}
		}
	}	
	
	
	if(res>global_res)	global_res = res;

	return global_res;

}

/*!Solving grid points marked black in gauss-seidel algorithm.*/
double gauss_seidel_black(Octree* node, Field* f, double global_res) {

	double res, err;
	double dx, dy, dz;
	double ue, uw, un, us, ut, ub, up;
	double prev_val, curr_val;
	double force = 100.0;
	double relax = 1.5;
	res = 0.0;
	
	Field temp(nx_block+2*pad,ny_block+2*pad,nz_block+2*pad, "temp");
	//temp.set_field(0.0);	
	
	dx = node->get_block_data()->dx;
	dy = node->get_block_data()->dy;
	dz = node->get_block_data()->dz;
	
	
	//jocobi work - add here
	for(int i=pad; i<nx_block+pad; i++) {
		for(int j=pad; j<ny_block+pad; j++) {
			for(int k=pad; k<nz_block+pad; k++) {
			
				if((i+j+k)%2!=0) {

					ue = f->val[i+1][j][k];		
					uw = f->val[i-1][j][k];		
					un = f->val[i][j+1][k];		
					us = f->val[i][j-1][k];		
					ut = f->val[i][j][k+1];		
					ub = f->val[i][j][k-1];		
					up = f->val[i][j][k];		
					
	
					if(i==pad)		temp.val[i][j][k] += (2.0*uw - up)/pow(dx,2); 	
					else			temp.val[i][j][k] += uw/pow(dx,2);	
					if(i==nx_block+pad-1)	temp.val[i][j][k] += (2.0*ue - up)/pow(dx,2); 	
					else			temp.val[i][j][k] += ue/pow(dx,2);	
					if(j==pad)		temp.val[i][j][k] += (2.0*us - up)/pow(dy,2); 	
					else			temp.val[i][j][k] += us/pow(dy,2);	
					if(j==ny_block+pad-1)	temp.val[i][j][k] += (2.0*un - up)/pow(dy,2); 	
					else			temp.val[i][j][k] += un/pow(dy,2);	
					if(k==pad)		temp.val[i][j][k] += (2.0*ub - up)/pow(dz,2); 	
					else			temp.val[i][j][k] += ub/pow(dz,2);	
					if(k==nz_block+pad-1)	temp.val[i][j][k] += (2.0*ut - up)/pow(dz,2); 	
					else			temp.val[i][j][k] += ut/pow(dz,2);	
	
					temp.val[i][j][k] -= force;
					temp.val[i][j][k] /= 2.0;
					temp.val[i][j][k] /= (1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
	
					prev_val = f->val[i][j][k];
					
					f->val[i][j][k] = (1.0 - relax) * f->val[i][j][k] + relax * temp.val[i][j][k];
					
					curr_val = f->val[i][j][k];
	
					/*magnitude to be considered*/
					err = std::abs(curr_val - prev_val);
					if(err>res) 	res = err;	
	
				}	
			}
		}
	}	
	
	
	if(res>global_res)	global_res = res;

	return global_res;

}


}	
