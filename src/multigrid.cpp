#include <stdlib.h>
#include <string>
#include "octreegrid.h"
#include <iostream>
//#include <cmath>
#include "direction.h"
#include "poisson.h"
#include "adapt.h"

using myOctree::Field;
using myOctree::VecField;
using myOctree::Octree;
using myOctree::Block;
using myOctree::pad;
using myOctree::level_nodes;


using myOctree::XDIR;
using myOctree::YDIR;
using myOctree::ZDIR;
using myOctree::RIGHT;
using myOctree::LEFT;

using myOctree::EAST_GHOST;
using myOctree::WEST_GHOST;
using myOctree::NORTH_GHOST;
using myOctree::SOUTH_GHOST;
using myOctree::TOP_GHOST;
using myOctree::BOTTOM_GHOST;
using myOctree::DOMAIN;



namespace amrsolver {



/*!Trilinear interpolation at given points

Input:

val0 @ (x0,y0,z1)
val1 @ (x0,y1,z1)
val2 @ (x1,y0,z1)
val3 @ (x1,y1,z1)
val4 @ (x0,y0,z0)
val5 @ (x0,y1,z0)
val6 @ (x1,y0,z0)
val7 @ (x1,y1,z0)

Output:

result @ (x2,y2,z2)

*/
double Trilinear_interpolate(double x0,double x1, double x2, double y0, double y1, double y2, double z0, double z1, double z2, \
				double val0, double val1, double val2, double val3, double val4, double val5, double val6, double val7) {

	double N0, N1, N2, N3, N4, N5, N6, N7;
	double result;

	N0 = ((x1-x2)*(y1-y2)*(z2-z0))/((x1-x0)*(y1-y0)*(z1-z0));
	N1 = ((x1-x2)*(y2-y0)*(z2-z0))/((x1-x0)*(y1-y0)*(z1-z0));
	N2 = ((x2-x0)*(y1-y2)*(z2-z0))/((x1-x0)*(y1-y0)*(z1-z0));
	N3 = ((x2-x0)*(y2-y0)*(z2-z0))/((x1-x0)*(y1-y0)*(z1-z0));
	N4 = ((x1-x2)*(y1-y2)*(z1-z2))/((x1-x0)*(y1-y0)*(z1-z0));
	N5 = ((x1-x2)*(y2-y0)*(z1-z2))/((x1-x0)*(y1-y0)*(z1-z0));
	N6 = ((x2-x0)*(y1-y2)*(z1-z2))/((x1-x0)*(y1-y0)*(z1-z0));
	N7 = ((x2-x0)*(y2-y0)*(z1-z2))/((x1-x0)*(y1-y0)*(z1-z0));

	result = val0*N0 + val1*N1 + val2*N2 + val3*N3 + val4*N4 + val5*N5 + val6*N6 + val7*N7;

	return result;	
}


/*!If this node is a left child, add 0, else if its a right child, add 1.
 
Parameters: octree node, field index, indices i j k, additions addx addy addz.
 */
void prolongate_ghost_for_child(Octree* node, int l, int i, int j, int k, int addx, int addy, int addz) {
	//indices
	int i0, i1, j0, j1, k0, k1;	
	

	//grid and parent's grid
	VecField* mesh = node->get_block_data()->mesh;
	VecField* meshp = node->get_parent()->get_block_data()->mesh; 

	//field and parent's field
       	Field *f = node->get_block_data()->scalarfields[l];
       	Field *fp = node->get_parent()->get_block_data()->scalarfields[l];


	if((node->get_block_data()->flag[i][j][k]==EAST_GHOST && node->neighbour[XDIR][RIGHT]==NULL)||
	   (node->get_block_data()->flag[i][j][k]==WEST_GHOST && node->neighbour[XDIR][LEFT]==NULL)||
	   (node->get_block_data()->flag[i][j][k]==NORTH_GHOST && node->neighbour[YDIR][RIGHT]==NULL)||
	   (node->get_block_data()->flag[i][j][k]==SOUTH_GHOST && node->neighbour[YDIR][LEFT]==NULL)||
	   (node->get_block_data()->flag[i][j][k]==TOP_GHOST && node->neighbour[ZDIR][RIGHT]==NULL)||
	   (node->get_block_data()->flag[i][j][k]==BOTTOM_GHOST && node->neighbour[ZDIR][LEFT]==NULL)) {
		if((((i-pad)%2==0)&&((i-pad)!=0))||((i-pad)==(f->Nx-1-2*pad))) {
			if((((j-pad)%2==0)&&((j-pad)!=0))||((j-pad)==(f->Ny-1-2*pad))) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2 ;
					k1=addz+(k+pad)/2 + 1;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
			}
			else if((((j-pad)%2!=0)&&((j-pad)!=(f->Ny-1-2*pad)))||((j-pad)==0)) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2;
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2;
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2; 
					k1=addz+(k+pad)/2 + 1;
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
		else if((((i-pad)%2!=0)&&((i-pad)!=(f->Nx-1-2*pad)))||((i-pad)==0)) {
			if((((j-pad)%2==0)&&((j-pad)!=0))||((j-pad)==(f->Ny-1-2*pad))) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2; 
					k1=addz+(k+pad)/2 + 1;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
			}
			else if((((j-pad)%2!=0)&&((j-pad)!=(f->Ny-1-2*pad)))||((j-pad)==0)) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2; 
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2; 
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2; 
					k1=addz+(k+pad)/2 + 1;
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
}

/*!If this node is a left child, add 0, else if its a right child, add 1.
 
Parameters: octree node, field index, indices i j k, additions addx addy addz.
 */
void prolongate_domain_for_child(Octree* node, int l, int i, int j, int k, int addx, int addy, int addz) {
	//indices
	int i0, i1, j0, j1, k0, k1;	
	

	//grid and parent's grid
	VecField* mesh = node->get_block_data()->mesh;
	VecField* meshp = node->get_parent()->get_block_data()->mesh; 

	//field and parent's field
       	Field *f = node->get_block_data()->scalarfields[l];
       	Field *fp = node->get_parent()->get_block_data()->scalarfields[l];

	//domain
	if(node->get_block_data()->flag[i][j][k]==DOMAIN) {
		//change conditions
//		if((i-pad)%2==0) {
//			if((j-pad)%2==0) {
//				if((k-pad)%2==0) {
		if((((i-pad)%2==0)&&((i-pad)!=0))||((i-pad)==(f->Nx-1-2*pad))) {
			if((((j-pad)%2==0)&&((j-pad)!=0))||((j-pad)==(f->Ny-1-2*pad))) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
//				else if((k-pad)%2!=0) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2 ;
					k1=addz+(k+pad)/2 + 1;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
			}
//			else if((j-pad)%2!=0) {
//				if((k-pad)%2==0) {
			else if((((j-pad)%2!=0)&&((j-pad)!=(f->Ny-1-2*pad)))||((j-pad)==0)) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2;
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
//				else if((k-pad)%2!=0) {
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
					i0=addx+(i+pad)/2 - 1;
					i1=addx+(i+pad)/2;
					j0=addy+(j+pad)/2;
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2; 
					k1=addz+(k+pad)/2 + 1;
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
//		else if((i-pad)%2!=0) {
//			if((j-pad)%2==0) {
//				if((k-pad)%2==0) {
		else if((((i-pad)%2!=0)&&((i-pad)!=(f->Nx-1-2*pad)))||((i-pad)==0)) {
			if((((j-pad)%2==0)&&((j-pad)!=0))||((j-pad)==(f->Ny-1-2*pad))) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
//				else if((k-pad)%2!=0) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2 - 1;
					j1=addy+(j+pad)/2;
					k0=addz+(k+pad)/2; 
					k1=addz+(k+pad)/2 + 1;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
			}
//			else if((j-pad)%2!=0) {
//				if((k-pad)%2==0) {
			else if((((j-pad)%2!=0)&&((j-pad)!=(f->Ny-1-2*pad)))||((j-pad)==0)) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1-2*pad))) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2; 
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2 - 1;
					k1=addz+(k+pad)/2;
					f->val[i][j][k] = Trilinear_interpolate(meshp->x[i0][j][k],meshp->x[i1][j][k],mesh->x[i][j][k], \
										meshp->y[i][j0][k],meshp->y[i][j1][k],mesh->y[i][j][k], \
										meshp->z[i][j][k0],meshp->z[i][j][k1],mesh->z[i][j][k], \
										fp->val[i0][j0][k1],fp->val[i0][j1][k1], \
										fp->val[i1][j0][k1],fp->val[i1][j1][k1], \
										fp->val[i0][j0][k0],fp->val[i0][j1][k0], \
										fp->val[i1][j0][k0],fp->val[i1][j1][k0]);
				}
//				else if((k-pad)%2!=0) {
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1-2*pad)))||((k-pad)==0)) {
					i0=addx+(i+pad)/2; 
					i1=addx+(i+pad)/2 + 1;
					j0=addy+(j+pad)/2; 
					j1=addy+(j+pad)/2 + 1;
					k0=addz+(k+pad)/2; 
					k1=addz+(k+pad)/2 + 1;
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
}

/*!If restricting from left child, add 0, else if from a right child add -N(domain cells).
 
Parameters: octree node, field index, indices i j k, additions addx addy addz.
 */
double restrict_from_child(Octree* node, Octree* child, int l, int i, int j, int k, int addx, int addy, int addz) {
	//indices
	int i0, i1, j0, j1, k0, k1;	
	double result;	

	//grid and child's grid
	VecField* mesh = node->get_block_data()->mesh;
	VecField* meshc = child->get_block_data()->mesh; 

	//field and child's field
       	Field *f = node->get_block_data()->scalarfields[l];
       	Field *fc = child->get_block_data()->scalarfields[l];

	//domain
	if(node->get_block_data()->flag[i][j][k]==DOMAIN) {
		std::cout << node->get_block_data()->flag[i][j][k] << std::endl;
		i0=addx+2*i-pad;
		i1=addx+2*i-pad+1;
		j0=addy+2*j-pad; 
		j1=addy+2*j-pad+1;
		k0=addz+2*k-pad;
		k1=addz+2*k-pad+1;
		result = Trilinear_interpolate(meshc->x[i0][j][k],meshc->x[i1][j][k],mesh->x[i][j][k], \
							meshc->y[i][j0][k],meshc->y[i][j1][k],mesh->y[i][j][k], \
							meshc->z[i][j][k0],meshc->z[i][j][k1],mesh->z[i][j][k], \
							fc->val[i0][j0][k1],fc->val[i0][j1][k1], \
							fc->val[i1][j0][k1],fc->val[i1][j1][k1], \
							fc->val[i0][j0][k0],fc->val[i0][j1][k0], \
							fc->val[i1][j0][k0],fc->val[i1][j1][k0]);
	}

	return result;
}


double restrict_at(Octree* node, int l, int i, int j, int k) {
			
	//field and block
       	Field *f = node->get_block_data()->scalarfields[l];
       	Block *b = node->get_block_data();

	double result;

	if(i>=f->Nx/2) {
		if(j>=f->Ny/2) {
			if(k>=f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(1,1,1),l,i,j,k,-b->iNx,-b->iNy,-b->iNz);
			}
			else if(k<f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(1,1,0),l,i,j,k,-b->iNx,-b->iNy,0);
			}
		}
		else if(j<f->Ny/2) {
			if(k>=f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(1,0,1),l,i,j,k,-b->iNx,0,-b->iNz);
			}
			else if(k<f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(1,0,0),l,i,j,k,-b->iNx,0,0);
			}
		}
	}
	if(i<f->Nx/2) {
		if(j>=f->Ny/2) {
			if(k>=f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(0,1,1),l,i,j,k,0,-b->iNy,-b->iNz);
			}
			else if(k<f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(0,1,0),l,i,j,k,0,-b->iNy,0);
			}
		}
		else if(j<f->Ny/2) {
			if(k>=f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(0,0,1),l,i,j,k,0,0,-b->iNz);
			}
			else if(k<f->Nz/2) {
				result = restrict_from_child(node,node->get_child_at(0,0,0),l,i,j,k,0,0,0);
			}
		}
	}

	return result;
}



/*!Prolongation of the given field's ghost cells from level-1 to level.*/
void prolongate_ghost(int level, std::string name) {

	//relative locations wrt siblings
	int *loc_x, *loc_y, *loc_z;
	loc_x = new int;
	loc_y = new int;
	loc_z = new int;

	myOctree::create_lists_of_level_nodes();
	
	if(myOctree::level_nodes[level].empty()) {
                std::cerr << "Error! No blocks in level " << level << std::endl;
                exit(1);
        }
	
	if(level==0) {
                std::cerr << "Error! Level cannot be zero" << std::endl;
                exit(1);
        }

	//loop over nodes
	for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {

        	//get location relative to siblings
            	(*it)->get_relative_location(loc_x,loc_y,loc_z);

            	//loop over scalar fields
            	for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {


			//field and block
       			Field *f = (*it)->get_block_data()->scalarfields[l];
       			Block *b = (*it)->get_block_data();

               		if( f->name == name ) {

  				for(int i=0; i<f->Nx; i++) {
         				for(int j=0; j<f->Ny; j++) {
                         			for(int k=0; k<f->Nz; k++) {
								

							if(*loc_x==1) {
								if(*loc_y==1) {
									if(*loc_z==1) {
										prolongate_ghost_for_child((*it),l,i,j,k,b->iNx/2,b->iNy/2,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_ghost_for_child((*it),l,i,j,k,b->iNx/2,b->iNy/2,0);
									}
								}
								else if(*loc_y==0) {
									if(*loc_z==1) {
										prolongate_ghost_for_child((*it),l,i,j,k,b->iNx/2,0,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_ghost_for_child((*it),l,i,j,k,b->iNx/2,0,0);
									}
								}
							}
							if(*loc_x==0) {
								if(*loc_y==1) {
									if(*loc_z==1) {
										prolongate_ghost_for_child((*it),l,i,j,k,0,b->iNy/2,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_ghost_for_child((*it),l,i,j,k,0,b->iNy/2,0);
									}
								}
								else if(*loc_y==0) {
									if(*loc_z==1) {
										prolongate_ghost_for_child((*it),l,i,j,k,0,0,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_ghost_for_child((*it),l,i,j,k,0,0,0);
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


/*!Prolongation of the given field's domain cells from level-1 to level.
 
Parameters: level, field name.
  */
void prolongate_domain(int level, std::string name) {

	//relative locations wrt siblings
	int *loc_x, *loc_y, *loc_z;
	loc_x = new int;
	loc_y = new int;
	loc_z = new int;

	myOctree::create_lists_of_level_nodes();
	
	if(myOctree::level_nodes[level].empty()) {
                std::cerr << "Error! No blocks in level " << level << std::endl;
                exit(1);
        }
	
	if(level==0) {
                std::cerr << "Error! Level cannot be zero" << std::endl;
                exit(1);
        }

	//loop over nodes
	for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {

        	//get location relative to siblings
            	(*it)->get_relative_location(loc_x,loc_y,loc_z);

            	//loop over scalar fields
            	for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {


			//field and block
       			Field *f = (*it)->get_block_data()->scalarfields[l];
       			Block *b = (*it)->get_block_data();

               		if( f->name == name ) {

  				for(int i=0; i<f->Nx; i++) {
         				for(int j=0; j<f->Ny; j++) {
                         			for(int k=0; k<f->Nz; k++) {
								
							if(*loc_x==1) {
								if(*loc_y==1) {
									if(*loc_z==1) {
										prolongate_domain_for_child((*it),l,i,j,k,b->iNx/2,b->iNy/2,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_domain_for_child((*it),l,i,j,k,b->iNx/2,b->iNy/2,0);
									}
								}
								else if(*loc_y==0) {
									if(*loc_z==1) {
										prolongate_domain_for_child((*it),l,i,j,k,b->iNx/2,0,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_domain_for_child((*it),l,i,j,k,b->iNx/2,0,0);
									}
								}
							}
							if(*loc_x==0) {
								if(*loc_y==1) {
									if(*loc_z==1) {
										prolongate_domain_for_child((*it),l,i,j,k,0,b->iNy/2,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_domain_for_child((*it),l,i,j,k,0,b->iNy/2,0);
									}
								}
								else if(*loc_y==0) {
									if(*loc_z==1) {
										prolongate_domain_for_child((*it),l,i,j,k,0,0,b->iNz/2);
									}
									else if(*loc_z==0) {
										prolongate_domain_for_child((*it),l,i,j,k,0,0,0);
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

/*!Restriction of the given field's domain cells from level+1 to level.
 
Parameters: level, field name.
  */
void restrict(int level, std::string name) {

	myOctree::create_lists_of_level_nodes();
	
	if(myOctree::level_nodes[level].empty()) {
                std::cerr << "Error! No blocks in level " << level << std::endl;
                exit(1);
        }
	

	//loop over nodes
	for (std::list<Octree*>::iterator it = level_nodes[level].begin(), end = level_nodes[level].end(); it != end; ++it) {


		//skip if its a leaf node
		if((*it)->isLeafNode()) {	
			continue;
		}
		
            	//loop over scalar fields
            	for(int l = 0; l<myOctree::scalar_fields.size() ; l++) {

			//field and block
       			Field *f = (*it)->get_block_data()->scalarfields[l];
               		if( f->name == name ) {

  				for(int i=0; i<f->Nx; i++) {
         				for(int j=0; j<f->Ny; j++) {
                         			for(int k=0; k<f->Nz; k++) {
			
					//	std::cout << i << " " << j << " " << k << std::endl;
				
							f->val[i][j][k] = restrict_at((*it),l,i,j,k);

						}
          				}
                 		}
       			}
         	}
        }
}



void multigrid(std::string name) {

     	gauss_seidel(0, name);


	for(int i=1; i<=myOctree::max_level; i++) {
        	prolongate_ghost(i,name);
       		prolongate_domain(i,name);
     		gauss_seidel(i, name);
	}

	for(int i=myOctree::max_level-1; i>=0; i--) {
		restrict(i,name);
		gauss_seidel(i,name);		
	}
}


}
