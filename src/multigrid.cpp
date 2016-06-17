#include <stdlib.h>
#include <string>
#include "octreegrid.h"
#include <iostream>
//#include <cmath>
#include "multigrid.h"

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
		if((((i-pad)%2==0)&&((i-pad)!=0))||((i-pad)==(f->Nx-1))) {
			if((((j-pad)%2==0)&&((j-pad)!=0))||((j-pad)==(f->Ny-1))) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1))) {
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
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1)))||((k-pad)==0)) {
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
			else if((((j-pad)%2!=0)&&((j-pad)!=(f->Ny-1)))||((j-pad)==0)) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1))) {
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
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1)))||((k-pad)==0)) {
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
		else if((((i-pad)%2!=0)&&((i-pad)!=(f->Nx-1)))||((i-pad)==0)) {
			if((((j-pad)%2==0)&&((j-pad)!=0))||((j-pad)==(f->Ny-1))) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1))) {
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
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1)))||((k-pad)==0)) {
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
			else if((((j-pad)%2!=0)&&((j-pad)!=(f->Ny-1)))||((j-pad)==0)) {
				if((((k-pad)%2==0)&&((k-pad)!=0))||((k-pad)==(f->Nz-1))) {
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
				else if((((k-pad)%2!=0)&&((k-pad)!=(f->Nz-1)))||((k-pad)==0)) {
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
		if((i-pad)%2==0) {
			if((j-pad)%2==0) {
				if((k-pad)%2==0) {
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
				else if((k-pad)%2!=0) {
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
			else if((j-pad)%2!=0) {
				if((k-pad)%2==0) {
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
				else if((k-pad)%2!=0) {
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
		else if((i-pad)%2!=0) {
			if((j-pad)%2==0) {
				if((k-pad)%2==0) {
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
				else if((k-pad)%2!=0) {
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
			else if((j-pad)%2!=0) {
				if((k-pad)%2==0) {
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
				else if((k-pad)%2!=0) {
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


/*!Prolongation of the given field's ghost cells from level-1 to level.
 
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

void multigrid(std::string name) {

     	gauss_seidel(0, name);
        prolongate_ghost(1,name);
       	prolongate_domain(1,name);
	

	//not a problem with exchange of ghost, its about prolongation
	//exchange_ghost_val(1,name);
     	//gauss_seidel(1, name);
        prolongate_ghost(2,name);
       	prolongate_domain(2,name);
     	//gauss_seidel(2, name);
        //prolongate_ghost(3,name);
       	//prolongate_domain(3,name);
     	//gauss_seidel(3, name);

}


}
