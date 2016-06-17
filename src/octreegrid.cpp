#include "octree.h"
#include "vtk.h"
#include "boundary.h"
#include <iostream>
#include "direction.h"
#include "ghost.h"

///AMR grid stuff
/*!Namespace containing things related to grid*/
namespace myOctree {

extern int max_level;

std::list<Octree*> nodes;
std::list<Octree*> leaf_nodes;
std::list<Octree*> root_nodes;
std::list<Octree*> level_nodes[20];


/*!Creates a list of leaf nodes. Clears the old list before creating the new one.*/ 
void create_list_of_leaf_nodes() {

	leaf_nodes.clear();

	for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {
		if((*i)->isLeafNode()) 
			leaf_nodes.push_back(*i);
	}
}

/*!Creates a new list of root nodes. Clears the old list before creating the new one.*/
void create_list_of_root_nodes() {

	root_nodes.clear();

	for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {
		if((*i)->isRootNode()) 
			root_nodes.push_back(*i);
	}
}

/*!Creates a vector of lists of nodes on each level. Clears the old lists before creating the new one.*/
void create_lists_of_level_nodes() {

	//clearing all lists
	for (unsigned level=0; level <= max_level; level++) {
		level_nodes[level].clear();
	}
	

	//pushing nodes to respective lists	
	for (std::list<Octree*>::iterator j = nodes.begin(), end = nodes.end(); j != end; ++j) {
                int level = (*j)->get_level();
		level_nodes[level].push_back(*j);
	}

}

/*!Creates an octree node*/
void create_node(int blocknumber, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int level, NodeBc **bc) {

	//memory allocation to new node
	Octree r(xmin,xmax,ymin,ymax,zmin,zmax,level);
	nodes.pop_back();
	Octree* root = new Octree(r);

	//boundaries are assigned to this root node
//	root->east_bc = east_bc;
//	root->west_bc = west_bc;
//	root->north_bc = north_bc;
//	root->south_bc = south_bc;
//	root->top_bc = top_bc;
//	root->bottom_bc = bottom_bc;

	for(int i=0;i<3;i++) {
		for(int j=0;j<2;j++) {
			root->bc[i][j] = bc[i][j];
		}
	}

	root->number = blocknumber;
		
}

/*!Reassigns neighbours after coarsening or refining the grid and after creating lists of level_nodes. 
  
Neighbours are assigned during refinement but since coarsening and refinement is a sequential process, neighbours might not have been created and a NULL would have been assigned.

So this function is to be called after all refinement or coarsening is complete in every time step to ensure correct assignment of neighbours.*/ 
void reassign_neighbours() {

	for (unsigned level=0; level <= max_level; level++) {


		for (std::list<Octree*>::iterator iter = level_nodes[level].begin(), end = level_nodes[level].end(); iter != end; ++iter) { 
			//if it is not leaf node	
			if(!((*iter)->isLeafNode())) {
	
				//setting neighbours to children
			        for(int k=0; k<2; k++) {
			                for(int j=0; j<2; j++) {
			                        for(int i=0; i<2; i++) {
			
			                               // if(i==0) {
			                               //         if((*iter)->west != NULL)  (*iter)->get_child_at(i, j, k)->west = (*iter)->west->get_child_at(i+1, j, k);
			                               // }
			                               // if(j==0) {
			                               //         if((*iter)->south != NULL) (*iter)->get_child_at(i, j, k)->south = (*iter)->south->get_child_at(i, j+1, k);
			                               // }
			                               // if(k==0) {
			                               //         if((*iter)->bottom != NULL)	(*iter)->get_child_at(i, j, k)->bottom = (*iter)->bottom->get_child_at(i, j, k+1);
			                               // }
			                               // if(i==1) {
			                               //         if((*iter)->east != NULL)  (*iter)->get_child_at(i, j, k)->east = (*iter)->east->get_child_at(i-1, j, k);
			                               // }
			                               // if(j==1) {
			                               //         if((*iter)->north != NULL) (*iter)->get_child_at(i, j, k)->north = (*iter)->north->get_child_at(i, j-1, k);
			                               // }
			                               // if(k==1) {
			                               //         if((*iter)->top != NULL)   (*iter)->get_child_at(i, j, k)->top = (*iter)->top->get_child_at(i, j, k-1);
			                               // }
			                                if(i==0) {
			                                        if((*iter)->neighbour[XDIR][LEFT] != NULL)  (*iter)->get_child_at(i, j, k)->neighbour[XDIR][LEFT] = (*iter)->neighbour[XDIR][LEFT]->get_child_at(i+1, j, k);
			                                }
			                                if(j==0) {
			                                        if((*iter)->neighbour[YDIR][LEFT] != NULL) (*iter)->get_child_at(i, j, k)->neighbour[YDIR][LEFT] = (*iter)->neighbour[YDIR][LEFT]->get_child_at(i, j+1, k);
			                                }
			                                if(k==0) {
			                                        if((*iter)->neighbour[ZDIR][LEFT] != NULL)	(*iter)->get_child_at(i, j, k)->neighbour[ZDIR][LEFT] = (*iter)->neighbour[ZDIR][LEFT]->get_child_at(i, j, k+1);
			                                }
			                                if(i==1) {
			                                        if((*iter)->neighbour[XDIR][RIGHT] != NULL)  (*iter)->get_child_at(i, j, k)->neighbour[XDIR][RIGHT] = (*iter)->neighbour[XDIR][RIGHT]->get_child_at(i-1, j, k);
			                                }
			                                if(j==1) {
			                                        if((*iter)->neighbour[YDIR][RIGHT] != NULL) (*iter)->get_child_at(i, j, k)->neighbour[YDIR][RIGHT] = (*iter)->neighbour[YDIR][RIGHT]->get_child_at(i, j-1, k);
			                                }
			                                if(k==1) {
			                                        if((*iter)->neighbour[ZDIR][RIGHT] != NULL)   (*iter)->get_child_at(i, j, k)->neighbour[ZDIR][RIGHT] = (*iter)->neighbour[ZDIR][RIGHT]->get_child_at(i, j, k-1);
			                                }
			                        }
			                }
			        }
			}
		}
	}
}

/*!Assigns neighours to root nodes*/
void set_root_neighbours() {

	create_list_of_root_nodes();

	for (std::list<Octree*>::iterator i = root_nodes.begin(), end = root_nodes.end(); i != end; ++i) {

		double xmax = (*i)->x_max;			
		double ymax = (*i)->y_max;			
		double zmax = (*i)->z_max;			
		double xmin = (*i)->x_min;			
		double ymin = (*i)->y_min;			
		double zmin = (*i)->z_min;

		for (std::list<Octree*>::iterator j = root_nodes.begin(), end = root_nodes.end(); j != end; ++j) {
			/*change it to xmax_nbr*/
			double x_max = (*j)->x_max;			
			double y_max = (*j)->y_max;			
			double z_max = (*j)->z_max;			
			double x_min = (*j)->x_min;			
			double y_min = (*j)->y_min;			
			double z_min = (*j)->z_min;
	
		//	if(xmax==x_min&&ymax==y_max&&zmax==z_max)
		//		(*i)->east = (*j);	

		//	if(xmin==x_max&&ymax==y_max&&zmax==z_max)
		//		(*i)->west = (*j);	
		//	
		//	if(xmax==x_max&&ymax==y_min&&zmax==z_max)
		//		(*i)->north = (*j);	
		//		
		//	if(xmax==x_max&&ymin==y_max&&zmax==z_max)
		//		(*i)->south = (*j);	
		//	
		//	if(xmax==x_max&&ymax==y_max&&zmax==z_min)
		//		(*i)->top = (*j);	
		//	
		//	if(xmax==x_max&&ymax==y_max&&zmin==z_max)
		//		(*i)->bottom = (*j);	
			if(xmax==x_min&&ymax==y_max&&zmax==z_max)
				(*i)->neighbour[XDIR][RIGHT] = (*j);	

			if(xmin==x_max&&ymax==y_max&&zmax==z_max)
				(*i)->neighbour[XDIR][LEFT] = (*j);	
			
			if(xmax==x_max&&ymax==y_min&&zmax==z_max)
				(*i)->neighbour[YDIR][RIGHT] = (*j);	
				
			if(xmax==x_max&&ymin==y_max&&zmax==z_max)
				(*i)->neighbour[YDIR][LEFT] = (*j);	
			
			if(xmax==x_max&&ymax==y_max&&zmax==z_min)
				(*i)->neighbour[ZDIR][RIGHT] = (*j);	
			
			if(xmax==x_max&&ymax==y_max&&zmin==z_max)
				(*i)->neighbour[ZDIR][LEFT] = (*j);	
				
				
		}			
	}
}	

/*!Prints the centre coordinates of the neighbours of the given list of nodes. 
 
  This can be used to test correctness in the assignment of neighbours.*/
 
void print_neighbour_information(std::list<Octree*>& nodes) {

	for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {
	
		printf("I am at %g %g %g\n",(*i)->x_centre,(*i)->y_centre,(*i)->z_centre);
	//	if((*i)->east!=NULL)	printf("East neighbour at %g %g %g\n",(*i)->east->x_centre,(*i)->east->y_centre,(*i)->east->z_centre);		
	//	if((*i)->west!=NULL)	printf("West neighbour at %g %g %g\n",(*i)->west->x_centre,(*i)->west->y_centre,(*i)->west->z_centre);		
	//	if((*i)->north!=NULL)	printf("North neighbour at %g %g %g\n",(*i)->north->x_centre,(*i)->north->y_centre,(*i)->north->z_centre);		
	//	if((*i)->south!=NULL)	printf("South neighbour at %g %g %g\n",(*i)->south->x_centre,(*i)->south->y_centre,(*i)->south->z_centre);		
	//	if((*i)->top!=NULL)	printf("Top neighbour at %g %g %g\n",(*i)->top->x_centre,(*i)->top->y_centre,(*i)->top->z_centre);		
	//	if((*i)->bottom!=NULL)	printf("Bottom neighbour at %g %g %g\n",(*i)->bottom->x_centre,(*i)->bottom->y_centre,(*i)->bottom->z_centre);		
		if((*i)->neighbour[XDIR][RIGHT]!=NULL)	printf("East neighbour at %g %g %g\n",(*i)->neighbour[XDIR][RIGHT]->x_centre,(*i)->neighbour[XDIR][RIGHT]->y_centre,(*i)->neighbour[XDIR][RIGHT]->z_centre);		
		if((*i)->neighbour[XDIR][LEFT]!=NULL)	printf("West neighbour at %g %g %g\n",(*i)->neighbour[XDIR][LEFT]->x_centre,(*i)->neighbour[XDIR][LEFT]->y_centre,(*i)->neighbour[XDIR][LEFT]->z_centre);		
		if((*i)->neighbour[YDIR][RIGHT]!=NULL)	printf("North neighbour at %g %g %g\n",(*i)->neighbour[YDIR][RIGHT]->x_centre,(*i)->neighbour[YDIR][RIGHT]->y_centre,(*i)->neighbour[YDIR][RIGHT]->z_centre);		
		if((*i)->neighbour[YDIR][LEFT]!=NULL)	printf("South neighbour at %g %g %g\n",(*i)->neighbour[YDIR][LEFT]->x_centre,(*i)->neighbour[YDIR][LEFT]->y_centre,(*i)->neighbour[YDIR][LEFT]->z_centre);		
		if((*i)->neighbour[ZDIR][RIGHT]!=NULL)	printf("Top neighbour at %g %g %g\n",(*i)->neighbour[ZDIR][RIGHT]->x_centre,(*i)->neighbour[ZDIR][RIGHT]->y_centre,(*i)->neighbour[ZDIR][RIGHT]->z_centre);		
		if((*i)->neighbour[ZDIR][LEFT]!=NULL)	printf("Bottom neighbour at %g %g %g\n",(*i)->neighbour[ZDIR][LEFT]->x_centre,(*i)->neighbour[ZDIR][LEFT]->y_centre,(*i)->neighbour[ZDIR][LEFT]->z_centre);		

	}

}

/*!Exchanges ghost values of all the user defined Scalar and Vector fields at the given level*/
void exchange_ghost_values_of_level(int level) {

	for(int l = 0; l<scalar_fields.size() ; l++) {
		amrsolver::exchange_ghost_val(level, scalar_fields[l]);
	}
	
	for(int l = 0; l<vector_fields.size() ; l++) {
		amrsolver::exchange_ghost_val(level, vector_fields[l]);
	}

}

/*!Sets up octree AMR grid*/
void OctreeGrid() {

	std::cerr << "\n"  <<"Setting up grid" << std::endl;

	set_root_neighbours();

	//This is to be done after neighbours are set at that particular level
	exchange_ghost_values_of_level(0);

	//prints neighbours information
	//print_neighbour_information(nodes);

}

}
