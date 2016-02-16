#include "octree.h"
#include "output.h"

namespace myOctree {


#define MAX_LEVEL 4

std::list<Octree*> nodes;
std::list<Octree*> leaf_nodes;
std::list<Octree*> root_nodes;
std::list<Octree*> level_nodes[MAX_LEVEL + 1];
int Block::iNx = NX_BLOCK;
int Block::iNy = NY_BLOCK;
int Block::iNz = NZ_BLOCK;


//creates a new list of leaf nodes 
void create_list_of_leaf_nodes() {

	leaf_nodes.clear();

	for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {
		if((*i)->isLeafNode()) 
			leaf_nodes.push_back(*i);
	}
}

//creates a new list of root nodes
void create_list_of_root_nodes() {

	root_nodes.clear();

	for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {
		if((*i)->isRootNode()) 
			root_nodes.push_back(*i);
	}
}

//creates a vector of lists of nodes on each level
void create_lists_of_level_nodes() {

	//clearing all lists
	for (unsigned level=0; level <= MAX_LEVEL; level++) {
	
		level_nodes[level].clear();
	}

	//pushing nodes to respective lists	
	for (std::list<Octree*>::iterator j = nodes.begin(), end = nodes.end(); j != end; ++j) {
                int level = (*j)->get_level();
		level_nodes[level].push_back(*j);
	}

}

//sets refine criterion for all the leaf nodes
void set_refinement_criteria() {

	create_list_of_leaf_nodes();

	for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {
		
		if((*i)->get_level() >= MAX_LEVEL)
			continue;

        	if((*i)->contains(1.1,0.9,0.99))
 			(*i)->setToRefine = true;
		if((*i)->contains(0.9,1.1,0.99))
 			(*i)->setToRefine = true;
        	if((*i)->contains(1.1,1.1,0.99))
 			(*i)->setToRefine = true;
        	if((*i)->contains(0.9,0.9,0.99))
 			(*i)->setToRefine = true;
	
	}
}

//refines the leaf nodes based on the criteria
void refine_nodes() {

	create_list_of_leaf_nodes();

	for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {
       
		if((*i)->get_level() >= MAX_LEVEL)
			continue;
 
        	if((*i)->setToRefine)
            		(*i)->refine();
    
	}

}

//sets coarse criterion for all the leaf nodes
void set_coarse_criteria() {

	create_list_of_leaf_nodes();

	for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {
	
		if(!((*i)->isRootNode())) {

		//setting coarsening criteria to siblings	
        		if((*i)->contains(1.1,0.9,0.99)) {
				//setting criteria to siblings
				for(int n=0; n<2; n++) {
					for(int m=0; m<2; m++) {
				        	for(int l=0; l<2; l++) {
							(*i)->get_parent()->get_child_at(l, m, n)->setToCoarsen = true;
						}
					}
				}
			}

			if((*i)->contains(0.9,1.1,0.99)) {
				//setting criteria to siblings
				for(int n=0; n<2; n++) {
					for(int m=0; m<2; m++) {
				        	for(int l=0; l<2; l++) {
							(*i)->get_parent()->get_child_at(l, m, n)->setToCoarsen = true;
						}
					}
				}
			}

        		if((*i)->contains(1.1,1.1,0.99)) {
				//setting criteria to siblings
				for(int n=0; n<2; n++) {
					for(int m=0; m<2; m++) {
				        	for(int l=0; l<2; l++) {
							(*i)->get_parent()->get_child_at(l, m, n)->setToCoarsen = true;
						}
					}
				}
			}

        		if((*i)->contains(0.9,0.9,0.99)) {
				//setting criteria to siblings
				for(int n=0; n<2; n++) {
					for(int m=0; m<2; m++) {
				        	for(int l=0; l<2; l++) {
							(*i)->get_parent()->get_child_at(l, m, n)->setToCoarsen = true;
						}
					}
				}
			}

		}
	
	}
}

//refines the leaf nodes based on the criteria
void coarsen_nodes() {

	create_list_of_leaf_nodes();

	for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end;) {
       
        	if((*i)->setToCoarsen && !((*i)->isRootNode()) && ((*i)->isLeafNode())) {
	
			//printf("deleting\n");
			//setting the children of its parent to NULL
			for(int n=0; n<2; n++) {
				for(int m=0; m<2; m++) {
			        	for(int l=0; l<2; l++) {
						(*i)->get_parent()->set_child_null_at(l, m, n);
					}
				}
			}

            		i = nodes.erase(i);	
		}
		else { 
			++i;
		}   
	}
}

void set_field() {

        create_list_of_leaf_nodes();

        for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {

                Field* mesh = (*i)->get_block_data()->mesh;
                int level = (*i)->get_level();

                printf("level = %d\n",level);
                mesh->set_field((double)level);

        }

}

//creates an octree node
void create_node(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int level) {

	Octree* root = new Octree;
	Octree r(xmin,xmax,ymin,ymax,zmin,zmax,level);
	nodes.pop_back();
	*root = r;

}

//reassign neighbours after coarsening or refining the grid and after creating lists of level_nodes
//this is done coz, while assigning neighbours during refinement, neighbours might not have been created and a NULL would be assigned in these situations
//So this is done to ensure that neighbours are assigned after all refinement or coarsening is done for this time step 
void reassign_neighbours() {

	for (unsigned level=0; level <= MAX_LEVEL; level++) {


		for (std::list<Octree*>::iterator iter = level_nodes[level].begin(), end = level_nodes[level].end(); iter != end; ++iter) { 
			//if it is not leaf node	
			if(!((*iter)->isLeafNode())) {
	
				//setting neighbours to children
			        for(int k=0; k<2; k++) {
			                for(int j=0; j<2; j++) {
			                        for(int i=0; i<2; i++) {
			
			                                if(i==0) {
			                                        if((*iter)->west != NULL)  (*iter)->get_child_at(i, j, k)->west = (*iter)->west->get_child_at(i+1, j, k);
			                                }
			                                if(j==0) {
			                                        if((*iter)->south != NULL) (*iter)->get_child_at(i, j, k)->south = (*iter)->south->get_child_at(i, j+1, k);
			                                }
			                                if(k==0) {
			                                        if((*iter)->bottom != NULL)	(*iter)->get_child_at(i, j, k)->bottom = (*iter)->bottom->get_child_at(i, j, k+1);
			                                }
			                                if(i==1) {
			                                        if((*iter)->east != NULL)  (*iter)->get_child_at(i, j, k)->east = (*iter)->east->get_child_at(i-1, j, k);
			                                }
			                                if(j==1) {
			                                        if((*iter)->north != NULL) (*iter)->get_child_at(i, j, k)->north = (*iter)->north->get_child_at(i, j-1, k);
			                                }
			                                if(k==1) {
			                                        if((*iter)->top != NULL)   (*iter)->get_child_at(i, j, k)->top = (*iter)->top->get_child_at(i, j, k-1);
			                                }
			                        }
			                }
			        }
			}
		}
	}
}

//sets neighours to root nodes
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

			double x_max = (*j)->x_max;			
			double y_max = (*j)->y_max;			
			double z_max = (*j)->z_max;			
			double x_min = (*j)->x_min;			
			double y_min = (*j)->y_min;			
			double z_min = (*j)->z_min;
	
			if(xmax==x_min&&ymax==y_max&&zmax==z_max)
				(*i)->east = (*j);	

			if(xmin==x_max&&ymax==y_max&&zmax==z_max)
				(*i)->west = (*j);	
			
			if(xmax==x_max&&ymax==y_min&&zmax==z_max)
				(*i)->north = (*j);	
				
			if(xmax==x_max&&ymin==y_max&&zmax==z_max)
				(*i)->south = (*j);	
			
			if(xmax==x_max&&ymax==y_max&&zmax==z_min)
				(*i)->top = (*j);	
			
			if(xmax==x_max&&ymax==y_max&&zmin==z_max)
				(*i)->bottom = (*j);	
				
		}			
	}
}	

void print_neighbour_information(std::list<Octree*>& nodes) {

	for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {
	
		printf("I am at %g %g %g\n",(*i)->x_centre,(*i)->y_centre,(*i)->z_centre);
		if((*i)->east!=NULL)	printf("East neighbour at %g %g %g\n",(*i)->east->x_centre,(*i)->east->y_centre,(*i)->east->z_centre);		
		if((*i)->west!=NULL)	printf("West neighbour at %g %g %g\n",(*i)->west->x_centre,(*i)->west->y_centre,(*i)->west->z_centre);		
		if((*i)->north!=NULL)	printf("North neighbour at %g %g %g\n",(*i)->north->x_centre,(*i)->north->y_centre,(*i)->north->z_centre);		
		if((*i)->south!=NULL)	printf("South neighbour at %g %g %g\n",(*i)->south->x_centre,(*i)->south->y_centre,(*i)->south->z_centre);		
		if((*i)->top!=NULL)	printf("Top neighbour at %g %g %g\n",(*i)->top->x_centre,(*i)->top->y_centre,(*i)->top->z_centre);		
		if((*i)->bottom!=NULL)	printf("Bottom neighbour at %g %g %g\n",(*i)->bottom->x_centre,(*i)->bottom->y_centre,(*i)->bottom->z_centre);		

	}

}

void OctreeGrid() {


	create_node(0.0,1.0,0.0,1.0,0.0,1.0,0);
	create_node(1.0,2.0,0.0,1.0,0.0,1.0,0);
	create_node(1.0,2.0,1.0,2.0,0.0,1.0,0);
	create_node(0.0,1.0,1.0,2.0,0.0,1.0,0);

	set_root_neighbours();

	//refinement	
	for(int i=0;i<=MAX_LEVEL;i++) {

		set_refinement_criteria();
		refine_nodes();
	}
	
	//coarsening
	for(int i=0;i<=MAX_LEVEL;i++) {

		//printf("coarsening\n");
		set_coarse_criteria();
		coarsen_nodes();
	}

	//do this after refinement or coarsening
	create_lists_of_level_nodes();
	reassign_neighbours();

	//prints neighbours information
	//print_neighbour_information(nodes);
	

}

}
