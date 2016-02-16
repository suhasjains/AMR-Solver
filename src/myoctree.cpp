#include "octree.h"
#include "output.h"

namespace myOctree {

std::list<Octree*> nodes;
std::list<Octree*> leaf_nodes;
std::list<Octree*> root_nodes;
std::vector<std::list<Octree*> > level_nodes;
int Block::iNx = NX_BLOCK;
int Block::iNy = NY_BLOCK;
int Block::iNz = NZ_BLOCK;

#define MAX_LEVEL 4

//creates a new list of leaf nodes 
void create_list_of_leaf_nodes() {

	leaf_nodes.clear();

	for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {
		if((*iterator)->isLeafNode()) 
			leaf_nodes.push_back(*iterator);
	}
}

//creates a new list of root nodes
void create_list_of_root_nodes() {

	root_nodes.clear();

	for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {
		if((*iterator)->isRootNode()) 
			root_nodes.push_back(*iterator);
	}
}

//creates a vector of lists of nodes on each level
void create_lists_of_level_nodes() {

	
	

}

//sets refine criterion for all the leaf nodes
void set_refine_criteria() {

	create_list_of_leaf_nodes();

	for (std::list<Octree*>::iterator iterator = leaf_nodes.begin(), end = leaf_nodes.end(); iterator != end; ++iterator) {
		
		if((*iterator)->get_level() >= MAX_LEVEL)
			continue;

        	if((*iterator)->contains(1.1,0.9,0.99))
 			(*iterator)->setToRefine = true;
		if((*iterator)->contains(0.9,1.1,0.99))
 			(*iterator)->setToRefine = true;
        	if((*iterator)->contains(1.1,1.1,0.99))
 			(*iterator)->setToRefine = true;
        	if((*iterator)->contains(0.9,0.9,0.99))
 			(*iterator)->setToRefine = true;
	
	}
}

//refines all the leaf nodes
void refine_nodes() {

	create_list_of_leaf_nodes();

	for (std::list<Octree*>::iterator iterator = leaf_nodes.begin(), end = leaf_nodes.end(); iterator != end; ++iterator) {
       
		if((*iterator)->get_level() >= MAX_LEVEL)
			continue;
 
        	if((*iterator)->setToRefine)
            		(*iterator)->refine();
    
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

	
	for(int i=0;i<=10;i++) {

		set_refine_criteria();
		refine_nodes();
	}

	//create_list_of_leaf_nodes();
	//print_neighbour_information(nodes);
	

}

}
