#include "octree.h"
#include "output.h"

namespace myOctree {

std::list<Octree*> nodes;
std::list<Octree*> leaf_nodes;
std::list<Octree*> root_nodes;
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

void create_node(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int level) {

	Octree* root = new Octree;
	Octree r(xmin,xmax,ymin,ymax,zmin,zmax,level);
	nodes.pop_back();
	*root = r;

}	

void OctreeGrid() {


	create_node(0.0,1.0,0.0,1.0,0.0,1.0,0);
	create_node(1.0,2.0,0.0,1.0,0.0,1.0,0);
	create_node(1.0,2.0,1.0,2.0,0.0,1.0,0);
	create_node(0.0,1.0,1.0,2.0,0.0,1.0,0);

	for(int i=0;i<=10;i++) {

		set_refine_criteria();
		refine_nodes();
	}

}

}
