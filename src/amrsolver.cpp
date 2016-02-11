#include "myoctree.h"
#include <list>

namespace amrsolver {

void set_field() {

	myOctree::create_list_of_leaf_nodes();

    	for (std::list<myOctree::Octree*>::iterator i = myOctree::leaf_nodes.begin(), end = myOctree::leaf_nodes.end(); i != end; ++i) {

   		myOctree::Field* mesh = (*i)->get_block_data()->mesh;
		int level = (*i)->get_level();
		
		mesh->set_field((double)level);	

    	}

}



}


using namespace amrsolver;


int main(int argc, char **argv) {

	myOctree::OctreeGrid();

	set_field();


}
