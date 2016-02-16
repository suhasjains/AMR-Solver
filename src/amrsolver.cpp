#include "myoctree.h"
#include "output.h"
using myOctree::Octree;
using myOctree::Field;

namespace amrsolver {


void set_field() {

	myOctree::create_list_of_leaf_nodes();

    	for (std::list<Octree*>::iterator i = myOctree::leaf_nodes.begin(), end = myOctree::leaf_nodes.end(); i != end; ++i) {

   		Field* mesh = (*i)->get_block_data()->mesh;
		int level = (*i)->get_level();
	
		//printf("level = %d\n",level);	
		mesh->set_field((double)level);	

    	}

}



}


using namespace std;


int main(int argc, char **argv) {

	myOctree::OctreeGrid();

	amrsolver::set_field();

	myOctree::create_list_of_leaf_nodes();

	myOctree::write_vtk(myOctree::level_nodes[3]);

}
