#ifndef MYOCTREE_MYOCTREE_H
#define MYOCTREE_MYOCTREE_H

#include "octree.h"

namespace myOctree {

void create_list_of_leaf_nodes();
void create_list_of_root_nodes();
void create_lists_of_level_nodes();
void create_node(int, double, double, double, double, double, double, int, NodeBc **bc);
void OctreeGrid();
void reassign_neighbours();

extern std::list<Octree*> nodes;
extern std::list<Octree*> leaf_nodes;
extern std::list<Octree*> root_nodes;
extern std::list<Octree*> level_nodes[];


}
#endif
