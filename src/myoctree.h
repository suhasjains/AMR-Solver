#ifndef MYOCTREE_MYOCTREE_H
#define MYOCTREE_MYOCTREE_H

#include "octree.h"

namespace myOctree {

void create_list_of_leaf_nodes();
void create_list_of_root_nodes();
void set_refine_criteria();
void refine_nodes();
void create_node(double, double, double, double, double, double, int);
void OctreeGrid();

extern std::list<Octree*> nodes;
extern std::list<Octree*> leaf_nodes;
extern std::list<Octree*> root_nodes;


}
#endif
