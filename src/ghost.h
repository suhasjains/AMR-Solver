#ifndef GHOST_H
#define GHOST_H


namespace amrsolver {


void exchange_ghost_val(int, std::string);
void prolongate_multilevel_ghost_exchange_for_child(myOctree::Octree*, myOctree::Octree*, int, int, int, int, int, int, int, int, int, int);
void prolongate_for_multilevel_ghost_exchange_at(myOctree::Octree*, myOctree::Octree*, int, int, int, int, int, int, int, int, int, int);
void exchange_multilevel_ghost_val(int, std::string);
}

#endif
