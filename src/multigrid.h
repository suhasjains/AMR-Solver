#ifndef MULTIGRID_H
#define MULTIGRID_H


namespace amrsolver {

void prolongate_ghost(int, std::string);
void prolongate_domain(int, std::string);
void restrict(int, std::string);
double Trilinear_interpolate(double, double, double, double, double, double, double, double, double, \
                                 double, double, double, double, double, double, double, double);
void prolongate_ghost_for_child(myOctree::Octree*, int, int, int, int, int, int, int);
void prolongate_domain_for_child(myOctree::Octree*, int, int, int, int, int, int, int);
double restrict_from_child(myOctree::Octree*, myOctree::Octree*, int, int, int, int, int, int, int);
double restrict_at(myOctree::Octree*, int, int, int, int);
void multigrid(std::string);

}

#endif
