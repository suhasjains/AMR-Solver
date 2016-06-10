#ifndef MYOCTREE_ADAPT_H
#define MYOCTREE_ADAPT_H

#include "octree.h"

namespace myOctree {

//void set_refinement_criteria();
void refine_nodes();
void set_refine_flag_based_on_gradient();
void reset_refine_flags();
//void set_coarsening_criteria();
void coarsen_nodes();
void set_coarsen_flag_based_on_gradient();
void reset_coarsen_flags();
void recheck_siblings_coarsen_flags();


extern int max_level;
}

#endif
