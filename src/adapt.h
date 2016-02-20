#ifndef MYOCTREE_ADAPT_H
#define MYOCTREE_ADAPT_H

#include "octree.h"

namespace myOctree {

void set_refine_criteria();
void refine_nodes();
void set_refine_flag_based_on_gradient();
void reset_refine_flags();

}

#endif
