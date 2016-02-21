#include "octreegrid.h"
#include "adapt.h"
#include "output.h"
using myOctree::Octree;
using myOctree::Field;
using myOctree::VecField;

namespace amrsolver {


void set_initial_field() {

	myOctree::create_list_of_leaf_nodes();

    	for (std::list<Octree*>::iterator i = myOctree::leaf_nodes.begin(), end = myOctree::leaf_nodes.end(); i != end; ++i) {

   		Field* field = (*i)->get_block_data()->field;
   		VecField* location = (*i)->get_block_data()->mesh;
		//int level = (*i)->get_level();
	
		//field->set_field((double)level);	
		for(int i=0;i<field->Nx;i++) {
                	for(int j=0;j<field->Ny;j++) {
                        	for(int k=0;k<field->Nz;k++) {
                                	
					double x = location->x[i][j][k];
					double y = location->y[i][j][k];
					double z = location->z[i][j][k];
 					
					(x*x + y*y >= 1.0)?(field->val[i][j][k] = 1.0):(field->val[i][j][k] = 100.0);		
					if(x*x + y*y >= 3.0)	{field->val[i][j][k] = 100.0; }		
                        	}
                	}
        	}

    	}

}

void adapt_gradient() {

    	//refinement
	/*change number to a variable*/    
    	for(int i=0;i<=myOctree::max_level;i++) {


		set_initial_field();
      		myOctree::set_refine_flag_based_on_gradient();
		myOctree::refine_nodes();

            	//reassigning neighbours after every level of refine call
            	myOctree::create_lists_of_level_nodes();
            	myOctree::reassign_neighbours();
    		myOctree::reset_refine_flags();	
//    	}

    	
	//coaresning
	/*change number to a variable*/    
  //  	for(int i=0;i<5;i++) {


		set_initial_field();
      		myOctree::set_coarsen_flag_based_on_gradient();
		myOctree::recheck_siblings_coarsen_flags();
		myOctree::coarsen_nodes();

            	//reassigning neighbours after every level of refine call
            	myOctree::create_lists_of_level_nodes();
            	myOctree::reassign_neighbours();
    		myOctree::reset_coarsen_flags();	
    	}

	amrsolver::set_initial_field();
	
}

}


using namespace std;


int main(int argc, char **argv) {

	myOctree::OctreeGrid();

	amrsolver::adapt_gradient();
	

	myOctree::create_list_of_leaf_nodes();

	myOctree::write_vtk(myOctree::leaf_nodes);

}
