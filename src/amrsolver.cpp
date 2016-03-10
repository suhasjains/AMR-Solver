#include "octreegrid.h"
#include "adapt.h"
#include "vtk.h"
#include "input.h"
#include "output.h"
#include <iostream>

using myOctree::Octree;
using myOctree::Field;
using myOctree::VecField;

namespace amrsolver {


void set_initial_field() {

	myOctree::create_list_of_leaf_nodes();
		std::cerr << myOctree::leaf_nodes.size()  << std::endl;

    	for (std::list<Octree*>::iterator i = myOctree::leaf_nodes.begin(), end = myOctree::leaf_nodes.end(); i != end; ++i) {

   		Field* field = (*i)->get_block_data()->field;
   		VecField* location = (*i)->get_block_data()->mesh;
		//std::cerr << field->Nz  << std::endl;
	
		for(int i=0;i<field->Nx;i++) {
                	for(int j=0;j<field->Ny;j++) {
                        	for(int k=0;k<field->Nz;k++) {
		std::cout << "Hi here at " << i << j << k << std::endl;
                                	
					double x = location->x[i][j][k];
					double y = location->y[i][j][k];
					double z = location->z[i][j][k];
 					
					//if(!field->val[i][j][k]) std::cout << "Not present at" << i << j << k << std::endl; 	

					((x-1.0)*(x-1.0) + (y-1.0)*(y-1.0) + (z-1.0)*(z-1.0) >= 0.421875)?(field->val[i][j][k] = 1.0):(field->val[i][j][k] = 100.0);		
					//((x-1.0)*(x-1.0) + (y-1.0)*(y-1.0) >= 0.5625)?(field->val[i][j][k] = 1.0):(field->val[i][j][k] = 100.0);		
					//(x*x + y*y >= 1.0)?(field->val[i][j][k] = 1.0):(field->val[i][j][k] = 100.0);		
					//if(x*x + y*y >= 3.0)	{field->val[i][j][k] = 100.0; }		
                        	}
                	}
        	}
    	}
}

void set_field() {

	myOctree::create_list_of_leaf_nodes();

    	for (std::list<Octree*>::iterator i = myOctree::leaf_nodes.begin(), end = myOctree::leaf_nodes.end(); i != end; ++i) {

   		Field* field = (*i)->get_block_data()->scalarfields[0];
   		VecField* location = (*i)->get_block_data()->mesh;
	
		for(int i=0;i<field->Nx;i++) {
                	for(int j=0;j<field->Ny;j++) {
                        	for(int k=0;k<field->Nz;k++) {
                                	
					double x = location->x[i][j][k];
					double y = location->y[i][j][k];
					double z = location->z[i][j][k];
					std::cout << "Present at" << i << j << k << std::endl; 	
 					
					((x-1.0)*(x-1.0) + (y-1.0)*(y-1.0) >= 0.5625)?(field->val[i][j][k] = 1.0):(field->val[i][j][k] = 100.0);		
					//if(x*x + y*y >= 3.0)	{field->val[i][j][k] = 100.0; }		
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

		std::cerr << "refine done" << std::endl;

            	//reassigning neighbours after every level of refine call
            	myOctree::create_lists_of_level_nodes();
            	myOctree::reassign_neighbours();
    		myOctree::reset_refine_flags();	
//    	}

    	
	//coaresning
	/*change number to a variable*/    
  //  	for(int i=0;i<5;i++) {

		std::cerr << "coarsening yet to be done" << std::endl;

		set_initial_field();
		std::cerr << "coarsening yet to be done" << std::endl;

      		myOctree::set_coarsen_flag_based_on_gradient();
		std::cerr << "coarsening yet to be done" << std::endl;

		myOctree::recheck_siblings_coarsen_flags();
		std::cerr << "coarsening yet to be done" << std::endl;

		myOctree::coarsen_nodes();

            	//reassigning neighbours after every level of refine call
            	myOctree::create_lists_of_level_nodes();
		std::cerr << "creating list" << std::endl;
            	myOctree::reassign_neighbours();
    		myOctree::reset_coarsen_flags();	
		std::cerr << "coarsening done" << std::endl;
    	}

	amrsolver::set_initial_field();
		std::cerr << "all done" << std::endl;
	
}

}


using namespace std;


int main(int argc, char **argv) {

	read_input_file();

	myOctree::OctreeGrid();

	cerr << "\n" << "adapting" << endl;

	amrsolver::adapt_gradient();
	cerr << "\n" << "setting field" << endl;
	
	amrsolver::set_field();
	
	myOctree::create_list_of_leaf_nodes();

	myOctree::write_vtk(myOctree::leaf_nodes);
	
	write_output_file();

	cerr << "\n" << "Finishing" << endl;
}
