#include "octreegrid.h"
#include "vtk.h"
#include <complex>
#include <iostream>

namespace myOctree {

int  max_level;


//sets refine criterion for all the leaf nodes
void set_refinement_criteria() {

        create_list_of_leaf_nodes();

        for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {

                if((*i)->get_level() >= max_level)
                        continue;

                if((*i)->contains(0.1,0.9,0.99))
                        (*i)->set_to_refine_with_nesting();
                if((*i)->contains(1.1,1.1,0.99))
                        (*i)->set_to_refine_with_nesting();
        //      if((*i)->contains(2.1,0.9,0.99))
        //              (*i)->set_to_refine_with_nesting();
        //      if((*i)->contains(5.9,0.1,0.99))
        //              (*i)->set_to_refine_with_nesting();
        //      if((*i)->contains(12.123,1.789,0.99))
        //              (*i)->set_to_refine_with_nesting();
        //      if((*i)->contains(18.9,2.9,0.99))
        //              (*i)->set_to_refine_with_nesting();

        }
}


//refines the leaf nodes based on the criteria
void refine_nodes() {

        create_list_of_leaf_nodes();

        for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {

                if((*i)->get_level() >= max_level)
                        continue;

                if((*i)->setToRefine)
                        (*i)->refine();

        }

}


//sets coarse criterion for all the leaf nodes
void set_coarse_criteria() {

        create_list_of_leaf_nodes();

        for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {

                if(!((*i)->isRootNode())) {

                        //setting coarsening criteria to siblings       
                        if((*i)->contains(1.1,0.9,0.99)) 
                                (*i)->set_to_coarsen_with_nesting();
                                
                        if((*i)->contains(0.9,1.1,0.99))
                                (*i)->set_to_coarsen_with_nesting();
                                
                        if((*i)->contains(1.1,1.1,0.99))
                                (*i)->set_to_coarsen_with_nesting();
                                
                        if((*i)->contains(0.9,0.9,0.99))
                                (*i)->set_to_coarsen_with_nesting();
                                
                }               

        }
}


//refines the leaf nodes based on the criteria
void coarsen_nodes() {

        create_list_of_leaf_nodes();

        for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end;) {

                if((*i)->setToCoarsen && !((*i)->isRootNode()) && ((*i)->isLeafNode())) {

                        //printf("deleting\n");
                        //setting the children of its parent to NULL
                        for(int n=0; n<2; n++) {
                                for(int m=0; m<2; m++) {
                                        for(int l=0; l<2; l++) {
                                                (*i)->get_parent()->set_child_null_at(l, m, n);
                                        }       
                                }       
                        }       

                        i = nodes.erase(i);
                }
                else {
                        ++i;
                }       
        }

}


//resets refine flags 
//Hence to be called after refinement
void reset_refine_flags() {

        for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {

                (*i)->setToRefine = false;
        }

}

//resets coarsen flags 
//Hence to be called after coarsening
void reset_coarsen_flags() {

        for (std::list<Octree*>::iterator i = nodes.begin(), end = nodes.end(); i != end; ++i) {

                (*i)->setToCoarsen = false;
        }
}


//adapts grid based on the field gradient
void set_refine_flag_based_on_gradient() {

        double x_grad, y_grad, z_grad;
        double dx, dy, dz;
        double mg;
        double tg;
        double max_of_max_gradient = 0.0;


        create_list_of_leaf_nodes();


        //calculation of max_total_gradient
        for (std::list<Octree*>::iterator it = leaf_nodes.begin(), end = leaf_nodes.end(); it != end; ++it) {

                //field based on which gradient is to calculated
                Field* f = (*it)->get_block_data()->field;
                dx = (*it)->get_block_data()->dx;
                dy = (*it)->get_block_data()->dy;
                dz = (*it)->get_block_data()->dz;

                mg = 0.0;

                for(int i=pad;i<(f->Nx-pad);i++) {
                        for(int j=pad;j<(f->Ny-pad);j++) {
                                for(int k=pad;k<(f->Nz-pad);k++) {

                                        //central difference gradient calculation 
                                        x_grad = (f->val[i+1][j][k] - f->val[i-1][j][k]) / (2.0*dx);
                                        y_grad = (f->val[i][j+1][k] - f->val[i][j-1][k]) / (2.0*dy);
                                        z_grad = (f->val[i][j][k+1] - f->val[i][j][k-1]) / (2.0*dz);

                                        tg = std::abs(x_grad) + std::abs(y_grad) + std::abs(z_grad);
                                        if(tg>mg)       mg = tg;
                                        //if((*it)->get_block_data()->total_gradient!=0.0)      printf("Yes not zero\n");                              
                                }
                        }
                }

                (*it)->get_block_data()->max_gradient = mg;
                if(mg>max_of_max_gradient)      max_of_max_gradient = mg;
        }

        //setting flags
        for (std::list<Octree*>::iterator it = leaf_nodes.begin(), end = leaf_nodes.end(); it != end; ++it) {

                if((*it)->get_block_data()->max_gradient > max_of_max_gradient*0.2) {
                        (*it)->set_to_refine_with_nesting();
			//printf("yes refine\n");
		}	
        }
}


//adapts grid based on the field gradient
void set_coarsen_flag_based_on_gradient() {

        double x_grad, y_grad, z_grad;
        double dx, dy, dz;
        double mg;
        double tg;
        double max_of_max_gradient = 0.0;


        create_list_of_leaf_nodes();

        //calculation of max_total_gradient
        for (std::list<Octree*>::iterator it = leaf_nodes.begin(), end = leaf_nodes.end(); it != end; ++it) {

                //field based on which gradient is to calculated
                Field* f = (*it)->get_block_data()->field;
                dx = (*it)->get_block_data()->dx;
                dy = (*it)->get_block_data()->dy;
                dz = (*it)->get_block_data()->dz;

                mg = 0.0;

                for(int i=pad;i<(f->Nx-pad);i++) {
                        for(int j=pad;j<(f->Ny-pad);j++) {
                                for(int k=pad;k<(f->Nz-pad);k++) {

                                        //central difference gradient calculation 
                                        x_grad = (f->val[i+1][j][k] - f->val[i-1][j][k]) / (2.0*dx);
                                        y_grad = (f->val[i][j+1][k] - f->val[i][j-1][k]) / (2.0*dy);
                                        z_grad = (f->val[i][j][k+1] - f->val[i][j][k-1]) / (2.0*dz);

                                        tg = std::abs(x_grad) + std::abs(y_grad) + std::abs(z_grad);
                                        if(tg>mg)       mg = tg;
                                        //if((*it)->get_block_data()->total_gradient!=0.0)      printf("Yes not zero\n");                              
                                }
                        }
                }

                (*it)->get_block_data()->max_gradient = mg;
                if(mg>max_of_max_gradient)      max_of_max_gradient = mg;
        }

        //setting flags
        for (std::list<Octree*>::iterator it = leaf_nodes.begin(), end = leaf_nodes.end(); it != end; ++it) {

		//printf("mg mmg = %lf %lf\n",(*it)->get_block_data()->max_gradient,max_of_max_gradient);

                if((*it)->get_block_data()->max_gradient < max_of_max_gradient*0.2) {
                        (*it)->set_to_coarsen_with_nesting();
			//printf("yes coarsen\n");
		}	
        }

}

//checks if siblings of a node are also set to coarsening. Else the flag is removed.
void recheck_siblings_coarsen_flags() {

        create_list_of_leaf_nodes();

        for (std::list<Octree*>::iterator it = leaf_nodes.begin(), end = leaf_nodes.end(); it != end; ++it) {
	
		if((*it)->setToCoarsen == true)	{
			
                        for(int n=0; n<2; n++) {
                                for(int m=0; m<2; m++) {
                                        for(int l=0; l<2; l++) {
						if((*it)->get_parent()->get_child_at(l,m,n)->setToCoarsen == false)
							(*it)->setToCoarsen = false;
					}
				}
			}
		}
	}		

}	

}
