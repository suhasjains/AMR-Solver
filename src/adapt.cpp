#include "octreegrid.h"
#include "output.h"
#include <complex>


namespace myOctree {

#define MAX_LEVEL 4


//sets refine criterion for all the leaf nodes
void set_refinement_criteria() {

        create_list_of_leaf_nodes();

        for (std::list<Octree*>::iterator i = leaf_nodes.begin(), end = leaf_nodes.end(); i != end; ++i) {

                if((*i)->get_level() >= MAX_LEVEL)
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

                if((*i)->get_level() >= MAX_LEVEL)
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

                for(int i=PAD;i<(f->Nx-PAD);i++) {
                        for(int j=PAD;j<(f->Ny-PAD);j++) {
                                for(int k=PAD;k<(f->Nz-PAD);k++) {

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

                //printf("max grad = %lf grad = %lf \n",max_total_gradient,(*it)->get_block_data()->total_gradient);            

                if((*it)->get_block_data()->max_gradient > max_of_max_gradient*0.2)
                //if((*it)->get_block_data()->total_gradient >= max_total_gradient/2.0)
                        (*it)->setToRefine = true;
        }

}


}
