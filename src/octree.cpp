#include "octree.h"
#include <iostream>

namespace myOctree {



//Member functions
bool Octree::isLeafNode() {
    
    return children[0][0][0] == NULL;
}

bool Octree::isRootNode() {
    
    return parent == NULL;
}

Block* Octree::get_block_data() {
    
    return block_data;
}

bool Octree::contains(double x, double y, double z) {

	/*change it to (x-x_min)*(x-x_max)<=....*/
	return x<=x_max && x>=x_min && y<=y_max && y>=y_min && z<=z_max && z>=z_min;

}

int Octree::get_level() {

	return this->level;

}

Octree* Octree::get_child_at(int i, int j, int k) {

	return this->children[i][j][k];	

}

void Octree::set_child_null_at(int i, int j, int k) {

	this->children[i][j][k] = NULL;

}

Octree* Octree::get_parent() {

	return this->parent;

}

void Octree::set_to_refine_with_nesting() {

	//set refinement criteria
	this->setToRefine = true;
	
	//check nesting
	if(!(this->isRootNode())) {

	
		if(this->east==NULL && this->east_bc==NONE) {

			if(!(this->parent->east->setToRefine==true)) { 
				this->parent->east->set_to_refine_with_nesting();	
			}

		}

		if(this->west==NULL && this->west_bc==NONE) {
			
			if(!(this->parent->west->setToRefine==true)) { 
				this->parent->west->set_to_refine_with_nesting();	
			}

		}

		if(this->north==NULL && this->north_bc==NONE) {
			
			if(!(this->parent->north->setToRefine==true)) { 
				this->parent->north->set_to_refine_with_nesting();	
			}

		}

		if(this->south==NULL && this->south_bc==NONE) {
			
			if(!(this->parent->south->setToRefine==true)) { 
				this->parent->south->set_to_refine_with_nesting();	
			}

		}
		
		if(this->top==NULL && this->top_bc==NONE) {
			
			if(!(this->parent->top->setToRefine==true)) { 
				this->parent->top->set_to_refine_with_nesting();	
			}

		}
		
		if(this->bottom==NULL && this->bottom_bc==NONE) {
			
			if(!(this->parent->bottom->setToRefine==true)) { 
				this->parent->bottom->set_to_refine_with_nesting();	
			}

		}

	}
}	

//checks if all of its neighbours are leaf nodes before setting the criteria
void Octree::set_to_coarsen_with_nesting() {
	
	bool criteria;	
	criteria = true;

	if(this->east) {
		if(!(this->east->isLeafNode())) {
			criteria = false;
		}
	}		
	if(this->west) {
		if(!(this->west->isLeafNode())) {
			criteria = false;
		}
	}		
	if(this->north) {
		if(!(this->north->isLeafNode())) {
			criteria = false;
		}
	}		
	if(this->south) {
		if(!(this->south->isLeafNode())) {
			criteria = false;
		}
	}		
	if(this->top) {
		if(!(this->top->isLeafNode())) {
			criteria = false;
		}
	}		
	if(this->bottom) {
		if(!(this->bottom->isLeafNode())) {
			criteria = false;
		}
	}
	/*change loop indices*/
	//setting criteria to siblings
//      	for(int n=0; n<2; n++) {
//           	for(int m=0; m<2; m++) {
//                    	for(int l=0; l<2; l++) {
//                             	this->get_parent()->get_child_at(l, m, n)->setToCoarsen = criteria;
//                 	}
//         	}
// 	}	

	this->setToCoarsen = criteria;

}

/*comment*/
void Octree::coarsen() {
    
    
        
        //delete this;
   //     nodes.erase(this);
    
	delete this;
    
    //        moving to its parent
    //        this->current = this->parent;
    //
    //        deleting all children
    //        for(int i=0; i<2; i++) {
    //                for(int j=0; j<2; j++) {
    //                        for(int k=0; k<2; k++)
    //                                this->current->children[i][j][k] = NULL;
    //                }
    //        }
    //
    //        delete this;
    
}

void Octree::refine() {
    
    int i, j, k;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double lev = this->level + 1;
    
    //creating children nodes
    for(k=0; k<2; k++) {
        for(j=0; j<2; j++) {
            for(i=0; i<2; i++) {
                
                if(i==0) {
                    xmin = this->x_min;
                    xmax = this->x_centre;
                }
                else {
                    xmin = this->x_centre;
                    xmax = this->x_max;
                }
                
                if(j==0) {
                    ymin = this->y_min;
                    ymax = this->y_centre;
                }
                else {
                    ymin = this->y_centre;
                    ymax = this->y_max;
                }
                
                if(k==0) {
                    zmin = this->z_min;
                    zmax = this->z_centre;
                }
                else {
                    zmin = this->z_centre;
                    zmax = this->z_max;
                }
                
                //creating new memory locations to children
                this->children[i][j][k] = new Octree;
                
                //creating new child object
                Octree child(xmin, xmax, ymin, ymax, zmin, zmax, lev);
                
                //deleting the pushed node to the list
                nodes.pop_back();
                
                //assigning parent to the child
                child.parent = this;
                
                //invoking copy constructor
                *(this->children[i][j][k]) = child;


            }
        }
    }
    
	//setting neighbours to children
    	for(k=0; k<2; k++) {
        	for(j=0; j<2; j++) {
            		for(i=0; i<2; i++) {

				if(i==0) {
					this->children[i][j][k]->east = this->children[i+1][j][k];
					if(this->west != NULL) 	this->children[i][j][k]->west = this->west->children[i+1][j][k];
				}	
				if(j==0) {
					this->children[i][j][k]->north = this->children[i][j+1][k];
					if(this->south != NULL)	this->children[i][j][k]->south = this->south->children[i][j+1][k];
				}	
				if(k==0) {
					this->children[i][j][k]->top = this->children[i][j][k+1];
					if(this->bottom != NULL)	this->children[i][j][k]->bottom = this->bottom->children[i][j][k+1];
				}	
				if(i==1) {
					this->children[i][j][k]->west = this->children[i-1][j][k];
					if(this->east != NULL)	this->children[i][j][k]->east = this->east->children[i-1][j][k];
				}	
				if(j==1) {
					this->children[i][j][k]->south = this->children[i][j-1][k];
					if(this->north != NULL)	this->children[i][j][k]->north = this->north->children[i][j-1][k];
				}	
				if(k==1) {
					this->children[i][j][k]->bottom = this->children[i][j][k-1];
					if(this->top != NULL)	this->children[i][j][k]->top = this->top->children[i][j][k-1];
				}	
			}
		}	
	}		

	//setting boundaries to children
    	for(k=0; k<2; k++) {
        	for(j=0; j<2; j++) {
            		for(i=0; i<2; i++) {

				if(i==0) {
					this->children[i][j][k]->east_bc = NONE;
					this->children[i][j][k]->west_bc = this->west_bc;
				}	
				if(j==0) {
					this->children[i][j][k]->north_bc = NONE;
					this->children[i][j][k]->south_bc = this->south_bc;
				}	
				if(k==0) {
					this->children[i][j][k]->top_bc = NONE;
					this->children[i][j][k]->bottom_bc = this->bottom_bc;
				}	
				if(i==1) {
					this->children[i][j][k]->west_bc = NONE;
					this->children[i][j][k]->east_bc = this->east_bc;
				}	
				if(j==1) {
					this->children[i][j][k]->south_bc = NONE;
					this->children[i][j][k]->north_bc = this->north_bc;
				}	
				if(k==1) {
					this->children[i][j][k]->bottom_bc = NONE;
					this->children[i][j][k]->top_bc = this->top_bc;
				}	
			}
		}	
	}		
		

}

/*check what is happening - destroy is necessary*/
Octree::~Octree() {
    
    //	for(int i=0; i<2; i++) {
    //                for(int j=0; j<2; j++) {
    //                        for(int k=0; k<2; k++)
    //                                delete children[i][j][k];
    //                }
    //        }
    //
    //
    //
    //        delete [] children;
    //
    //        delete block_data;
    //
    //        delete parent;
    //
    //        delete current;
    //
    //        delete [] siblings;
}

/*copy constructor not working*/
Octree::Octree(const Octree &obj) {
    
	std::cerr << " copy constructor of octree is working " << std::endl;

    x_min = obj.x_min;
    y_min = obj.y_min;
    z_min = obj.z_min;
    x_max = obj.x_max;
    y_max = obj.y_max;
    z_max = obj.z_max;
    x_centre = obj.x_centre;
    y_centre = obj.y_centre;
    z_centre = obj.z_centre;
    level = obj.level;
    parent = obj.parent;
    block_data = obj.block_data;
	setToRefine = obj.setToRefine;
	setToCoarsen = obj.setToCoarsen;

	number = obj.number;

	//neighbors
	east = obj.east;
	west = obj.west;
	north = obj.north;
	south = obj.south;
	top = obj.top;
	bottom = obj.bottom; 
	
	//boundary conditions
	east_bc = obj.east_bc;	   
	west_bc = obj.west_bc;	   
	north_bc = obj.north_bc;	   
	south_bc = obj.south_bc;	   
	top_bc = obj.top_bc;	   
	bottom_bc = obj.bottom_bc;	   
	
	//boundary values
	east_bc_val = obj.east_bc_val;	   
	west_bc_val = obj.west_bc_val;	   
	north_bc_val = obj.north_bc_val;	   
	south_bc_val = obj.south_bc_val;	   
	top_bc_val = obj.top_bc_val;	   
	bottom_bc_val = obj.bottom_bc_val;	   
 
    memcpy(children,obj.children,sizeof(Octree***)*2);
    for(int i=0;i<2;i++) {
        memcpy(children[i],obj.children[i],sizeof(Octree**)*2);
        for(int j=0;j<2;j++) {
            memcpy(children[i][j],obj.children[i][j],sizeof(Octree*)*2);
        }
    }
}

Octree::Octree( double x1, double x2, double y1, double y2, double z1, double z2, int l ) : x_min(x1), x_max(x2), y_min(y1), y_max(y2), z_min(z1), z_max(z2), level(l)   {
    
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            for(int k=0; k<2; k++)
                children[i][j][k] = NULL;
        }
    }
    
	parent = NULL;
   	setToRefine = false;
	setToCoarsen = false; 
	
	//neighbors
	east = NULL;
	west = NULL;
	north = NULL;
	south = NULL;
	top = NULL;
	bottom = NULL; 
   
	//boundary conditions
	east_bc = NONE;
	west_bc = NONE;
	north_bc = NONE;
	south_bc = NONE;
	top_bc = NONE;
	bottom_bc = NONE; 

	number = 0;
 
    x_centre = (x_min + x_max ) / 2.0;
    y_centre = (y_min + y_max ) / 2.0;
    z_centre = (z_min + z_max ) / 2.0;
    
    
    //creating block to assign it to the data
    block_data = new Block;
    Block new_block(x_min, x_max, y_min, y_max, z_min, z_max);
    *block_data = new_block;
    
    //printf("dx=%g, dy=%g, dz=%g \n", block_data->dx, block_data->dy, block_data->dz);
    
    //make current pointer point to the current object
    //      current = this;
    
    nodes.push_back(this);
    
}

Octree::Octree() {
    
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            for(int k=0; k<2; k++)
                children[i][j][k] = NULL;
        }
    }
    
 	parent = NULL;
	setToRefine = false;
	setToCoarsen = false;   

	//neighbors
	east = NULL;
	west = NULL;
	north = NULL;
	south = NULL;
	top = NULL;
	bottom = NULL; 
	
	//boundary conditions
	east_bc_val = 1.0;
	west_bc_val = 1.0;
	north_bc_val = 1.0;
	south_bc_val = 1.0;
	top_bc_val = 1.0;
	bottom_bc_val = 1.0; 


	//boundary conditions
	east_bc = NONE;
	west_bc = NONE;
	north_bc = NONE;
	south_bc = NONE;
	top_bc = NONE;
	bottom_bc = NONE; 

	number = 0;

    //creating block to assign it to the data
    block_data = new Block;
    Block new_block;
    *block_data = new_block;
    
    //make current pointer point to the current object
    //      current = this;
    nodes.push_back(this);
}


}
