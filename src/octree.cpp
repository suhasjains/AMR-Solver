#include "octree.h"
#include <iostream>
#include "direction.h"

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

	
	//	if(this->east==NULL && this->east_bc==NONE) {

	//		if(!(this->parent->east->setToRefine==true)) { 
	//			this->parent->east->set_to_refine_with_nesting();	
	//		}

	//	}

	//	if(this->west==NULL && this->west_bc==NONE) {
	//		
	//		if(!(this->parent->west->setToRefine==true)) { 
	//			this->parent->west->set_to_refine_with_nesting();	
	//		}

	//	}

	//	if(this->north==NULL && this->north_bc==NONE) {
	//		
	//		if(!(this->parent->north->setToRefine==true)) { 
	//			this->parent->north->set_to_refine_with_nesting();	
	//		}

	//	}

	//	if(this->south==NULL && this->south_bc==NONE) {
	//		
	//		if(!(this->parent->south->setToRefine==true)) { 
	//			this->parent->south->set_to_refine_with_nesting();	
	//		}

	//	}
	//	
	//	if(this->top==NULL && this->top_bc==NONE) {
	//		
	//		if(!(this->parent->top->setToRefine==true)) { 
	//			this->parent->top->set_to_refine_with_nesting();	
	//		}

	//	}
	//	
	//	if(this->bottom==NULL && this->bottom_bc==NONE) {
	//		
	//		if(!(this->parent->bottom->setToRefine==true)) { 
	//			this->parent->bottom->set_to_refine_with_nesting();	
	//		}

	//	}

		for(int m=0;m<3;m++) {
                	for(int n=0;n<2;n++) {	
				if(this->neighbour[m][n]==NULL && this->bc[m][n]==NONE) {
                      			if(!(this->parent->neighbour[m][n]->setToRefine==true)) { 
                              			this->parent->neighbour[m][n]->set_to_refine_with_nesting();     
                      			}
              			}
			}
		}	
	}
}	

//checks if all of its neighbours are leaf nodes before setting the criteria
void Octree::set_to_coarsen_with_nesting() {
	
	bool criteria;	
	criteria = true;

//	if(this->east) {
//		if(!(this->east->isLeafNode())) {
//			criteria = false;
//		}
//	}		
//	if(this->west) {
//		if(!(this->west->isLeafNode())) {
//			criteria = false;
//		}
//	}		
//	if(this->north) {
//		if(!(this->north->isLeafNode())) {
//			criteria = false;
//		}
//	}		
//	if(this->south) {
//		if(!(this->south->isLeafNode())) {
//			criteria = false;
//		}
//	}		
//	if(this->top) {
//		if(!(this->top->isLeafNode())) {
//			criteria = false;
//		}
//	}		
//	if(this->bottom) {
//		if(!(this->bottom->isLeafNode())) {
//			criteria = false;
//		}
//	}

	for(int m=0;m<3;m++) {
		for(int n=0;n<2;n++) {
			if(this->neighbour[m][n]) {
		//std::cerr << "\n" << "hi" << std::endl;
              			if(!(this->neighbour[m][n]->isLeafNode())) {
                      			criteria = false;
              			}
      			}     	
		}
	}	

	this->setToCoarsen = criteria;

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

                
                //creating new child object
                Octree child(xmin, xmax, ymin, ymax, zmin, zmax, lev);
                
                //deleting the pushed node to the list
                nodes.pop_back();
                
                //assigning parent to the child
                child.parent = this;
                
                //invoking copy constructor
                this->children[i][j][k] = new Octree(child);
                //*(this->children[i][j][k]) = child;


            }
        }
    }
    
	//setting neighbours to children
    	for(k=0; k<2; k++) {
        	for(j=0; j<2; j++) {
            		for(i=0; i<2; i++) {

			//	if(i==0) {
			//		this->children[i][j][k]->east = this->children[i+1][j][k];
			//		if(this->west != NULL) 	this->children[i][j][k]->west = this->west->children[i+1][j][k];
			//	}	
			//	if(j==0) {
			//		this->children[i][j][k]->north = this->children[i][j+1][k];
			//		if(this->south != NULL)	this->children[i][j][k]->south = this->south->children[i][j+1][k];
			//	}	
			//	if(k==0) {
			//		this->children[i][j][k]->top = this->children[i][j][k+1];
			//		if(this->bottom != NULL)	this->children[i][j][k]->bottom = this->bottom->children[i][j][k+1];
			//	}	
			//	if(i==1) {
			//		this->children[i][j][k]->west = this->children[i-1][j][k];
			//		if(this->east != NULL)	this->children[i][j][k]->east = this->east->children[i-1][j][k];
			//	}	
			//	if(j==1) {
			//		this->children[i][j][k]->south = this->children[i][j-1][k];
			//		if(this->north != NULL)	this->children[i][j][k]->north = this->north->children[i][j-1][k];
			//	}	
			//	if(k==1) {
			//		this->children[i][j][k]->bottom = this->children[i][j][k-1];
			//		if(this->top != NULL)	this->children[i][j][k]->top = this->top->children[i][j][k-1];
			//	}	
				if(i==0) {
					this->children[i][j][k]->neighbour[XDIR][RIGHT] = this->children[i+1][j][k];
					if(this->neighbour[XDIR][LEFT] != NULL) 	this->children[i][j][k]->neighbour[XDIR][LEFT] = this->neighbour[XDIR][LEFT]->children[i+1][j][k];
				}	
				if(j==0) {
					this->children[i][j][k]->neighbour[YDIR][RIGHT] = this->children[i][j+1][k];
					if(this->neighbour[YDIR][LEFT] != NULL)	this->children[i][j][k]->neighbour[YDIR][LEFT] = this->neighbour[YDIR][LEFT]->children[i][j+1][k];
				}	
				if(k==0) {
					this->children[i][j][k]->neighbour[ZDIR][RIGHT] = this->children[i][j][k+1];
					if(this->neighbour[ZDIR][LEFT] != NULL)	this->children[i][j][k]->neighbour[ZDIR][LEFT] = this->neighbour[ZDIR][LEFT]->children[i][j][k+1];
				}	
				if(i==1) {
					this->children[i][j][k]->neighbour[XDIR][LEFT] = this->children[i-1][j][k];
					if(this->neighbour[XDIR][RIGHT] != NULL)	this->children[i][j][k]->neighbour[XDIR][RIGHT] = this->neighbour[XDIR][RIGHT]->children[i-1][j][k];
				}	
				if(j==1) {
					this->children[i][j][k]->neighbour[YDIR][LEFT] = this->children[i][j-1][k];
					if(this->neighbour[YDIR][RIGHT] != NULL)	this->children[i][j][k]->neighbour[YDIR][RIGHT] = this->neighbour[YDIR][RIGHT]->children[i][j-1][k];
				}	
				if(k==1) {
					this->children[i][j][k]->neighbour[ZDIR][LEFT] = this->children[i][j][k-1];
					if(this->neighbour[ZDIR][RIGHT] != NULL)	this->children[i][j][k]->neighbour[ZDIR][RIGHT] = this->neighbour[ZDIR][RIGHT]->children[i][j][k-1];
				}	
			}
		}	
	}		

	//setting boundaries to children
    	for(k=0; k<2; k++) {
        	for(j=0; j<2; j++) {
            		for(i=0; i<2; i++) {

			//	if(i==0) {
			//		this->children[i][j][k]->east_bc = NONE;
			//		this->children[i][j][k]->west_bc = this->west_bc;
			//	}	
			//	if(j==0) {
			//		this->children[i][j][k]->north_bc = NONE;
			//		this->children[i][j][k]->south_bc = this->south_bc;
			//	}	
			//	if(k==0) {
			//		this->children[i][j][k]->top_bc = NONE;
			//		this->children[i][j][k]->bottom_bc = this->bottom_bc;
			//	}	
			//	if(i==1) {
			//		this->children[i][j][k]->west_bc = NONE;
			//		this->children[i][j][k]->east_bc = this->east_bc;
			//	}	
			//	if(j==1) {
			//		this->children[i][j][k]->south_bc = NONE;
			//		this->children[i][j][k]->north_bc = this->north_bc;
			//	}	
			//	if(k==1) {
			//		this->children[i][j][k]->bottom_bc = NONE;
			//		this->children[i][j][k]->top_bc = this->top_bc;
			//	}	
				if(i==0) {
					this->children[i][j][k]->bc[XDIR][RIGHT] = NONE;
					this->children[i][j][k]->bc[XDIR][LEFT] = this->bc[XDIR][LEFT];
				}	
				if(j==0) {
					this->children[i][j][k]->bc[YDIR][RIGHT] = NONE;
					this->children[i][j][k]->bc[YDIR][LEFT] = this->bc[YDIR][LEFT];
				}	
				if(k==0) {
					this->children[i][j][k]->bc[ZDIR][RIGHT] = NONE;
					this->children[i][j][k]->bc[ZDIR][LEFT] = this->bc[ZDIR][LEFT];
				}	
				if(i==1) {
					this->children[i][j][k]->bc[XDIR][LEFT] = NONE;
					this->children[i][j][k]->bc[XDIR][RIGHT] = this->bc[XDIR][RIGHT];
				}	
				if(j==1) {
					this->children[i][j][k]->bc[YDIR][LEFT] = NONE;
					this->children[i][j][k]->bc[YDIR][RIGHT] = this->bc[YDIR][RIGHT];
				}	
				if(k==1) {
					this->children[i][j][k]->bc[ZDIR][LEFT] = NONE;
					this->children[i][j][k]->bc[ZDIR][RIGHT] = this->bc[ZDIR][RIGHT];
				}	
			}
		}	
	}		

	//setting boundary conditions to children fields 
    	for(k=0; k<2; k++) {
        	for(j=0; j<2; j++) {
            		for(i=0; i<2; i++) {
		
				for(int l = 0; l<scalar_fields.size() ; l++) {
					for (int m=0; m<3; m++) {
						for (int n=0; n<2; n++) {
					
							if(this->children[i][j][k]->bc[m][n] == BOUNDARY) {
								this->children[i][j][k]->get_block_data()->scalarfields[l]->bc[m][n] = \
								this->get_block_data()->scalarfields[l]->bc[m][n];  			
							}
							
							if(this->children[i][j][k]->bc[m][n] == NONE) {
								this->children[i][j][k]->get_block_data()->scalarfields[l]->bc[m][n] = none;
							}
							
							if(this->children[i][j][k]->bc[m][n] == MPI_BOUNDARY) {
								this->children[i][j][k]->get_block_data()->scalarfields[l]->bc[m][n] = mpi_boundary;
							}
						}	
					}
				}
				
				for(int l = 0; l<vector_fields.size() ; l++) {
					for (int m=0; m<3; m++) {
						for (int n=0; n<2; n++) {
					
							if(this->children[i][j][k]->bc[m][n] == BOUNDARY) {
								this->children[i][j][k]->get_block_data()->vectorfields[l]->xbc[m][n] = \
								this->get_block_data()->vectorfields[l]->xbc[m][n];  			
								this->children[i][j][k]->get_block_data()->vectorfields[l]->ybc[m][n] = \
								this->get_block_data()->vectorfields[l]->ybc[m][n];  			
								this->children[i][j][k]->get_block_data()->vectorfields[l]->zbc[m][n] = \
								this->get_block_data()->vectorfields[l]->zbc[m][n];  			
							}
							
							if(this->children[i][j][k]->bc[m][n] == NONE) {
								this->children[i][j][k]->get_block_data()->vectorfields[l]->xbc[m][n] = none;
								this->children[i][j][k]->get_block_data()->vectorfields[l]->ybc[m][n] = none;
								this->children[i][j][k]->get_block_data()->vectorfields[l]->zbc[m][n] = none;
							}
							
							if(this->children[i][j][k]->bc[m][n] == MPI_BOUNDARY) {
								this->children[i][j][k]->get_block_data()->vectorfields[l]->xbc[m][n] = mpi_boundary;
								this->children[i][j][k]->get_block_data()->vectorfields[l]->ybc[m][n] = mpi_boundary;
								this->children[i][j][k]->get_block_data()->vectorfields[l]->zbc[m][n] = mpi_boundary;
							}
						}	
					}
				}
			}
		}	
	}		
}

/*check what is happening - destroy is necessary*/
Octree::~Octree() {
//	std::cerr << " destructor of octree is working " << std::endl;
    
            delete block_data;
}

/*copy constructor not working*/
Octree::Octree(const Octree &obj) {
    
//	std::cerr << " copy constructor of octree is working " << std::endl;

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

    	for(int i=0; i<2; i++) {
        	for(int j=0; j<2; j++) {
            		for(int k=0; k<2; k++)
                		children[i][j][k] = obj.children[i][j][k];
        	}
    	}

    	parent = obj.parent;
    	block_data = new Block(*(obj.block_data));
	setToRefine = obj.setToRefine;
	setToCoarsen = obj.setToCoarsen;

	number = obj.number;

	//neighbors
//	east = obj.east;
//	west = obj.west;
//	north = obj.north;
//	south = obj.south;
//	top = obj.top;
//	bottom = obj.bottom; 
	
	for(int m=0;m<3;m++) {
		for(int n=0;n<2;n++) {
			neighbour[m][n] = obj.neighbour[m][n];
		}
	}	
	
	//boundary conditions
//	east_bc = obj.east_bc;	   
//	west_bc = obj.west_bc;	   
//	north_bc = obj.north_bc;	   
//	south_bc = obj.south_bc;	   
//	top_bc = obj.top_bc;	   
//	bottom_bc = obj.bottom_bc;	   

	for(int m=0;m<3;m++) {
		for(int n=0;n<2;n++) {
			bc[m][n] = obj.bc[m][n];
		}
	}	
    	

	nodes.push_back(this);
 
}

Octree::Octree( double x1, double x2, double y1, double y2, double z1, double z2, int l ) : x_min(x1), x_max(x2), y_min(y1), y_max(y2), z_min(z1), z_max(z2), level(l)   {
//	std::cerr << " parametrized constructor of octree is working " << std::endl;
    
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
//	east = NULL;
//	west = NULL;
//	north = NULL;
//	south = NULL;
//	top = NULL;
//	bottom = NULL; 
	
	for(int m=0;m<3;m++) {
		for(int n=0;n<2;n++) {
			neighbour[m][n] = NULL;
		}
	}	
   
	//boundary conditions
//	east_bc = NONE;
//	west_bc = NONE;
//	north_bc = NONE;
//	south_bc = NONE;
//	top_bc = NONE;
//	bottom_bc = NONE; 
	
	for(int m=0;m<3;m++) {
		for(int n=0;n<2;n++) {
			bc[m][n] = NONE;
		}
	}	

	number = 0;
 
    	x_centre = (x_min + x_max ) / 2.0;
    	y_centre = (y_min + y_max ) / 2.0;
    	z_centre = (z_min + z_max ) / 2.0;
    
    
    	//creating block to assign it to the data
    	Block new_block(x_min, x_max, y_min, y_max, z_min, z_max);
    	block_data = new Block(new_block);
    
    	nodes.push_back(this);
    
}

Octree::Octree() {
//	std::cerr << " default constructor of octree is working " << std::endl;
    
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
//	east = NULL;
//	west = NULL;
//	north = NULL;
//	south = NULL;
//	top = NULL;
//	bottom = NULL; 
	
	for(int m=0;m<3;m++) {
		for(int n=0;n<2;n++) {
			neighbour[m][n] = NULL;
		}
	}	
	
	x_centre = 0.0;
	y_centre = 0.0;
	z_centre = 0.0;
	x_min = 0.0;
	y_min = 0.0;
	z_min = 0.0;
	x_max = 0.0;
	y_max = 0.0;
	z_max = 0.0;

	level = 0;

	//boundary conditions
//	east_bc = NONE;
//	west_bc = NONE;
//	north_bc = NONE;
//	south_bc = NONE;
//	top_bc = NONE;
//	bottom_bc = NONE; 
	

	for(int m=0;m<3;m++) {
		for(int n=0;n<2;n++) {
			bc[m][n] = NONE;
		}
	}	

	number = 0;

    //creating block to assign it to the data
    Block new_block;
    block_data = new Block(new_block);
    
    nodes.push_back(this);
}


}
