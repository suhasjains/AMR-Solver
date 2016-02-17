#include "octree.h"


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
    
    //	//assigning each child its siblings
    //	int local_sibling_count;
    //	for(i=0; i<8; i++) {
    //		local_sibling_count = 0;
    //		for(j=0; j<8; j++) {
    //			if(i==j) continue;
    //			node[i].siblings[local_sibling_count] = &node[j];
    //			local_sibling_count++;
    //		}
    //	}
    //
    //

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

Octree::Octree(const Octree &obj) {
    
    
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
	east_bc = NONE;
	west_bc = NONE;
	north_bc = NONE;
	south_bc = NONE;
	top_bc = NONE;
	bottom_bc = NONE; 

    //creating block to assign it to the data
    block_data = new Block;
    Block new_block;
    *block_data = new_block;
    
    //make current pointer point to the current object
    //      current = this;
    nodes.push_back(this);
}


}
