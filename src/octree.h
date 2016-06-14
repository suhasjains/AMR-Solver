#ifndef MYOCTREE_OCTREE_H_
#define MYOCTREE_OCTREE_H_
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include "block.h"
#include "boundary.h"

namespace myOctree {


//class forward declaration
class Octree;
extern std::list<Octree*> nodes;

///Class to store octree datastructure as nodes of the tree.
/*!
Objects of this class are nodes of the tree.
When a new object is created, children are assigned, parent is assigned, data block is assigned and finally the pointer of the object is pushed into the "nodes" list.  

Usage:
Octree object_name(x min, x max, y min, y max, z min, z max, level); 
This creates a node of the tree at the specified level.
This object can be further refined or coarsened using member functions of this class.
*/
class Octree {

	public:
	Octree(); /*!<Default constructor*/
	
	Octree( double x1, double x2, double y1, double y2, double z1, double z2, int l ); /*!<Parametrized constructor with initialization fields*/

        Octree(const Octree &obj); /*!<Copy constructor*/

	~Octree(); /*!<Destructor*/

	void refine(); /*!<Function to add 8 children*/

	bool isLeafNode(); /*!</Function to test if leaf node*/

	bool isRootNode(); /*!<Function to test if root node*/

	/*change it to containsPoint*/
	bool contains(double x, double y, double z); /*!</Function to test if the point lies inside the block*/

	int get_level(); /*!<Function that returns the level of the node*/	

	Octree* get_child_at(int, int, int); /*!<Function that returns pointer to the child at indices (relative)*/

	void get_relative_location(int*, int*, int*); /*!Returns relative location w.r.t siblings*/

	Octree* get_parent(); /*!<Function that returns pointer to the parent*/

	void set_child_null_at(int, int, int); /*!<Function that sets child pointer at indices (relative) to NULL*/	

	Block* get_block_data(); /*!<Function to access block data*/

	void set_to_refine_with_nesting(); /*!<Sets on refine flag with nesting*/ 
	
	void set_to_coarsen_with_nesting(); /*!<Sets on coarsen flag with nesting*/ 

	bool setToRefine; /*!<Refinement flag*/

	bool setToCoarsen; /*!<Coarsen flag*/
	
	//neighbors
	/*change this to neighbour[direction][position]*/
//	Octree *east;
//	Octree *west;
//	Octree *north;
//	Octree *south;
//	Octree *top;
//	Octree *bottom;

	Octree *neighbour[3][2];

	//boundary conditions
	/*same as above */
//	NodeBc east_bc;
//	NodeBc west_bc;
//	NodeBc north_bc;
//	NodeBc south_bc;
//	NodeBc top_bc;
//	NodeBc bottom_bc;

	NodeBc bc[3][2]; /*!<Node boundary conditions*/
	
	//boundary values
	/*same as above */
//	double east_bc_val;
//	double west_bc_val;
//	double north_bc_val;
//	double south_bc_val;
//	double top_bc_val;
//	double bottom_bc_val;
	
	/*! \name Node's block dimensions*/
	//@{
	double x_centre, y_centre, z_centre;
	double x_min, x_max;
	double y_min, y_max;
	double z_min, z_max;
	//@}


	int number;

	
	private:
	Octree *children[2][2][2]; /*!<Each node has upto 8 children (2^3 for 3 dimensions)*/
	
	Block *block_data; /*!<Pointer to the block data*/

	Octree *parent; /*!<Pointer to the parent*/

	int level; /*!<Level in the tree*/

	
	protected:	

};


}
#endif
