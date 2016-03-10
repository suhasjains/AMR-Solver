#include <string.h>
#include <string>
#include <list>
#include <vector>
#include "boundary.h"
#include "field.h"
#include "vecfield.h"

namespace myOctree {

extern int pad;
extern std::vector<std::string> scalar_fields;
extern std::vector<std::string> vector_fields;
extern int nx_block;
extern int ny_block;
extern int nz_block;

//BLOCK CLASS
//This class is a generic data block of a octree node in the block-based AMR mesh.
//This class contains all the fields making the domain.
//
//Usage:
//Block object(xmin, xmax, ymin, ymax, zmin, zmax );.
//Parameters are the boundaries of the block.  
class Block {

	public:
	
	VecField *mesh;
	Field *field;

	/*complete this*/	
	Field **scalarfields;
	VecField **vectorfields;

	//VecField *gradient;
	double max_gradient;
	double x_centre, y_centre, z_centre;	
	double x_min, x_max;
	double y_min, y_max;
	double z_min, z_max;
	static int iNx;
	static int iNy;
	static int iNz;
	double dx, dy, dz;

	//parametrized constructor with initialization fields
	Block( double x1, double x2, double y1, double y2, double z1, double z2 );
		
	//default constructor
	Block();

	//Copy constructor
        Block(const Block &obj);

	//Destructor
	~Block();

};


}
