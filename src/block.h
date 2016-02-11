#include <string.h>

namespace myOctree {

#define PAD 2
#define NX_BLOCK 10
#define NY_BLOCK 10
#define NZ_BLOCK 10



//FIELD CLASS
//This is a template class for any scalar field variable in the domain.
//val - 3d array for storing field variable values at the cells.
//Use three objects of this class to build a vector field.
//
//Usage:
//Field object(nx,ny,nx);
//Parameters are the number of cells along x, y and z (including the padding) which usually remains  same for all the blocks in the mesh.
//Padding represents extra layer of cells at all the sides of the block, which acts as ghost cells at the block boundaries or buffer cells at the processor boundaries.
//New objects of this class are defined in the constructors of Block class and the pointers to these objects are members of the Block class.
class Field {

        public:
        //Members
        int Nx,Ny,Nz;              //size
        int N;                  //size
        double*** val;            //values

        //Constructors
        //allocates memory to the field variables equal to the number of cells in the domain
        Field( int N_x, int N_y, int N_z );

	//default constructor
        Field();
 	
 	//Copy constructor
        Field(const Field &obj);

        //Destructor
        ~Field();


        private:
        protected:
	
};



//BLOCK CLASS
//This class is a generic data block of a octree node in the block-based AMR mesh.
//This class contains all the fields making the domain.
//
//Usage:
//Block object(xmin, xmax, ymin, ymax, zmin, zmax );.
//Parameters are the boundaries of the block.  
class Block {

	public:
	
	Field *mesh;
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
		
	//member function
	void calculate_grid_size();

	//member function
	void calculate_centre();
};


}
