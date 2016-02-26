#include <string.h>
#include <string>
#include <list>
#include <vector>

namespace myOctree {

extern int pad;
extern std::vector<std::string> scalar_fields;
extern std::vector<std::string> vector_fields;

//FIELD CLASS
//This is a template class for any scalar field variable in the domain.
//val - 3d array for storing field variable values at the cells.
//Use three objects of this class to build a vector field.
//
//Usage:
//Field object(nx,ny,nz);
//Parameters are the number of cells along x, y and z (including the padding) which usually remains  same for all the blocks in the mesh.
//Padding represents extra layer of cells at all the sides of the block, which acts as ghost cells at the block boundaries or buffer cells at the processor boundaries.
//New objects of this class are defined in the constructors of Block class and the pointers to these objects are members of the Block class.
class Field {

        public:
        //Members
        int Nx,Ny,Nz;              //size
        int N;                  //size
        double*** val;            //values
	std::string name;

        //Constructors
        //allocates memory to the field variables equal to the number of cells in the domain
        Field( int N_x, int N_y, int N_z, std::string info );

	//default constructor
        Field();
 	
 	//Copy constructor
        Field(const Field &obj);

        //Destructor
        ~Field();

	//sets field with the input value
	void set_field(double);


        private:
        protected:
	
};

//VECFIELD CLASS
//This is a template class for any vector field variable in the domain.
//x, y, z - 3d arrays for storing field variable values at the cells along x, y and z directions.
//Use three objects of this class to build a vector field.
//
//Usage:
//Field object(nx,ny,nz);
//Parameters are the number of cells along x, y and z (including the padding) which usually remains  same for all the blocks in the mesh.
//Padding represents extra layer of cells at all the sides of the block, which acts as ghost cells at the block boundaries or buffer cells at the processor boundaries.
//New objects of this class are defined in the constructors of Block class and the pointers to these objects are members of the Block class.
class VecField {

        public:
        //Members
        int Nx,Ny,Nz;              //size
        int N;                  //size
        double*** x;            //values
        double*** y;            //values
        double*** z;            //values
	std::string name;

        //Constructors
        //allocates memory to the field variables equal to the number of cells in the domain
        VecField( int N_x, int N_y, int N_z, std::string info );

	//default constructor
        VecField();
 	
 	//Copy constructor
        VecField(const VecField &obj);

        //Destructor
        ~VecField();

	//sets field with the input value
	void set_field(double);


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
