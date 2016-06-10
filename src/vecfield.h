namespace myOctree {

///Template class for any vector field variable in the domain.
/*!
x, y, z - 3d arrays for storing field variable values at the cells along x, y and z directions.

Usage:
Field object(nx,ny,nz);
Parameters are the number of cells along x, y and z (including the padding) which usually remains  same for al    l the blocks in the mesh.
Padding represents extra layer of cells at all the sides of the block, which acts as ghost cells at the block     boundaries or buffer cells at the processor boundaries.
New objects of this class are defined in the constructors of Block class and the pointers to these objects are     members of the Block class.
*/
class VecField {

        public:
	
	/*!\name Size*/
        //@{
	int Nx,Ny,Nz;              
        int N;                  
        //@}
	
	/*!\name 3D array to store x,y,z comonent of field values*/
	//@{
	double*** x;             
        double*** y;             
        double*** z;             
        //@}
	
	
	std::string name; 	/*!<Name of the field*/

	/*!\name Boundary conditions*/
	//@{
	FieldBc xbc[3][2];       
	FieldBc ybc[3][2];       
	FieldBc zbc[3][2];       
	//@}

        VecField( int N_x, int N_y, int N_z, std::string info ); /*!<Parametrized constructor*/ 

        VecField(); /*!<Default constructor*/

        VecField(const VecField &obj); /*!<Copy constructor*/

        ~VecField(); /*!<Destructor*/

        void set_field(double); /*!<Function to set the given value to the field*/


        private:
        protected:

};


}
