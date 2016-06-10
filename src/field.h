namespace myOctree {

///Template class for any scalar field variable in the domain.
/*!
Usage:
Field object(nx,ny,nz);
Parameters are the number of cells along x, y and z (including the padding) which usually remains  same for all the blocks in the mesh.
Padding represents extra layer of cells at all the sides of the block, which acts as ghost cells at the block boundaries or buffer cells at the processor boundaries.
New objects of this class are defined in the constructors of Block class and the pointers to these objects are members of the Block class.
*/
class Field {

        public:
        
	/*! \name Size*/	
	//@{
	int Nx,Ny,Nz;              
        int N;
	//@}	
        double*** val;            /*!<3D array for storing field varible values at the cells*/
        std::string name;       /*!<Field name*/
        FieldBc bc[3][2];       /*!<Field boundary conditions*/

        Field( int N_x, int N_y, int N_z, std::string info ); /*!<Parametrized constructor with initialization fields*/

        Field(); /*!<Default constructor*/

        Field(const Field &obj); /*!<Copy constructor*/

        ~Field(); /*!<Destructor*/

        void set_field(double); /*!<Function to set the given value in the field*/


        private:
        protected:

};


}
