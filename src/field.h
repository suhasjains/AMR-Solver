namespace myOctree {

//FIELD CLASS
//This is a template class for any scalar field variable in the domain.
//val - 3d array for storing field variable values at the cells.
//Use three objects of this class to build a vector field.
//
//Usage:
//Field object(nx,ny,nz);
//Parameters are the number of cells along x, y and z (including the padding) which usually remains  same for al    l the blocks in the mesh.
//Padding represents extra layer of cells at all the sides of the block, which acts as ghost cells at the block     boundaries or buffer cells at the processor boundaries.
//New objects of this class are defined in the constructors of Block class and the pointers to these objects are     members of the Block class.
class Field {

        public:
        //Members
        int Nx,Ny,Nz;              //size
        int N;                  //size
        double*** val;            //values
        std::string name;       //field name
        FieldBc bc[3][2];       //boundary conditions

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


}
