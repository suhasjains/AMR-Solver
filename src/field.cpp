#include "block.h"
namespace myOctree {

//parametrized constructor with initialization fields
Field::Field( int N_x, int N_y, int N_z, std::string info ) : Nx(N_x), Ny(N_y), Nz(N_z), name(info)  {
        //      std::cerr << "parametrized constructor of field is working" << std::endl;

        N = Nx*Ny*Nz;
        val = new double** [Nx];
        for(int i=0;i<Nx;i++) {
                val[i] = new double* [Ny];
                for(int j=0;j<Ny;j++) {
                        val[i][j] = new double [Nz];
                }
        }
}

//default constructor
Field::Field() {
        //      std::cerr << "default constructor of field is working" << std::endl;

        Nx = 0;
        Ny = 0;
        Nz = 0;
        N = Nx*Ny*Nz;
        val = new double** [Nx];
        for(int i=0;i<Nx;i++) {
                val[i] = new double* [Ny];
                for(int j=0;j<Ny;j++) {
                        val[i][j] = new double [Nz];
                }
        }
}

//Copy constructor
Field::Field(const Field &obj) {

        //      std::cerr << "copy constructor of field is working" << std::endl;

        Nx = obj.Nx;
        Ny = obj.Ny;
        Nz = obj.Nz;
        N = obj.N;
        name = obj.name;
        //      std::cerr << "field is" << name <<  std::endl;

        val = new double** [Nx];
        for(int i=0;i<Nx;i++) {
                val[i] = new double* [Ny];
                for(int j=0;j<Ny;j++) {
                        val[i][j] = new double [Nz];
                        memcpy(val[i][j],obj.val[i][j],sizeof(double)*Nz);
                }
        }

        for(int i = 0; i < 3; ++ i) {
                memcpy(&(bc[i][0]), &(obj.bc[i][0]), 2 * sizeof(FieldBc));
        }

}

//Destructor
 Field::~Field() {
        //      std::cerr << "delete constructor of field is working" << std::endl;

        for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j)
                        delete [] val[i][j];

                delete [] val[i];
                }

        delete [] val;

}

//member function
void Field::set_field(double value) {

        for(int i=0;i<this->Nx;i++) {
                for(int j=0;j<this->Ny;j++) {
                        for(int k=0;k<this->Nz;k++) {
                                this->val[i][j][k] = value;
                        }
                }
        }
}


}
