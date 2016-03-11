#include "block.h"
namespace myOctree {

//parametrized constructor with initialization fields
VecField::VecField( int N_x, int N_y, int N_z, std::string info ) : Nx(N_x), Ny(N_y), Nz(N_z), name(info) {
                //std::cerr << "parametrized constructor of vecfield is working" << std::endl; 
                
        N = Nx*Ny*Nz;
        x = new double** [Nx];
        y = new double** [Nx];
        z = new double** [Nx];
        for(int i=0;i<Nx;i++) {
                x[i] = new double* [Ny];
                y[i] = new double* [Ny];
                z[i] = new double* [Ny];
                for(int j=0;j<Ny;j++) {
                        x[i][j] = new double [Nz];
                        y[i][j] = new double [Nz];
                        z[i][j] = new double [Nz];
                }       
        }

	set_field(0.0);
       
}       

//default constructor
VecField::VecField() {
 //             std::cerr << "default constructor of vecfield is working" << std::endl;
 
        Nx = 0;
        Ny = 0;
        Nz = 0;
        N = Nx*Ny*Nz;
        x = new double** [Nx];
        y = new double** [Nx];
        z = new double** [Nx];
        for(int i=0;i<Nx;i++) {
                x[i] = new double* [Ny];
                y[i] = new double* [Ny];
                z[i] = new double* [Ny];
                for(int j=0;j<Ny;j++) {
                        x[i][j] = new double [Nz];
                        y[i][j] = new double [Nz];
                        z[i][j] = new double [Nz];
                }
        }
}

 //Copy constructor
VecField::VecField(const VecField &obj) {

//              std::cerr << "copy constructor of vecfield is working" << std::endl;

        Nx = obj.Nx;
        Ny = obj.Ny;
        Nz = obj.Nz;
        N = obj.N;
        name = obj.name;

        x = new double** [Nx];
        y = new double** [Nx];
        z = new double** [Nx];
        for(int i=0;i<Nx;i++) {
                x[i] = new double* [Ny];
                y[i] = new double* [Ny];
                z[i] = new double* [Ny];
                for(int j=0;j<Ny;j++) {
                        x[i][j] = new double [Nz];
                        y[i][j] = new double [Nz];
                        z[i][j] = new double [Nz];
                        memcpy(x[i][j],obj.x[i][j],sizeof(double)*Nz);
                        memcpy(y[i][j],obj.y[i][j],sizeof(double)*Nz);
                        memcpy(z[i][j],obj.z[i][j],sizeof(double)*Nz);
                }
        }
}

//Destructor
 VecField::~VecField() {
                //std::cerr << "destructor of vecfield is working" << std::endl;

        for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                        delete [] x[i][j];
                        delete [] y[i][j];
                        delete [] z[i][j];
                }
                delete [] x[i];
                delete [] y[i];
                delete [] z[i];
        }

        delete [] x;
        delete [] y;
        delete [] z;
}

//member function
void VecField::set_field(double value) {

        for(int i=0;i<this->Nx;i++) {
                for(int j=0;j<this->Ny;j++) {
                        for(int k=0;k<this->Nz;k++) {
                                this->x[i][j][k] = value;
                                this->y[i][j][k] = value;
                                this->z[i][j][k] = value;
                        }
                }
        }
}

}
