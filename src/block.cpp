#include "block.h"
#include <stdio.h>

namespace myOctree {

//parametrized constructor with initialization fields
Field::Field( int N_x, int N_y, int N_z ) : Nx(N_x), Ny(N_y), Nz(N_z) {

        N = Nx*Ny*Nz;
        //val = new double [10][10][10];
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


        Nx = obj.Nx;
        Ny = obj.Ny;
        Nz = obj.Nz;
        N = obj.N;
        memcpy(val,obj.val,sizeof(double**)*Nx);
        for(int i=0;i<Nx;i++) {
                memcpy(val[i],obj.val[i],sizeof(double*)*Ny);
                for(int j=0;j<Ny;j++) {
                        memcpy(val[i][j],obj.val[i][j],sizeof(double)*Nz);
                }
        }
}

//Destructor
 Field::~Field() {

//	for (int i = 0; i < Nx; ++i) {
//	        for (int j = 0; j < Ny; ++j)
//	        delete [] val[i][j];
//	
//	        delete [] val[i];
//	        }
//	
//	delete [] val;
}

//parametrized constructor with initialization fields
Block::Block( double x1, double x2, double y1, double y2, double z1, double z2 ) : x_min(x1), x_max(x2), y_min(y1), y_max(y2), z_min(z1), z_max(z2) {

        dx = ( x_max - x_min ) / iNx;
        dy = ( y_max - y_min ) / iNy;
        dz = ( z_max - z_min ) / iNz;

        //printf("dx=%g, dy=%g, dz=%g \n", dx, dy, dz);

        x_centre = (x_min + x_max ) / 2.0;
        y_centre = (y_min + y_max ) / 2.0;
        z_centre = (z_min + z_max ) / 2.0;

        //dynamical allocation of the object
        mesh = new Field;
        Field mesh_field(iNx+2*PAD,iNy+2*PAD,iNz+2*PAD);
        *mesh = mesh_field;

        //printf("N=%d\n",mesh->N);

}

//default constructor
Block::Block() {

        mesh = new Field;
        Field mesh_field(iNx+2*PAD,iNy+2*PAD,iNz+2*PAD);
        *mesh = mesh_field;
}

//Copy constructor
Block::Block(const Block &obj) {

        x_centre = obj.x_centre;
        y_centre = obj.y_centre;
        z_centre = obj.z_centre;
        x_min = obj.x_min;
        y_min = obj.y_min;
        x_min = obj.x_min;
        x_max = obj.x_max;
        y_max = obj.y_max;
        z_max = obj.z_max;
        dx = obj.dx;
        dy = obj.dy;
        dz = obj.dx;
        iNx = obj.iNx;
        iNy = obj.iNy;
        iNz = obj.iNz;
        mesh = obj.mesh;

}

//Destructor
Block::~Block() {

//        delete mesh;

}

//member function
void Block::calculate_grid_size() {

        this->dx = ( this->x_max - this->x_min ) / this->iNx;
        this->dy = ( this->y_max - this->y_min ) / this->iNy;
        this->dz = ( this->z_max - this->z_min ) / this->iNz;
}

//member function
void Block::calculate_centre() {
        this->x_centre = (this->x_min + this->x_max ) / 2.0;
        this->y_centre = (this->y_min + this->y_max ) / 2.0;
        this->z_centre = (this->z_min + this->z_max ) / 2.0;
}


}
