#include "block.h"
#include <stdio.h>
#include <iostream>

namespace myOctree {

int pad = 2;
std::vector<std::string> scalar_fields;
std::vector<std::string> vector_fields;	
int nx_block = 40;
int ny_block = 40;
int nz_block = 20;
int Block::iNx = nx_block;
int Block::iNy = ny_block;
int Block::iNz = nz_block;

//parametrized constructor with initialization fields
Block::Block( double x1, double x2, double y1, double y2, double z1, double z2 ) : x_min(x1), x_max(x2), y_min(y1), y_max(y2), z_min(z1), z_max(z2) {
//  		std::cerr << "Parametrized constructor of block is working" << std::endl;

	max_gradient = 0.0;

        dx = ( x_max - x_min ) / iNx;
        dy = ( y_max - y_min ) / iNy;
        dz = ( z_max - z_min ) / iNz;

	//std::cerr << dz << std::endl;

        //printf("dx=%g, dy=%g, dz=%g \n", dx, dy, dz);

        x_centre = (x_min + x_max ) / 2.0;
        y_centre = (y_min + y_max ) / 2.0;
        z_centre = (z_min + z_max ) / 2.0;

        //dynamical allocation of the objects
        VecField mesh_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, "mesh");
        mesh = new VecField(mesh_field);
        
        Field field_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, "field");
	field = new Field(field_field);

	scalarfields = new Field* [scalar_fields.size()];
	for(int i = 0; i<scalar_fields.size() ; i++) {
        	Field field_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, scalar_fields[i]);
	       	scalarfields[i] = new Field(field_field);	
	}       

	vectorfields = new VecField* [vector_fields.size()];
	for(int i = 0; i<vector_fields.size() ; i++) {
        	VecField vec_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, vector_fields[i]);
	       	vectorfields[i] = new VecField(vec_field);	
	}       

	//std::cerr << "Hi" << z_max << std::endl;
	
	//storing cell centre locations in mesh vector field
	for(int i=0;i<mesh->Nx;i++) {
		for(int j=0;j<mesh->Ny;j++) {
			for(int k=0;k<mesh->Nz;k++) {
				mesh->x[i][j][k] = x_min - pad*dx + dx*(i + 0.5);	 
				mesh->y[i][j][k] = y_min - pad*dy + dy*(j + 0.5);	 
				mesh->z[i][j][k] = z_min - pad*dz + dz*(k + 0.5);	 
			}
		}
	}

}

//default constructor
Block::Block() {
  //		std::cerr << "default constructor of block is working" << std::endl;

        VecField mesh_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, "mesh");
        mesh = new VecField(mesh_field);
	
        Field field_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, "field");
	field = new Field(field_field);
	
	scalarfields = new Field* [scalar_fields.size()];
	for(int i = 0; i<scalar_fields.size() ; i++) {
        	Field field_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, scalar_fields[i]);
	       	scalarfields[i] = new Field(field_field);	
	}       
	
	vectorfields = new VecField* [vector_fields.size()];
	for(int i = 0; i<vector_fields.size() ; i++) {
        	VecField vec_field(iNx+2*pad,iNy+2*pad,iNz+2*pad, vector_fields[i]);
	       	vectorfields[i] = new VecField(vec_field);	
	}       
}

//Copy constructor
Block::Block(const Block &obj) {
  //		std::cerr << "copy constructor of block is working" << std::endl;

        x_centre = obj.x_centre;
        y_centre = obj.y_centre;
        z_centre = obj.z_centre;
        x_min = obj.x_min;
        y_min = obj.y_min;
        z_min = obj.z_min;
        x_max = obj.x_max;
        y_max = obj.y_max;
        z_max = obj.z_max;
        dx = obj.dx;
        dy = obj.dy;
        dz = obj.dz;
        iNx = obj.iNx;
        iNy = obj.iNy;
        iNz = obj.iNz;
	mesh = new VecField(*(obj.mesh));
        field = new Field(*(obj.field));
	max_gradient = obj.max_gradient;

	scalarfields = new Field* [scalar_fields.size()];
        for(int i = 0; i<scalar_fields.size() ; i++) {
                scalarfields[i] = new Field(*((obj.scalarfields)[i]));
        }

        vectorfields = new VecField* [vector_fields.size()];
        for(int i = 0; i<vector_fields.size() ; i++) {
                vectorfields[i] = new VecField(*((obj.vectorfields)[i]));
        }

}

//Destructor
Block::~Block() {
  //		std::cerr << "destructor of block is working" << std::endl;

        delete mesh;
	delete field;

	for (int i = 0; i < scalar_fields.size(); ++i)
                        delete scalarfields[i];

   	delete scalarfields;

	for (int i = 0; i < vector_fields.size(); ++i)
                        delete vectorfields[i];

   	delete vectorfields;

}


}
