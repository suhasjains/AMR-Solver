#include "octree.h"
#include <stdio.h>

namespace myOctree {


long int get_point(int i, int j, int k, int Npx, int Npy) {

	int I = i+1;
	int J = j+1;
	int K = k+1;
	
	long int n;

	n = k * Npx * Npy + j * Npx + I; 	

	return n-1;
} 


void write_vtk(std::list<Octree*>& nodes) {



	int Npx = nx_block + 1;
	int Npy = ny_block + 1;
	int Npz = nz_block + 1;
	

	char filename[30];
        sprintf(filename, "output.vtk");
        FILE *fp = fopen(filename, "w");

        fprintf(fp,"# vtk DataFile Version 3.0\n");
        fprintf(fp,"particle point data\n");
        fprintf(fp,"ASCII\n");
        fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	long int nPoints = Npx * Npy * Npz * nodes.size();
	long int nCells = nx_block * ny_block * nz_block * nodes.size();
       	long int point[8]; 
	int node_count = 0;


	fprintf(fp,"POINTS %ld double\n", nPoints);

	for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {
    	
		Block* block_data = (*iterator)->get_block_data();
		int N = block_data->mesh->N;
		double dx = block_data->dx;
		double dy = block_data->dy;
		double dz = block_data->dz;
		double x_min = block_data->x_min;		
		double y_min = block_data->y_min;		
		double z_min = block_data->z_min;		

//		printf("dx=%g dy=%g, dz=%g N =%d\n", dx, dy, dz, N);
//		printf("east bc = %d\n",(*iterator)->east_bc);
//		printf("west bc = %d\n",(*iterator)->west_bc);
//		printf("north bc = %d\n",(*iterator)->north_bc);
//		printf("south bc = %d\n",(*iterator)->south_bc);
//		printf("top bc = %d\n",(*iterator)->top_bc);
//		printf("bottom bc = %d\n",(*iterator)->bottom_bc);

	        for(int k = 0; k<Npz; k++) {
	                for(int j = 0; j<Npy; j++) {
	                        for(int i = 0; i<Npx ; i++) {
	                                fprintf(fp,"%2.8lf %2.8lf %2.8lf\n",x_min + i*dx, y_min + j*dy, z_min + k*dz);
	                        }
	                }
	        }


	}

	fprintf(fp,"\n\nCELLS %ld %ld\n", nCells, 9*nCells);

        for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {

                for(int k = 0; k<nz_block; k++) {
                        for(int j = 0; j<ny_block; j++) {
                                for(int i = 0; i<nx_block ; i++) {

					point[0] = node_count*Npx*Npy*Npz + get_point(i,j,k,Npx,Npy);
					point[1] = node_count*Npx*Npy*Npz + get_point(i+1,j,k,Npx,Npy);
					point[2] = node_count*Npx*Npy*Npz + get_point(i+1,j+1,k,Npx,Npy);
					point[3] = node_count*Npx*Npy*Npz + get_point(i,j+1,k,Npx,Npy);
					point[4] = node_count*Npx*Npy*Npz + get_point(i,j,k+1,Npx,Npy);
					point[5] = node_count*Npx*Npy*Npz + get_point(i+1,j,k+1,Npx,Npy);
					point[6] = node_count*Npx*Npy*Npz + get_point(i+1,j+1,k+1,Npx,Npy);
					point[7] = node_count*Npx*Npy*Npz + get_point(i,j+1,k+1,Npx,Npy);
	
                                        fprintf(fp,"8 %ld %ld %ld %ld %ld %ld %ld %ld \n",point[0], point[1], point[2], point[3], point[4], point[5], point[6], point[7]);
                                }
                        }
                }

		node_count++;

        }


	fprintf(fp,"\n\nCELL_TYPES %ld\n", nCells);

        for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {

                for(int k = 0; k<nz_block; k++) {
                        for(int j = 0; j<ny_block; j++) {
                                for(int i = 0; i<nx_block; i++) {

                                        fprintf(fp,"12\n");
                                }
                        }
                }


        }	

	fprintf(fp,"\n\nCELL_DATA %ld\n", nCells);
	fprintf(fp,"SCALARS input_field_data double 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");

        for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {


		Block* block_data = (*iterator)->get_block_data();
	

                for(int k = pad; k<(nz_block+pad); k++) {
                        for(int j = pad; j<(ny_block+pad); j++) {
                                for(int i = pad; i<(nx_block+pad); i++) {
	
                                        fprintf(fp,"%lf \n",block_data->field->val[i][j][k]);
                                }
                        }
                }

        }


	for (int f = 0 ; f < scalar_fields.size() ; f++) {
		fprintf(fp,"\nSCALARS %s double 1\n",scalar_fields[f].c_str());
		fprintf(fp,"LOOKUP_TABLE default\n");
	
	        for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {
	
	
			Block* block_data = (*iterator)->get_block_data();
		
	
	                for(int k = pad; k<(nz_block+pad); k++) {
	                        for(int j = pad; j<(ny_block+pad); j++) {
	                                for(int i = pad; i<(nx_block+pad); i++) {
		
	                                        fprintf(fp,"%lf \n",block_data->scalarfields[f]->val[i][j][k]);
	                                }
	                        }
	                }
	
	        }
	}
	
	fprintf(fp,"\nVECTORS location_data double\n");
	//fprintf(fp,"LOOKUP_TABLE default\n");

        for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {


		Block* block_data = (*iterator)->get_block_data();
	

                for(int k = pad; k<(nz_block+pad); k++) {
                        for(int j = pad; j<(ny_block+pad); j++) {
                                for(int i = pad; i<(nx_block+pad); i++) {
	
                                        fprintf(fp,"%lf %lf %lf\n",block_data->mesh->x[i][j][k], block_data->mesh->y[i][j][k], block_data->mesh->z[i][j][k]);
                                }
                        }
                }

        }
	
	for (int f = 0 ; f < vector_fields.size() ; f++) {
		fprintf(fp,"\nVECTORS %s double\n",vector_fields[f].c_str());
		//fprintf(fp,"LOOKUP_TABLE default\n");
	
	        for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {
	
	
			Block* block_data = (*iterator)->get_block_data();
		
	
	                for(int k = pad; k<(nz_block+pad); k++) {
	                        for(int j = pad; j<(ny_block+pad); j++) {
	                                for(int i = pad; i<(nx_block+pad); i++) {
		
	                                        fprintf(fp,"%lf %lf %lf\n",block_data->vectorfields[f]->x[i][j][k], block_data->vectorfields[f]->y[i][j][k], block_data->vectorfields[f]->z[i][j][k]);
	                                }
	                        }
	                }
	
	        }
	}

}

}
