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



	int Npx = NX_BLOCK + 1;
	int Npy = NY_BLOCK + 1;
	int Npz = NZ_BLOCK + 1;
	

	char filename[30];
        sprintf(filename, "output.vtk");
        FILE *fp = fopen(filename, "w");

        fprintf(fp,"# vtk DataFile Version 3.0\n");
        fprintf(fp,"particle point data\n");
        fprintf(fp,"ASCII\n");
        fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	long int nPoints = Npx * Npy * Npz * nodes.size();
	long int nCells = NX_BLOCK * NY_BLOCK * NZ_BLOCK * nodes.size();
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

		printf("dx=%g dy=%g, dz=%g N =%d\n", dx, dy, dz, N);

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

                for(int k = 0; k<NZ_BLOCK; k++) {
                        for(int j = 0; j<NY_BLOCK; j++) {
                                for(int i = 0; i<NX_BLOCK ; i++) {

					point[0] = node_count*Npx*Npy*Npz + get_point(i,j,k,Npx,Npy);
					point[1] = node_count*Npx*Npy*Npz + get_point(i+1,j,k,Npx,Npy);
					point[2] = node_count*Npx*Npy*Npz + get_point(i,j+1,k,Npx,Npy);
					point[3] = node_count*Npx*Npy*Npz + get_point(i+1,j+1,k,Npx,Npy);
					point[4] = node_count*Npx*Npy*Npz + get_point(i,j,k+1,Npx,Npy);
					point[5] = node_count*Npx*Npy*Npz + get_point(i+1,j,k+1,Npx,Npy);
					point[6] = node_count*Npx*Npy*Npz + get_point(i,j+1,k+1,Npx,Npy);
					point[7] = node_count*Npx*Npy*Npz + get_point(i+1,j+1,k+1,Npx,Npy);
	
                                        fprintf(fp,"8 %ld %ld %ld %ld %ld %ld %ld %ld \n",point[0], point[1], point[2], point[3], point[4], point[5], point[6], point[7]);
                                }
                        }
                }

		node_count++;

        }


	fprintf(fp,"\n\nCELL_TYPES %ld\n", nCells);

        for (std::list<Octree*>::iterator iterator = nodes.begin(), end = nodes.end(); iterator != end; ++iterator) {

                for(int k = 0; k<Npz; k++) {
                        for(int j = 0; j<Npy; j++) {
                                for(int i = 0; i<Npx ; i++) {

                                        fprintf(fp,"11\n");
                                }
                        }
                }


        }	



}

}
