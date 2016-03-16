
namespace amrsolver {

void jacobi(int, std::string);
void gauss_seidel(int, std::string);
double jacobi_for_field(myOctree::Octree*, myOctree::Field*, double);
double gauss_seidel_red(myOctree::Octree*, myOctree::Field*, double);
double gauss_seidel_black(myOctree::Octree*, myOctree::Field*, double);


}
