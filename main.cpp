#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

    double dimensions[3] = {1.0, 1.0, 1.0};
    int mesh_size[3] = {2, 2, 2};

    lattice L(dimensions, mesh_size);

for(auto& node_ : L.get_nodes() ){ //  int i=0; i<mesh_size[0]*mesh_size[1]*mesh_size[2]; i++){
    //std::cout << std::fixed << std::setprecision(1) << node_->get_pos()[0] << ", " << node_->get_pos()[1] << ", " << node_->get_pos()[2] << std::endl;

    node_->generate_moment(&L);
for(int i=0; i<19; i++){
    std::cout << node_->get_m()[i] << ", "; // test->generate_moment(L.evector**);
  }
    std::cout << std::endl;
}
    return 0;

}


