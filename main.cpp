#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

    double dimensions[3] = {1.0, 1.0, 1.0};
    int mesh_size[3] = {3, 3, 3};

    lattice L(dimensions, mesh_size);

for(auto& node_ : L.get_nodes() ){ //  int i=0; i<mesh_size[0]*mesh_size[1]*mesh_size[2]; i++){
    std::cout << std::fixed << std::setprecision(1) << node_->get_pos()[0] << ", " << node_->get_pos()[1] << ", " << node_->get_pos()[2] << std::endl;

    //std::cout << std::fixed << std::setprecision(1) << L.get_node(i)->get()->get_f()[0] << ", " << L.get_node(i)->get()->get_f()[1] << ", " << L.get_node(i)->get()->get_f()[2] << std::endl;
}
    return 0;

}


