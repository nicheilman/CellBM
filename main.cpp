#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

    std::vector<float> dimensions = {1.0, 1.0, 1.0};
    std::vector<int> mesh_size = {3, 3, 3};

    lattice L(dimensions, mesh_size);

for(int i=0; i<mesh_size[0]*mesh_size[1]*mesh_size[2]; i++){
    auto testnode = L.get_node(i);
    std::cout << testnode.get_f()[0] << testnode.get_f()[1] << testnode.get_f()[2] << std::endl;
}
    return 0;

}


