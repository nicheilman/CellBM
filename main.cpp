#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

    std::vector<float> dimensions = {1.0, 1.0, 1.0};
    std::vector<int> mesh_size = {3, 3, 3};

    lattice L(dimensions, mesh_size);

    //std::cout << L.get_f << std::endl;

    return 0;

}


