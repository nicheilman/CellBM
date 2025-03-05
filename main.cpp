#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

 std::cout << std::fixed << std::setprecision(3);

    double dimensions[3] = {1.0, 1.0, 1.0};
    int mesh_size[3] = {2, 2, 2};

    lattice L(dimensions, mesh_size);

for(auto& node_ : L.get_nodes() ){ //  

    //node_->ftom(&L, 1);


for(int i=0; i<19; i++){
    std::cout << node_->get_m()[i] << ", "; // 
  }
    std::cout << std::endl;
for(int i=0; i<19; i++){
    std::cout << node_->get_f()[i] << ", "; // 
  }
    std::cout << std::endl;

}

    return 0;

}


