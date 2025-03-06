#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

 std::cout << std::fixed << std::setprecision(3);

    double dimensions[3] = {1.0, 1.0, 1.0};
    int mesh_size[3] = {3, 3, 3};

    lattice L(dimensions, mesh_size);

for(int t=0; t<1; t++){

for(auto& node_ : L.get_nodes() ){ //  

    node_->ftom(&L, 0);

    node_->calc_eq();

    node_->collision();

    node_->ftom(&L, 1);
}

    L.stream();

for(auto& node_ : L.get_nodes()){

    node_->update_f();
}

for(auto& node_ : L.get_nodes()){
 
if(!node_->get_wallflag()){
    for(int i=0; i<3; i++){
        std::cout << node_->get_pos()[i] << ", " << ", "; // 
      }
        std::cout << node_->get_wallflag() << ", " << std::endl;
      }
    }
  }
    
return 0;

}


