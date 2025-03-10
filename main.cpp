#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

 std::cout << std::fixed << std::setprecision(3);

    double dimensions[3] = {1.0, 1.0, 1.0};
    int mesh[3] = {11, 11, 11};
    double dt = 0.1;

    lattice L(dimensions, mesh);

    for(auto& node_ : L.get_nodes()){
	node_->ftom(&L, 1);
}

for(double t=0.; t<50; t+=dt){

for(auto& node_ : L.get_nodes() ){ //  

    node_->ftom(&L, 0);

    node_->calc_eq();

    node_->collision(dt);

    node_->ftom(&L, 1);
}

    L.stream();

for(auto& node_ : L.get_nodes()){
    if(node_->get_wallflag() == 0) node_->update_f();
}

for(auto& node_ : L.get_nodes()){
 
if( node_->get_pos() == std::vector<double>{0.5, 0.5, 0.5} ){
    for(int i=0; i<19; i++){
	std::cout << node_->get_m()[i] << ", ";
	}
        std::cout << std::endl;
      }
/*if( node_->get_pos() == std::vector<double>{0.0, 0.5, 0.5} ){
    for(int i=0; i<19; i++){
        std::cout << node_->get_m()[i] << ", ";
        }
        std::cout << std::endl;
      }
if( node_->get_pos() == std::vector<double>{0.5, 0.8, 0.8} ){
    for(int i=0; i<19; i++){
        std::cout << node_->get_m()[i] << ", ";
        }
        std::cout << std::endl;
      }
if( node_->get_pos() == std::vector<double>{0.0, 0.8, 0.8} ){
    for(int i=0; i<19; i++){
        std::cout << node_->get_m()[i] << ", ";
        }
        std::cout << std::endl;
      }
*/
    }
  }
    
return 0;

}

