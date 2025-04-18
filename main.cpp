#include "header.hpp"
#include "lattice.hpp"
#include "IB.hpp"
#include "mesh_writer.hpp"

int main(int argc, char** argv){

 std::cout << std::fixed << std::setprecision(1);

    double dimensions[3] = {1.0, 1.0, 1.0};
    int mesh[3] = {21, 21, 21};
    double dt = 1.;

    lattice L(dimensions, mesh);
    IB flag;

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
    if(node_->get_wallflag() != 2) node_->update_f();
    //else if( (node_->get_pos()[0] == 0.) || (node_->get_pos()[0] == dimensions[0]) )node_->update_f();

if(node_->get_pos()[1] == 0.5){
//    for(int i=0; i<19; i++){
        std::cout << node_->get_m()[3]*1000 << " " ;
//        }
        if(node_->get_pos()[2] == dimensions[2])std::cout << std::endl;
      }
}
if(int(t) % 10 == 0){
MeshWriter::writeVTK("test/test_mesh_"+std::to_string(int(t))+".vtk", mesh, L.get_nodes() );
    }
}

//MeshWriter::writeVTK("test/test_mesh.vtk", mesh, L.get_nodes() );

return 0;

}

