#include "header.hpp"
#include "lattice.hpp"
#include "IB.hpp"
#include "mesh_writer.hpp"
using namespace std;

int main(int argc, char** argv){

 cout << fixed << setprecision(1);

    int mesh_size = 21;
    double dimensions[3] = {1.0, 1.0, 1.0};
    int mesh[3] = {int(mesh_size*dimensions[0]), int(mesh_size*dimensions[1]), int(mesh_size*dimensions[2])};
    int flag_mesh[3] = {2*mesh[0], 2*mesh[1], 0}; // Must be ~2x the fluid mesh spacing
    double dt = 0.1;
    double force[3] = {0., 0., 0.};
    double mesh_space = 1.0/mesh_size ; 
    double dist, dist_x, dist_y, dist_z, kernel;
 
    lattice L(dimensions, mesh, dt);
    IB flag(flag_mesh);

    for(auto& node_ : L.get_nodes()){
	node_->calc_eq();
	node_->ftom(&L, 1);
	node_->update_f();
}

/* --- Main Loop --- */
for(double t=0.; t<150; t+=dt){

	cout << t << endl;

for(auto& node_ : L.get_nodes() ){ //  

    node_->ftom(&L, 0);

    node_->calc_eq();

force[0] = 0.; force[1] = 0.; force[2] = 0.;

//----------------------------------------------------------------------
if( (node_->get_pos()[2] > 0.3) && (node_->get_pos()[2] < 0.7) ){
for(auto& IB_node_ : flag.get_IB_nodes() ) {
    dist_x = abs(node_->get_pos()[0] - IB_node_->get_IB_pos()[0]) ;
    dist_y = abs(node_->get_pos()[1] - IB_node_->get_IB_pos()[1]) ;
    dist_z = abs(node_->get_pos()[2] - IB_node_->get_IB_pos()[2]) ;


    dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
    dist = sqrt(dist);
    //dist = sqrt( pow(node_->get_pos()[0] - IB_node_->get_pos()[0], 2) + pow(node_->get_pos()[1] - IB_node_->get_pos()[1], 2) + pow(node_->get_pos()[2] - IB_node_->get_pos()[2], 2) );
    if(dist <= 2*mesh_space ){
	kernel = (1 - abs( dist_x/mesh_space/2)) * (1 - abs( dist_y/mesh_space/2)) * (1 - abs( dist_z/mesh_space/2));
	for(int i=0; i<3; i++) force[i] = -2*node_->get_m()[i+1]/node_->get_m()[0] * kernel * dt / mesh_space ;
	}
    }
}
//----------------------------------------------------------------------

    node_->collision(dt, force);

    node_->ftom(&L, 1);
}

    L.stream();

for(auto& node_ : L.get_nodes()){
    node_->update_f();


/* --- Visualize --- */
if(node_->get_pos()[0] == 0.5){
//    for(int i=0; i<19; i++){
        cout << node_->get_m()[3]*1000 << " " ;
//        }
        if(node_->get_pos()[2] == dimensions[2]) cout << endl;
      }
}

/* --- Output --- */
if(int(t) == 0){
MeshWriter::writeflagVTK("test2/test_flag.vtk", flag_mesh, flag.get_IB_nodes() );
    }
if(int(t/dt) % 10 == 0){
MeshWriter::writeVTK("test2/test_mesh_"+to_string(int(t))+".vtk", mesh, L.get_nodes() );
    }
}


return 0;

}

