#include "header.hpp"
#include "lattice.hpp"
#include "IB.hpp"
#include "mesh_writer.hpp"
using namespace std;

int main(int argc, char** argv){

 cout << fixed << setprecision(1);

 // Need a function to read in these parameters // 
    int mesh_size = 31;
    double dimensions[3] = {1.0, 1.0, 1.0};
    double dt = 0.1;
 // ------------------------------------------- //

    int mesh[3] = {int(mesh_size*dimensions[0]), int(mesh_size*dimensions[1]), int(mesh_size*dimensions[2])};
    int flag_mesh[3] = {2*mesh[0], 2*mesh[1], 0}; // Must be ~2x the fluid mesh spacing
    double force[3] = {0., 0., 0.};
    double mesh_space = 1.0/mesh_size ; 
    double dist, dist_x, dist_y, dist_z, kernel;
    bool internal;
 
    lattice L(dimensions, mesh, dt);
    IB flag(flag_mesh);

    for(auto& node_ : L.get_nodes()){
	node_->calc_eq();
	node_->ftom(&L, 1);
	node_->update_f();
}

/* --- Main Loop --- */
for(double t=0.; t<100; t+=dt){

	cout << t << endl;

for(auto& node_ : L.get_nodes() ){ //  

    node_->ftom(&L, 0);

    node_->calc_eq();

force[0] = 0.; force[1] = 0.; force[2] = 0.; kernel = 0.; internal = false;

//----------------Needs to be converted to own function------------------------

// Only fluid nodes near the IB need to be looked at //
// Will be best to use Bounding Box //
if( (node_->get_pos()[2] > 0.35 - 2*mesh_space) && (node_->get_pos()[2] < 0.65 + 2*mesh_space) 
	&& (node_->get_pos()[1] > 0.35 - 2*mesh_space) && (node_->get_pos()[1] < 0.65 + 2*mesh_space) ){
for(auto& IB_node_ : flag.get_IB_nodes() ) {
    dist_x = (node_->get_pos()[0] - IB_node_->get_IB_pos()[0]) ;
    dist_y = (node_->get_pos()[1] - IB_node_->get_IB_pos()[1]) ;
    dist_z = (node_->get_pos()[2] - IB_node_->get_IB_pos()[2]) ;

// Internal nodes are treated as solid, Do not need IB consideration //
    if(pow( node_->get_pos()[1] - 0.5, 2 ) + pow( node_->get_pos()[2] - 0.5, 2 ) <= 0.15*0.15) internal = true;
    if(internal) continue;

    dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
    dist = sqrt(dist);
    if(dist <= 2*mesh_space ){
	kernel = (1 - abs( dist_x/mesh_space/2)) * (1 - abs( dist_y/mesh_space/2)) * (1 - abs( dist_z/mesh_space/2));
	for(int i=0; i<3; i++) force[i] = -1*node_->get_m()[i+1]/node_->get_m()[0] * kernel; 
	}
    }
}
//----------------------------------------------------------------------

    node_->collision(dt, force, internal);

    node_->ftom(&L, 1);
}

    L.stream();

for(auto& node_ : L.get_nodes()){
    node_->update_f();


/* --- Visualize --- */
if(node_->get_pos()[0] == 0.5){
//    for(int i=0; i<19; i++){
        cout << node_->get_m()[3]/node_->get_m()[0]*10000 << " " ;
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

