
#include "lattice.hpp"
//#include "node.hpp"

lattice::lattice(

	double domain_size[3],
        int mesh_size[3]

)

{
    //node_lst_.reserve(domain_size[0]*domain_size[1]*domain_size[2]);

    double mesh_x = domain_size[0]/mesh_size[0];
    double mesh_y = domain_size[1]/mesh_size[1];
    double mesh_z = domain_size[2]/mesh_size[2];

    for(int i=0; i<mesh_size[0]; i++){
        for(int j=0; j<mesh_size[1]; j++){
            for (int k=0; k<mesh_size[2]; k++){
		if(i==0 || j==0 || k==0 || i==(mesh_size[0]-1) || j==(mesh_size[1]-1) || k==(mesh_size[2]-1) )
		    node_lst_.push_back(std::make_unique<node>(i*mesh_x,j*mesh_y,k*mesh_z,1));
		else
		    node_lst_.push_back(std::make_unique<node>(i*mesh_x,j*mesh_y,k*mesh_z,0));
            }
        }
    }


    //node_lst_ = node_lst;
}
//----------------------------------------------------------------------------------------//

void lattice::stream(){

for(auto& node_ : get_nodes() ){
return;
}

return;
}


