
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
		//std::unique_ptr<node> node_ = std::make_unique<node>();
                //node_(i,j,k,0);
		node_lst_.push_back(std::make_unique<node>(i*mesh_x,j*mesh_y,k*mesh_z,0));
            }
        }
    }


    //node_lst_ = node_lst;
}


