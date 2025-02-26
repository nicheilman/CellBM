
#include "lattice.hpp"
//#include "node.hpp"

lattice::lattice(

	std::vector<float> domain_size,
        std::vector<int> mesh_size

)

{

    for(int i=0; i<domain_size[0]; i++){
	node_lst_.push_back(node(0,0,0,0));

}


    //node_lst_ = node_lst;
}


