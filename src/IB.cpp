#include "header.hpp"
#include "IB.hpp"
//#include "IB_node.hpp"

IB::IB(){

int mesh_size[3] = {5, 5, 0};

for(int i=0; i<mesh_size[0]; i++){
        for(int j=0; j<mesh_size[1]; j++){
            for (int k=0; k<mesh_size[2]; k++){
                IB_node_lst_.push_back(std::make_shared<IB_node>(i*0.1+0.5,j*0.1+0.5,k*0.1+0.5, i*mesh_size[1]*mesh_size[2]+j*mesh_size[2]+k));
            }
        }
    }

}


