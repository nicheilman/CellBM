//#include "header.hpp"
#include "IB.hpp"
//#include "IB_node.hpp"

IB::IB(

    int mesh_size[3]

){

mesh_size_ = mesh_size;
for(int i=0; i<mesh_size_[0]+1; i++){
        for(int j=0; j<mesh_size_[1]+1; j++){
            for (int k=0; k<mesh_size_[2]+1; k++){
		IB_node_lst_.push_back(std::make_shared<IB_node>( i*(0.2/mesh_size[0]) +0.4, 0.25*sin(j*(6.2831/mesh_size[1]) )+0.5, 0.25*cos(j*(6.2831/mesh_size[1]) )+0.5 /*+ 0.5/mesh_size[0]*/, i*mesh_size_[1]+j));
            }
        }
    }

}


