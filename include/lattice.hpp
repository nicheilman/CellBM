
#include "header.hpp"

class lattice{

    protected:

    std::vector<int> node_lst_;

    public:

    lattice() = delete;                         //default constructor
    lattice(const lattice& c) = default;           //copy constructor
    lattice(lattice&& c) = delete;                 //move constructor
    lattice& operator=(const lattice& c) = delete; //copy assignment operator
    lattice& operator=(lattice&& c) = delete;      //move assignment operator 


    lattice( std::vector<float> dimensions, 
	std::vector<int> mesh_size );

    auto get_node(int node_id){return node_lst_[node_id];};


};

