
#include "header.hpp"
#include "node.hpp"

class lattice{

    protected:

    std::vector<std::unique_ptr<node>> node_lst_ = {};

    public:

    lattice() = delete;                         //default constructor
    lattice(const lattice& c) = default;           //copy constructor
    lattice(lattice&& c) = delete;                 //move constructor
    lattice& operator=(const lattice& c) = delete; //copy assignment operator
    lattice& operator=(lattice&& c) = delete;      //move assignment operator 


    lattice( std::vector<double> domain_size, 
	std::vector<int> mesh_size );

    auto& get_nodes(){return node_lst_;};


};

