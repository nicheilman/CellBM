
#include "header.hpp"
#include "node.hpp"

class node;

class lattice{

    protected:

    std::vector<std::unique_ptr<node>> node_lst_ = {};
    static constexpr double evector[19][19] = { {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //
				    {0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1},  //
				    {0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1},  //
				    {0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0},  //
				   {-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},  //
				    {0, 2, 2,-1,-1,-1,-1,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1},  //
				    {0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0,-1,-1,-1,-1, 1, 1, 1, 1},  //
				    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1},  //
				    {0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},  //
				    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0},  //
				    {0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1},  //
				    {0, 0, 0,-2, 2, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1},  //
				    {0, 0, 0, 0, 0,-2, 2, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0},  //
				    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1},  //
				    {0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0,-1, 1, 1,-1},  //
				    {0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0},  //
				    {1,-2,-2,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},  //
				    {0,-2,-2, 1, 1, 1, 1,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1},  //
				    {0, 0, 0,-1,-1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1, 1, 1, 1, 1} };//

    friend class node;

    public:

    lattice() = delete;                         //default constructor
    lattice(const lattice& c) = default;           //copy constructor
    lattice(lattice&& c) = delete;                 //move constructor
    lattice& operator=(const lattice& c) = delete; //copy assignment operator
    lattice& operator=(lattice&& c) = delete;      //move assignment operator 


    lattice( double domain_size[], 
	int mesh_size[] );

    auto& get_nodes(){return node_lst_;};
    auto get_evector(int i, int j){return evector[i][j];};

};

