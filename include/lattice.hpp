
#include "header.hpp"
#include "node.hpp"

class node;

class lattice{

    protected:

    std::vector<std::unique_ptr<node>> node_lst_ = {};
    static constexpr double evector[19][19] = { {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, // Mass Density, \rho
					       {0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1},  // X-Velocity Density, v_x
					       {0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1},  // Y-Velocity Density, v_y
					       {0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0},  // Z-Velocity Density, v_z
					      {-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},  // Kinetic Energy, e
					       {0, 2, 2,-1,-1,-1,-1,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1},  // P_xx
					       {0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0,-1,-1,-1,-1, 1, 1, 1, 1},  // P_ww = P_yy - P_zz
					       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1},  // P_xy
					       {0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},  // P_yz
					       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0},  // P_zx
					       {0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1},  // Ghost Moments
					       {0, 0, 0,-2, 2, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1},  //      |
					       {0, 0, 0, 0, 0,-2, 2, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0},  //      |
					       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1},  //      |
					       {0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0,-1, 1, 1,-1},  //      |
					       {0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0},  //      |
					       {1,-2,-2,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},  //      |
					       {0,-2,-2, 1, 1, 1, 1,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1},  //      V
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

