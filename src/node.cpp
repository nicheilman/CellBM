
#include "node.hpp"

node::node(const double dx, const double dy, const double dz, int wall_flag){

f_ = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
m_ = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
dx_ = dx;
dy_ = dy;
dz_ = dz;

}

//--------------------------------------------------------------------------------//

std::vector<double> node::generate_moment(){

return {0};

}

