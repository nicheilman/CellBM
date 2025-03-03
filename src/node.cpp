
#include "lattice.hpp"

node::node(const double dx, const double dy, const double dz, int wall_flag){

f_ = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
m_ = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
dx_ = dx;
dy_ = dy;
dz_ = dz;

}

//--------------------------------------------------------------------------------//

double node::generate_moment(lattice* L){

auto test = L->evector;


return test[5][1];

}

