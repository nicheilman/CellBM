
#include "lattice.hpp"

node::node(const double dx, const double dy, const double dz, int wall_flag){

f_ = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
m_ = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
dx_ = dx;
dy_ = dy;
dz_ = dz;

}

//--------------------------------------------------------------------------------//

void node::generate_moment(lattice* L){

std::fill(std::begin(m_), std::end(m_), 0);

for(int i=0; i<19; i++){
	for(int j=0; j<19; j++){
		m_[i] += L->evector[i][j] * f_[j];
	}
}

return;

}

