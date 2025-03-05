
#include "lattice.hpp"

node::node(const double dx, const double dy, const double dz, int wall_flag){

f_ = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
m_ = {1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
dx_ = dx;
dy_ = dy;
dz_ = dz;

}

//--------------------------------------------------------------------------------//

void node::ftom(lattice* L, bool flag){
if(!flag){
std::fill(std::begin(m_), std::end(m_), 0);

for(int i=0; i<19; i++){
	for(int j=0; j<19; j++){
		m_[i] += L->evector[i][j] * f_[j];
	}
    }
}

else{
std::fill(std::begin(f_), std::end(f_), 0); 

for(int i=0; i<19; i++){
        for(int j=0; j<19; j++){
                f_[i] +=   (L->inv_evector[i][j] * m_[j]) / 144;
        }
    }
}

return;
}

