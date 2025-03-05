
#include "lattice.hpp"

node::node(const double dx, const double dy, const double dz, int wall_flag){

f_ = {1., 1., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
m_ = {1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
f_tmp = f_;
tau = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
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
//--------------------------------------------------------------------//

void node::calc_eq(){

            m2 = m_[1] * m_[1] + m_[2] * m_[2] + m_[3] * m_[3];
            meq[0] = m_[0];
            meq[1] = m_[1];
            meq[2] = m_[2];
            meq[3] = m_[3];
            meq[4] = m2 / meq[0];
            meq[5] = (3 * m_[1] * m_[1] - m2) / meq[0];
            meq[6] = (m_[2] * m_[2] - m_[3] * m_[3]) / meq[0];
            meq[7] = m_[1] * m_[2] / meq[0];
            meq[8] = m_[3] * m_[2] / meq[0];
            meq[9] = m_[3] * m_[1] / meq[0];
            for (int q = 10; q < velo_dim; ++q)
                meq[q] = 0.0;

return;
}
//--------------------------------------------------------------------//

void node::collision(){

for(int i=0; i<velo_dim; i++)
	m_[i] += -1.0 * (m_[i] - meq[i]) / tau[i];

return;
}



