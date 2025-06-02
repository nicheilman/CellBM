
#include "lattice.hpp"
#include "IB.hpp"

node::node(const double dx, const double dy, const double dz, int wall_flag, int idx){

f_ = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
m_ = {1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
f_tmp = f_;
lambda = {0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
dx_ = dx;
dy_ = dy;
dz_ = dz;
wallflag_ = wall_flag;
idx_ = idx;

}

//--------------------------------------------------------------------------------//

void node::ftom(lattice* L, bool flag){
if(!flag){
std::fill(std::begin(m_), std::end(m_), 0);

//#pragma omp parallel for schedule(static)
for(int i=0; i<19; i++){
	for(int j=0; j<19; j++){
		m_[i] += L->evector[i][j] * f_[j];
	}
    }
}

else{
std::fill(std::begin(f_), std::end(f_), 0); 
//#pragma omp parallel for
for(int i=0; i<19; i++){
        for(int j=0; j<19; j++){
                f_[i] +=   (L->inv_evector[i][j] * m_[j]) / 144;
if(f_[i] < 0.0) f_[i] = 0.0;
        }
    }
std::fill(std::begin(f_eq), std::end(f_eq), 0);
//#pragma omp parallel for
for(int i=0; i<19; i++){
        for(int j=0; j<19; j++){
                f_eq[i] +=   (L->inv_evector[i][j] * meq[j]) / 144;
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

void node::collision(double dt, double fext[3], bool internal){
//#pragma omp parallel for
for(int i=4; i<velo_dim; i++)
	m_[i] += -1.0 * (m_[i] - meq[i]) * dt * lambda[i]*20 ;

if(internal) for(int i=0; i<3; i++) m_[i+1] = 0.0; 

double f_extx = fext[0];
double f_exty = fext[1];
double f_extz = fext[2];

if(dz_ == 0.0) f_extz = 0.0002;

if(!internal){
      m_[1] += f_extx;
      m_[2] += f_exty;
      m_[3] += f_extz;
      m_[4] += 2.0*(m_[1]*f_extx+m_[2]*f_exty+m_[3]*f_extz)/meq[0];
      m_[5] += (4.0*m_[1]*f_extx-2.0*(m_[2]*f_exty+m_[3]*f_extz))/meq[0];
      m_[6] += 2.0*(m_[2]*f_exty-m_[3]*f_extz)/meq[0];
      m_[7] += (m_[1]*f_exty+m_[2]*f_extx)/meq[0];
      m_[8] += (m_[3]*f_exty+m_[2]*f_extz)/meq[0];
      m_[9] += (m_[3]*f_extx+m_[1]*f_extz)/meq[0];
}

return;

}






