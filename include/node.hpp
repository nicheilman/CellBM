
#include "header.hpp"

class lattice;
class IB;

class node: public std::enable_shared_from_this<node>{

    protected:
    
    const int velo_dim = 19;
    std::array<double, 19> f_, f_tmp, m_, meq, lambda, f_eq; 
    double dx_, dy_, dz_;
    double p_, u_[3]; 
    double m2;
    int wallflag_, idx_;

    static constexpr double wk[19] = {3., 18., 18., 18., 18., 18., 18., 36., 36., 36., 36., 36., 36., 36., 36., 36., 36., 36., 36.};

    public:

    node() = default;                         //default constructor
    node(const node& c) = default;           //copy constructor
    node(node&& c) = default;                 //move constructor
    node& operator=(const node& c) = default; //copy assignment operator
    node& operator=(node&& c) = default;      //move assignment operator 


    node( double dx, double dy, double dz, int wallflag, int idx);

    auto get_f(){return f_;};
    auto get_m(){return m_;};
    std::vector<double> get_pos(){return {dx_, dy_, dz_};};
    std::vector<double> get_velocity(){return {m_[1]/m_[0], m_[2]/m_[0], m_[3]/m_[0]};};
    auto get_wallflag(){return wallflag_;};
    int get_idx(){return idx_;};

    void set_f(std::shared_ptr<node> node_, int i){f_tmp[i] = node_->get_f()[i]; return;};
    void set_f_wall(int i, int j, double c_i){f_tmp[i] = f_[j] + (c_i / wk[i]); return;}; 
    void set_f_p(int i, int j, int c_i){f_tmp[i] = f_eq[i]; return;};
    void update_f(){for(int i=1; i<velo_dim; i++)f_[i] = f_tmp[i]; return;};

    void ftom(lattice* L, bool flag);
    void calc_eq();
    void collision(double dt, double fext[3], bool internal);

};



