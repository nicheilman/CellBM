
#include "header.hpp"

class lattice;

class node: public std::enable_shared_from_this<node>{

    protected:
    
    const int velo_dim = 19;
    std::array<double, 19> f_, f_tmp, m_, meq, tau; 
    //std::array<double, 19> ;
    double dx_, dy_, dz_; 
    double m2;
    int wallflag_;

    public:

    node() = default;                         //default constructor
    node(const node& c) = default;           //copy constructor
    node(node&& c) = default;                 //move constructor
    node& operator=(const node& c) = default; //copy assignment operator
    node& operator=(node&& c) = default;      //move assignment operator 


    node( double dx, double dy, double dz, int wallflag );

    auto get_f(){return f_;};
    auto get_m(){return m_;};
    std::vector<double> get_pos(){return {dx_, dy_, dz_};};
    auto get_wallflag(){return wallflag_;};

    void set_f(std::shared_ptr<node> node_, int i){f_tmp[i] = node_->get_f()[i]; return;};

    void ftom(lattice* L, bool flag);
    void calc_eq();
    void collision();

};



