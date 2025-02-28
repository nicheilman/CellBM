
#include "header.hpp"


class node: public std::enable_shared_from_this<node>{

    protected:
    
    const int velo_dim = 19;
    std::array<double, 3> f_ = {1., 2., 3.};
    std::array<double, 19> m_;

    public:

    node() = default;                         //default constructor
    node(const node& c) = default;           //copy constructor
    node(node&& c) = default;                 //move constructor
    node& operator=(const node& c) = default; //copy assignment operator
    node& operator=(node&& c) = default;      //move assignment operator 


    node( double dx, double dy, double dz, bool wallflag );

    auto get_f(){return f_;};
    auto get_m(){return m_;};

};



