
#include "header.hpp"


class node{

    protected:

    std::vector<float> f_, m_;

    public:

    node() = default;                         //default constructor
    node(const node& c) = default;           //copy constructor
    node(node&& c) = default;                 //move constructor
    node& operator=(const node& c) = default; //copy assignment operator
    node& operator=(node&& c) = default;      //move assignment operator 


    node( std::vector<int> node_lst );

    auto get_f(){return f_;};
    auto get_m(){return m_;};

};



