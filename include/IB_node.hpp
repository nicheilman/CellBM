
#include "header.hpp"


class IB_node: public std::enable_shared_from_this<IB_node>{

    protected:
    
    double dx_, dy_, dz_;
    double p_, u_[3]; 
    int idx_;


    public:

    IB_node() = default;                         //default constructor
    IB_node(const IB_node& c) = default;           //copy constructor
    IB_node(IB_node&& c) = default;                 //move constructor
    IB_node& operator=(const IB_node& c) = default; //copy assignment operator
    IB_node& operator=(IB_node&& c) = default;      //move assignment operator 


    IB_node( double dx, double dy, double dz, int idx);

    std::vector<double> get_pos(){return {dx_, dy_, dz_};};
    int get_idx(){return idx_;};


};



