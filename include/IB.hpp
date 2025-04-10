
#include "header.hpp"
#include "IB_node.hpp"

class IB{

    protected:
    
    std::vector<std::shared_ptr<IB_node>> IB_node_lst_ = {};

    public:

    IB() = default;                         //default constructor
    IB(const IB& c) = default;           //copy constructor
    IB(IB&& c) = default;                 //move constructor
    IB& operator=(const IB& c) = default; //copy assignment operator
    IB& operator=(IB&& c) = default;      //move assignment operator 



};



