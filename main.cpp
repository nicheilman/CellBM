#include "header.hpp"
#include "lattice.hpp"

int main(int argc, char** argv){

    const std::vector<int> node_lst = {0,0,0};
    lattice L(node_lst);

    std::cout << L.get_node(1) << std::endl;

    return 0;

}


