#include "src/Solver/Solver.hpp"

int main(int argc, char* argv[])
{
    if (argc < 2)
    {   
        std::cerr << "No window size given! " << std::endl;
        std::cerr << "Usage: " << argv[0] << " 800" << std::endl;
    };

    Solver(static_cast<uint>(std::stoi(argv[1]))).run();
    
    return 0;
}
