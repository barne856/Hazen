#include <vector> // std::vector
#include <iostream> // std::cout

int main()
{
    std::vector<double> vec{6.0, 2.0, 3.0, 4.0, 5.0};
    for(int i = 0; i < vec.size(); i++)
    {
        std::cout << vec.size() << std::endl;
    }
return 0;
}