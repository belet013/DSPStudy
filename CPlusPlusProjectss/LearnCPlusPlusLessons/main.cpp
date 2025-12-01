#include <iostream>

void doSomething(unsigned int x){

}

int main(){

    signed int s {-1};
    unsigned int u {1};

    if (s < u)
        std::cout << "-1 is less than 1 \n";
    else
        std::cout << "1 is less than -1 \n";
    return 0;
}