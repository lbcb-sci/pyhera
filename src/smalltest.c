// #include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{

std::string s1 = "Ctg2_RC";
std::string s2 = "_RC";

std::cout << endl << "FIND: " << s1.rfind(s2) << endl;
std::cout << "SIZE1: " << s1.size() << endl;
std::cout << "SIZE2: " << s2.size() << endl;
}
