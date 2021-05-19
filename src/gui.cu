#include "interface.cuh"

int main(int *argcp, char**argv)
{

	interface myInt;
	myInt.InitWindow(argcp, argv);
	return 0;
}