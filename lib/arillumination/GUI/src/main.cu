#include "mc.cuh"
#include "interface.cuh"

int main(int *argcp, char**argv)
{
	interface GUI;
	GUI.InitWindow(argcp, argv);

	// float rSpot = fiber.get_rSpot(1.33, field.get_zBound());

	return 0;
}