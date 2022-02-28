#include "../src/model.h"

int main()
{
	model testMod;
	testMod.set_dim(200, 101, 101, 500);

	model newMod;
	newMod = testMod; // assignment operator called

	return 0;
}