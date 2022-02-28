/* 
	main file from which we start the reconstruction GUI
	Author: Urs Hofmann
	Mail: mail@hofmannu.org
	Date: 18.01.2022
*/

#include <iostream>
#include "gui.h"

using namespace std;

int main(int *argcp, char**argv)
{
	gui GUI;
	GUI.init(argcp, argv);

	return 0;
}