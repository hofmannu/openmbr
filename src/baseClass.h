#include <cstring>
#include <iostream>
#include <cstdio>

#ifndef BASECLASS_H
#define BASECLASS_H

using namespace std;

class baseClass
{
protected:
	bool flagVerbose = 1; // enable or disable verbose output
	void vprintf(const string txtMsg) const;
	void vprintf(const string txtMsg, const bool flagName) const;
	string className = "unkown";
};

#endif