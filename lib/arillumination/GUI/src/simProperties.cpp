#include "simProperties.h"

// constant get properties
uint64_t simProperties::get_nPhotons() const {return nPhotons;}
uint64_t simProperties::get_nPPerThread() const {return nPPerThread;}
uint64_t simProperties::get_nPhotonsTrue() const {return nPhotonsTrue;}
uint64_t simProperties::get_threadsPerBlock() const {return threadsPerBlock;}
uint64_t simProperties::get_nBlocks() const {return nBlocks;}

simProperties::simProperties()
{
	calc_threads();
}

void simProperties::set_nPhotons(const uint64_t _nPhotons)
{
	nPhotons = _nPhotons;
	calc_threads();
	return;
}

void simProperties::set_nPPerThread(const uint64_t _nPPerThread)
{
	nPPerThread = _nPPerThread;
	calc_threads();
	return;
}

void simProperties::calc_threads()
{
	nBlocks = nPhotons / (threadsPerBlock * nPPerThread);
	nPhotonsTrue = nBlocks * threadsPerBlock * nPPerThread;
	return;
}