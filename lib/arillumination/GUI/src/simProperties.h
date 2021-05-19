#ifndef SIMPROPERTIES_H
#define SIMPROPERTIES_H

#include <cstdio>
#include <cstdint>

class simProperties
{
private:
	uint64_t nPhotons = 2e6; // wanted number of photons
	uint64_t nPhotonsTrue; // true number of simulated photons
	uint64_t nPPerThread = 10; // number of photons simulated in each thread
	uint64_t threadsPerBlock = 512; // threads per block on GPU
	uint64_t nBlocks; // dependent variable

	void calc_threads();
public:
	simProperties();

	void set_nPhotons(const uint64_t _nPhotons);
	void set_nPPerThread(const uint64_t _nPPerThread);

	uint64_t get_nPhotons() const;
	uint64_t get_nPPerThread() const;
	uint64_t get_nPhotonsTrue() const;
	uint64_t get_threadsPerBlock() const;
	uint64_t get_nBlocks() const;

};
#endif