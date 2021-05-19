/*
	class used to build simulations based on our forward model building
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 28.09.2020

	Todo:
		- implement different shapes of absorbers (started with struct definition)
		- implement fluence simulation fields and simulations
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "transducerProperties.h"
#include "transModel.cuh"
#include "volume.h"
#include "mbrecon.cuh"
#include "reconSettings.h"
#include "noise.h"

#include <SDL2/SDL.h>
#include <GL/glew.h>    // Initialize with gl3wInit()
#include "../lib/imgui/imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "ImGuiFileDialog.h"
#include "imgui_plot.h"
#include <thread>  
#include "color_mapper.cuh"
#include "../lib/arillumination/GUI/src/mc.cuh"
#include "../lib/arillumination/GUI/src/optProperties.h"
#include "../lib/arillumination/GUI/src/fiberProperties.h"
#include "../lib/arillumination/GUI/src/simProperties.h"
#include "../lib/arillumination/GUI/src/mcFieldProperties.h"


// point absorber
struct pointAbsorber
{
	float pos[3];
	float absorbance;
};

// line absorber
struct lineAbsorber
{
	float posStart[3];
	float posEnd[3];
	float absorbancePerLength;
};

// spherical absorber
struct sphericalAbsorber
{
	float centerPos[3];
	float radius;
	float absorbancePerVolume;
};

// tubular absorber
struct tubeAbsorber
{
	float posStart[3];
	float posEnd[3];
	float radius;
	float absorbancePerVolume;
};

class simulation
{
private:
	volume* absorberMatrix; // matrix representing absorber position
	volume* signalMatrix; // matrix representing signal
	transducerProperties* trans; // represents transducer properties
	transModel* model; // transducer model matrix builder
	modelMatrix* modMatrix; // transducer model matrix
	reconSettings* sett; // reconstruction settings
	mbrecon recon;
	noise noiser; // noise adder

	// class for fluence simulation
	mc mcsim;
	fiberProperties* arfiber; // multimode fiber properties
	simProperties* simprop; // monte carlo  simulation properties
	optProperties* tissueprop;
	mcFieldProperties* mcfieldprop;

	// to display transducer model
	color_mapper modelMapper;
	GLuint mod_texture;
	color_mapper modelCutMapper;
	GLuint modcut_texture;

	// to display input volume
	color_mapper myMapper;
	GLuint in_vol_texture;
	GLuint in_vol_texture_slice;

	int currSliceZ = 0;
	int currSliceY = 0;

	const char* windowTitle = "Simulation GUI";
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 0.10f);

	bool show_trans_window = 1;
	bool show_absorbance_window = 1;
	bool show_fluence_window = 1;

	bool isImpulseLoaded = 0; // did we already load transducer impulse response
	float zLvlModel = 0;
	float rLvlModel = 0;
	bool isModelBuilt = 0; 

	// fluence related settings
	bool flagFluenceModeling = 1;
	float zLevelSurfaceMm = 6; // z depth of simulated surface

	bool isSimulationDone = 0; // status flag if simulation finished
	bool isVolDefined = 0; // flag if volume is defined
	bool isAbsorberVolDefined = 0; // is the absorption volume defined?
	vector<pointAbsorber> pAbs; // vector representing point absorbers
	float fSampl = 250; // sampling frequency used for simulation in MHz
	float noiseLevel = 0.05; // noise level in factor of maximum amplitude

	void GenerateAbsorbanceVolume();
	void TransducerSettingsWindow();
	void AbsorbanceMatrixWindow();
	void SimulationSettings();
	void FluenceWindow();
	void PrintVolumeInformation(volume* vol);

	void ImImagesc(const float* data, const uint64_t sizex, const uint64_t sizey, 
		GLuint* out_texture, const color_mapper myCMap);


public:
	// class constructor and destructor
	simulation();
	~simulation();
	void InitWindow(int *argcp, char**argv);
	void MainDisplayCode();
	void Simulate();

};

#endif