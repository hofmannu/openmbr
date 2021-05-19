#ifndef INTERFACE_H
#define INTERFACE_H

#include <SDL2/SDL.h>
#include <GL/glew.h>    // Initialize with gl3wInit()
#include "../lib/imgui/imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "ImGuiFileDialog.h"
#include "imgui_plot.h"
#include "mc.cuh"
#include "fiberProperties.h"
#include "fieldProperties.h"
#include "optProperties.h"
#include "simProperties.h"
#include <cstdio>
#include "color_mapper.cuh"

using namespace std;

class interface
{
public:
	void InitWindow(int *argcp, char**argv);
	interface();
	~interface();

private:
	void MainDisplayCode();
	void Properties();
	void Result();
	void TissueProperties();
	void ImImagesc(
	const float* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap);


	const char* windowTitle = "Fluence GUI";
	// ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 0.10f);
	ImVec4 clear_color = ImVec4(0.60f, 0.55f, 0.45f, 0.10f);


	bool show_properties_window = 1;
	bool show_tissue_properties = 1;
	bool show_results = 1;

	bool is_output_defined = 0;

	mc sim;
	fiberProperties* arfiber;
	fieldProperties* field;
	simProperties* simprop;
	optProperties* tissue;
	optProperties* water;

	// elements required for plotting of resulting fluence
	color_mapper fluence_mapper;
	GLuint fluence_texture;

	bool flagLogPlot = 1;
	uint64_t lastNr;
	uint64_t lastNz;

};

#endif