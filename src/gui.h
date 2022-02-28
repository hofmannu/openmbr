
#ifndef GUI_H
#define GUI_H

#include <SDL2/SDL.h>
#include <GL/glew.h>    // Initialize with gl3wInit()

#include <thread>
#include <ctime>
#include <chrono>

#include "../lib/imgui/imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "ImGuiFileDialog.h"
#include "imgui_plot.h"
#include "recon.h"
#include "vector3.h"
#include "reconsett.h"
#include "model.h"
#include "color_mapper.h"
#include "../lib/CVolume/src/volume.h"

class gui
{
private:

	const char* windowTitle = "openmbr";
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 0.10f); // bg color

	recon rec;
	reconsett* sett;
	model* mod;
	volume* absMat;
	volume* sensField;

	volume* sigMat;
	bool isSigMatLoaded = 0;

	void MainDisplayCode();
	void SettingsWindow();
	void ModelWindow();
	void DataLoader();

	void ImImagesc(const float* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap);

	// all those variables we use for the preview

	// for input dataset vizualization
	int currSliceZ = 0;
	int currSliceY = 0;
	color_mapper inDataMapper;
	GLuint inDataTexture;
	GLuint inDataTextureSlice;

	// for model vizualization
	int modSliceZ = 0;
	int modSliceY = 0;
	color_mapper modDataMapper;
	GLuint modDataTexture;
	GLuint modDataTextureSlice;

public:
	gui();

	void init(int *argcp, char**argv);

};

#endif