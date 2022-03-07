
#ifndef GUI_H
#define GUI_H

#include <SDL2/SDL.h>
#include <GL/glew.h>    // Initialize with gl3wInit()

#include <thread>
#include <ctime>
#include <chrono>
#include <cuda_runtime.h>
#include <cuda.h>


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
#include "gpu_info.h"


class gui
{
private:

	const char* windowTitle = "openmbr";
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 0.10f); // bg color

	recon rec;
	reconsett* sett; // pointer to reconstruction settings
	model* mod; // model matrix used for iterative inversion
	volume* absMat; // absorbance matrix
	volume* sensField; // sensitivity field of transducer (for plotting)

	gpu_info ginfo;

	volume* sigMat;
	bool isSigMatLoaded = 0; // is the signal matrix loaded and ready?
	bool isReconRunning = 0; // is a reconstruction currently in progress?

	void MainDisplayCode();
	void SettingsWindow();
	void ModelWindow();
	void DataLoader();
	void ReconPreview();
	void SystemInformation();

	inline void AbortButton();
	inline void ReconButton();

	inline void LoadDataButton();
	inline void ReloadDataButton();

	inline void LoadModelButton();
	inline void ReloadModelButton();

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

	// for model vizualization
	int absSliceZ = 0;
	int absSliceY = 0;
	color_mapper absDataMapper;
	GLuint absDataTexture;
	GLuint absDataTextureSlice;

	std::thread reconThread; // separate thread in which we will run the reconstruction
	bool canRecon();
	void boolIndicator(const bool status);

public:
	gui();
	void init(int *argcp, char**argv);

};

#endif