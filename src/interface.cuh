/* 
	class representing GUI and handling user interfacing
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
 	Date: 24.08.2020	
*/

#ifndef INTERFACE_H
#define INTERFACE_H

#include "reconSettings.h"
#include "mbrecon.cuh"
#include "transducerProperties.h"
#include "transModel.cuh"

#include <SDL2/SDL.h>
#include <GL/glew.h>    // Initialize with gl3wInit()
#include "../lib/imgui/imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "ImGuiFileDialog.h"
#include "imgui_plot.h"
#include <thread>  
#include "color_mapper.cuh"
#include "../lib/vtkwriter/vtkwriter.h"

using namespace std;

class interface
{
	public:
		interface(); // class constructor
		~interface(); // class destructor
		void InitWindow(int *argcp, char**argv);	
		void MainDisplayCode();
		void ReconSettingWindow();
		void TransducerSettingsWindow();
		void Reconstructor();
		void ModelBuilder();
		void InputVolViz();
		void OutputVolViz();
		void ImImagesc(const float* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap);

	private:

		// variables used for plotting of input volume
		color_mapper myMapper;
		GLuint in_vol_texture;
		GLuint in_vol_texture_slice;
		int currSliceZ = 0;
		int currSliceY = 0;

		// variables used for plotting of output volume
		color_mapper reconVolMapper;
		GLuint rec_vol_texture;
		GLuint rec_vol_texture_slice;
		float currSliceZRec = 0; // specified in mm
		float currSliceYRec = 0; // specified in mm

		// variables used for plotting of model		
		color_mapper modelMapper;
		GLuint mod_texture;

		color_mapper modelCutMapper;
		GLuint modcut_texture;
		
		float zLvlModel = 0; // z slice to show of model matrix
		float rLvlModel = 0; // r slice to show of model matrix


		const char* windowTitle = "MB Recon GUI";
		ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 0.10f);

		bool show_recon_window = 0;
		bool show_trans_window = 1;
		bool show_data_loader = 1;
		bool show_recon_vol = 0;
		bool show_reconstructor = 1;
		
		bool isDataSetDefined = 0; 
		bool isTransducerValid = 0;
		bool isImpulseLoaded = 0; // is impulse response loaded
		bool isReconVolDefined = 0; // did we perform reconstruction
		bool isModelBuilt = 0; // did we build ourt transducer model

		uint8_t transducer_name_length = 64;
		uint8_t reconIter;
		
		// subclasses providing functionality
		mbrecon recon;
		reconSettings sett;
		transducerProperties* trans;
		transModel* model;
		modelMatrix* modMatrix;
		volume* inputDataVol;
		volume* reconVol;
		gaussian_beam* gauss;
		fieldProperties* fieldProps;
};

#endif
