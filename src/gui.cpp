#include "gui.h"


// displays a small help marker next to the text
static void HelpMarker(const char* desc)
{
	ImGui::TextDisabled("(?)");
	if (ImGui::IsItemHovered())
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
	return;
}

gui::gui()
{
	sett = rec.get_psett();
	mod = rec.get_pmodel();
	sigMat = rec.get_psigMat();
	absMat = rec.get_pabsMat();
	sensField = mod->get_psensField();
}

void gui::init(int *argcp, char**argv)
{
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
	{
	  printf("Error: %s\n", SDL_GetError());
		return;
	}
	// main_display_function goes somewhere here
	const char* glsl_version = "#version 140";
	// to find out which glsl version you are using, run glxinfo from terminal
	// and look for "OpenGL shading language version string"
	// https://en.wikipedia.org/wiki/OpenGL_Shading_Language

	SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
	
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

	SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
	SDL_Window* window = SDL_CreateWindow(windowTitle, 
		SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1900, 1080, window_flags);
	SDL_GLContext gl_context = SDL_GL_CreateContext(window);
	SDL_GL_MakeCurrent(window, gl_context);
	SDL_GL_SetSwapInterval(1); // Enable vsync

	bool err = glewInit() != GLEW_OK;
	if (err)
	{
		printf("Failed to initialize OpenGL loader!");
	  throw "FailedOpenGLLoader";
	}

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	ImGui::StyleColorsDark();
	ImGui_ImplSDL2_InitForOpenGL(window, gl_context);

	ImGui_ImplOpenGL3_Init(glsl_version);
	bool done = false;
	while (!done)
	{
		SDL_Event event;
		while (SDL_PollEvent(&event))
		{
			ImGui_ImplSDL2_ProcessEvent(&event);
			if (event.type == SDL_QUIT)
				done = true;
		}
		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplSDL2_NewFrame(window);
		ImGui::NewFrame();
		MainDisplayCode();
		// Rendering
		ImGui::Render();
		
		glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		//glUseProgram(0); // You may want this if using this code in an OpenGL 3+ context where shaders may be bound
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		SDL_GL_SwapWindow(window);
	}

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplSDL2_Shutdown();
	ImGui::DestroyContext();

	SDL_GL_DeleteContext(gl_context);
 	SDL_DestroyWindow(window);
 	SDL_Quit();
	return;
}

void gui::MainDisplayCode()
{
	SettingsWindow(); // settings for reconstruction
	ModelWindow(); // model loading and preview
	DataLoader(); // loading and previewing signal matrix
	ReconPreview(); // small preview window for reconstructed datasets
	SystemInformation();
	return;
}

// all the settings of the reconstruction
void gui::SettingsWindow()
{
	ImGui::Begin("Reconstruction");

	// here the user can define all important reconstruction settings
	if (ImGui::CollapsingHeader("Settings"))
	{
		float regStength = sett->get_regStrength();
		ImGui::InputFloat("Reg strength", &regStength);
		sett->set_regStrength(regStength);
		ImGui::SameLine();
		HelpMarker("Defines the regularization strength applied during inversion. Noise free datasets can run with 0 regularization, noisy datasets typically require a value around 5");

		int nIter = sett->get_nIter();
		ImGui::InputInt("Number of iterations", &nIter);
		sett->set_nIter(nIter);
		ImGui::SameLine();
		HelpMarker("The inversion runs iteratively. With each iteration of the LSQR, the result should converge a bit more towards the true solution. Typically more iterations will result in a more accuracte result but require more time to be calculated.");

		vector3<float> dr = sett->get_dr() * 1e6f;
		ImGui::InputFloat3("Resolution of output [microm]", &dr.x);
		sett->set_dr(dr * 1e-6);

		vector3<float> startPos = sett->get_startPos() * 1e3f;
		ImGui::InputFloat3("MinPos [mm]", &startPos.x);
		sett->set_startPos(startPos * 1e-3);

		vector3<float> endPos = sett->get_endPos() * 1e3f;
		ImGui::InputFloat3("MaxPos [mm]", &endPos.x);
		sett->set_endPos(endPos * 1e-3);

		float tLim[2] = {sett->get_tMin() * 1e6f, sett->get_tMax() * 1e6f}; 
		ImGui::InputFloat2("tLim [micros]", tLim);
		sett->set_tMin(tLim[0] * 1e-6f);
		sett->set_tMax(tLim[1] * 1e-6f);
	}

	if (ImGui::CollapsingHeader("Program status"))
	{

		ImGui::Columns(2);

		// check if the model was loaded
		ImGui::Text("Model loaded?"); 
		ImGui::SameLine();
		HelpMarker("To run a reconstruction, you first need to load a valid model from a file using the model matrix interface.");
		ImGui::NextColumn();
		boolIndicator(mod->get_isLoaded());
		ImGui::NextColumn();

		// check if the dataset was loaded from the harddrive
		ImGui::Text("Data loaded?"); 
		ImGui::SameLine();
		HelpMarker("To run a reconstruction, you first need to load a dataset from a h5 file using the data loader interface.");
		ImGui::NextColumn();
		boolIndicator(isSigMatLoaded);
		ImGui::NextColumn();

		// check if the GPU was detected properly
		ImGui::Text("GPU detected?");
		ImGui::SameLine();
		HelpMarker("The reconstruction runs on a GPU and we therefore need to find at least one CUDA capable device on your machine.");
		ImGui::NextColumn();
		boolIndicator(ginfo.get_nGpus() > 0);
		ImGui::NextColumn();

		ImGui::Text("Reconstruction running?");
		ImGui::SameLine();
		HelpMarker("This indicator will switch to true in the moment we start our reconstruction. The process will run in a separate thread so that the GUI can stay responsive.");
		ImGui::NextColumn();
		boolIndicator(rec.get_isRunning());
		ImGui::NextColumn();

		ImGui::Text("Abort request?");
		ImGui::SameLine();
		HelpMarker("If the user send a request to abort the execution of the current reconstruction, this will turn to true.");
		ImGui::NextColumn();
		boolIndicator(rec.get_flagAbort());
		ImGui::Columns(1);
	}

	if (rec.get_isRunning())
	{
		if (ImGui::CollapsingHeader("Reconstruction status"))
		{
			ImGui::Columns(2);
			// start time of reconstruction
			const time_t tt = system_clock::to_time_t(rec.get_tStart());
			ImGui::Text("Start time");
			ImGui::NextColumn(); 
			ImGui::Text("%s", ctime(&tt));
			ImGui::NextColumn();
			
			// elapsed time since start
			ImGui::Text("Elapsed time");
			ImGui::NextColumn();

			ImGui::NextColumn();

			// current iteration
			ImGui::Text("Iteration");
			ImGui::NextColumn();
			ImGui::Text("%d / %d", rec.get_iIter(), sett->get_nIter());
			ImGui::NextColumn();
			
			ImGui::Text("Convergence");
			ImGui::NextColumn();
			if (rec.get_iIter() > 1)
			{
				ImGui::Text("%.2f", rec.get_relConvergence());
			}
			else
			{
				ImGui::Text("-");
			}
			
			ImGui::NextColumn();
			ImGui::Text("Status"); 
			ImGui::NextColumn();
			ImGui::Text("%s", rec.get_statusVerbal());

			ImGui::Columns(1);
		}
		ImGui::Columns(2);
		
	}
	else
	{
		if (isReconRunning) 
		{
			// if we get here the reconstruction just finished
			reconThread.join();
			isReconRunning = 0;
		}
	}

	ImGui::Columns(2);
	ReconButton();
	ImGui::NextColumn();

	// add a button to stop a runing reconstruction
	AbortButton();
	ImGui::NextColumn();

	ImGui::End();
	return;
}

// a button to start a reconstruction which is either active or inactive
inline void gui::ReconButton()
{
	// add a button to start a reconstruction
	const bool flagCanRecon = canRecon();
	if (!flagCanRecon)
	{
		ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
	}

	if (ImGui::Button("Reconstruct"))
	{
		isReconRunning = 1;
		recon* recPtr = &rec;
		reconThread = recPtr->reconstruct2thread();
	}

	if (!flagCanRecon)
	{
		ImGui::PopItemFlag();
  	ImGui::PopStyleVar();
	}
	return;
}

// an abort button to stop reconstruction after the next iteration (clean exit)
inline void gui::AbortButton()
{
	const bool enableAbortButton = rec.get_isRunning() && (!rec.get_flagAbort());
	if (!enableAbortButton)
	{
		ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
	}

	if (ImGui::Button("Terminate"))
	{
		rec.set_flagAbort(1);	
	}

	if (!enableAbortButton)
	{
		ImGui::PopItemFlag();
  	ImGui::PopStyleVar();
	}

	ImGui::SameLine();
	HelpMarker("Abort the currently running reconstruction after finishing the next iteration through the LSQR. Could take a bit of time after requesting abort.");
	return;
}

// some infos about the workstation
void gui::SystemInformation()
{
	ImGui::Begin("System information");
	if (ImGui::CollapsingHeader("GPU"))
	{
		ImGui::Columns(2);
		for (int iGpu = 0; iGpu < ginfo.get_nGpus(); iGpu++)
		{
			const cudaDeviceProp* currProps = ginfo.get_props(iGpu);
			
			// device name
			ImGui::Text("Device name");
			ImGui::NextColumn();
			ImGui::Text("%s", currProps->name);
			ImGui::NextColumn();

			ImGui::Text("Total global memory [Gb]");
			ImGui::NextColumn();
			ImGui::Text("%lu", currProps->totalGlobalMem / (1024 * 1024 * 1024));
			ImGui::NextColumn();

			ImGui::Text("Clock rate [GHz]");
			ImGui::NextColumn();
			ImGui::Text("%.2f", (float) currProps->clockRate / (1024.0f * 1024.0f));
			ImGui::NextColumn();

			ImGui::Text("Number of multiprocessors [1]");
			ImGui::NextColumn();
			ImGui::Text("%d", currProps->multiProcessorCount);
			ImGui::NextColumn();

			ImGui::Text("Warp size [1]");
			ImGui::NextColumn();
			ImGui::Text("%d", currProps->warpSize);
			ImGui::NextColumn();

			ImGui::Text("Memory clock rate [GHz]");
			ImGui::NextColumn();
			ImGui::Text("%.2f", ((float) currProps->memoryClockRate) / (1024.0f * 1024.0f));
			ImGui::NextColumn();
		
			ImGui::Text("Memory Bus Width [bits]");
			ImGui::NextColumn();
			ImGui::Text("%d", currProps->memoryBusWidth);
			ImGui::NextColumn();


			const int pmb = 2.0*currProps->memoryClockRate*(currProps->memoryBusWidth/8)/1.0e6;
			ImGui::Text("Peak Memory Bandwidth [GB/s]");
			ImGui::NextColumn();
			ImGui::Text("%d", pmb);
			ImGui::NextColumn();


		}
		ImGui::Columns(1);
	}

	if (ImGui::CollapsingHeader("Host"))
	{

	}
	ImGui::End();
	return;
}

// button to load the model matrix from a file
inline void gui::LoadModelButton()
{
	if (ImGui::Button("Load from file"))
	{
		ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", 
			"Choose model matrix file", ".h5\0", ".");
	}

	if (ImGuiFileDialog::Instance()->FileDialog("ChooseFileDlgKey")) 
	{
		if (ImGuiFileDialog::Instance()->IsOk == true)
		{
			std::string inputFilePath = ImGuiFileDialog::Instance()->GetFilepathName();
			mod->load_from_file(inputFilePath);
			mod->calc_sensField();
			mod->calc_minmax();

			if (mod->get_tEnd() < sett->get_tMax())
			{
				sett->set_tMax(mod->get_tEnd());
			}

			if (mod->get_t0() > sett->get_tMin())
			{
				sett->set_tMin(mod->get_t0());
			}
		}
		ImGuiFileDialog::Instance()->CloseDialog("ChooseFileDlgKey");
	}
	return;
}

// make a button to reload model from the same file once we did the first loading
inline void gui::ReloadModelButton()
{
	if (!mod->get_isLoaded())
	{
		ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
	}

	if (ImGui::Button("Reload"))
	{
		mod->load_from_file();
		mod->calc_sensField();
		mod->calc_minmax();

		if (mod->get_tEnd() < sett->get_tMax())
		{
			sett->set_tMax(mod->get_tEnd());
		}

		if (mod->get_t0() > sett->get_tMin())
		{
			sett->set_tMin(mod->get_t0());
		}
	}

	if (!mod->get_isLoaded())
	{
		ImGui::PopItemFlag();
  	ImGui::PopStyleVar();
	}

	return;
}

// interfacing with the model matrix
void gui::ModelWindow()
{
	ImGui::Begin("Model matrix");
	ImGui::Columns(2);

	LoadModelButton();
	ImGui::NextColumn();

	ReloadModelButton();
	ImGui::NextColumn();
	ImGui::Columns(1);

	
	if (mod->get_isLoaded())
	{
		// some general informations about model size
		if (ImGui::CollapsingHeader("Model information"))
		{
			ImGui::Columns(2);
			ImGui::Text("File path");
			ImGui::NextColumn();
			ImGui::Text("%s", mod->get_filePath());
			ImGui::NextColumn();

			// display resolution of dataset
			ImGui::Text("Resolution [x, y, z, t]");
			ImGui::NextColumn();
			ImGui::Text("%.1f, %.1f, %.1f, %.1f", 
				mod->get_dx() * 1e6, mod->get_dy() * 1e6, 
				mod->get_dz() * 1e6, mod->get_dt() * 1e9);
			ImGui::NextColumn();

			// display origin of dataset
			ImGui::Text("Origin [x, y, z, t]");
			ImGui::NextColumn();
			ImGui::Text("%.1f, %.1f, %.1f, %.1f", 
				mod->get_x0() * 1e3, mod->get_y0() * 1e3, 
				mod->get_z0() * 1e3, mod->get_t0() * 1e6);
			ImGui::NextColumn();

			ImGui::Text("End [x, y, z, t]");
			ImGui::NextColumn();
			ImGui::Text("%.1f, %.1f, %.1f, %.1f", 
				mod->get_xEnd() * 1e3, mod->get_yEnd() * 1e3, 
				mod->get_zEnd() * 1e3, mod->get_tEnd() * 1e6);
			ImGui::NextColumn();
			
			ImGui::Text("Size [x, y, z, t]");
			ImGui::NextColumn();
			ImGui::Text("%lu, %lu, %lu, %lu", 
				mod->get_nx(), mod->get_ny(), mod->get_nz(), mod->get_nt());		
			ImGui::NextColumn();

			ImGui::Text("Data range");
			ImGui::NextColumn();
			ImGui::Text("%.2f ... %.2f", mod->get_minVal(), mod->get_maxVal());
			ImGui::NextColumn();

			ImGui::Columns(1);
		}

		// model slicing 
		if (ImGui::CollapsingHeader("Model Preview"))
		{
			ImGui::SliderInt("tSlice", &modSliceZ, 0, sensField->get_dim(2) - 1);
			ImImagesc(sensField->get_psliceZ((uint64_t) modSliceZ),
				sensField->get_dim(0), sensField->get_dim(1), &modDataTexture, modDataMapper);
			
			int width = 400;
			int height = (float) width / sensField->get_length(0) * sensField->get_length(1);
			ImGui::Image((void*)(intptr_t)modDataTexture, ImVec2(width, height));
			
			ImGui::SliderInt("xSlice", &modSliceY, 0, sensField->get_dim(0) - 1);

			ImImagesc(sensField->get_psliceY((uint64_t) modSliceY),
			 	sensField->get_dim(0), sensField->get_dim(2), &modDataTextureSlice, modDataMapper);
			height = (float) width / sensField->get_length(0) * sensField->get_length(2);
			ImGui::Image((void*)(intptr_t)modDataTextureSlice, ImVec2(width, height));
			
			ImGui::SliderFloat("MinVal", modDataMapper.get_pminVal(), sensField->get_minVal(), sensField->get_maxVal(), "%.1f");
			ImGui::SliderFloat("MaxVal", modDataMapper.get_pmaxVal(), sensField->get_minVal(), sensField->get_maxVal(), "%.1f");
			ImGui::ColorEdit4("Min color", modDataMapper.get_pminCol(), ImGuiColorEditFlags_Float);
			ImGui::ColorEdit4("Max color", modDataMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
		}
	}



	ImGui::End();
	return;
}

// small helper to load and preview raw data
void gui::DataLoader()
{
	ImGui::Begin("Data loader");
	ImGui::Columns(2);

	LoadDataButton(); // button to load a new dataset from a file
	ImGui::NextColumn();
	
	ReloadDataButton(); // button to reload dataset from disc
	ImGui::NextColumn();

	ImGui::Columns(1);

	if (isSigMatLoaded)
	{
		if (ImGui::CollapsingHeader("Dataset information"))
		{
			// print information about dataset if defined
			ImGui::Columns(2);	
			ImGui::Text("File"); ImGui::NextColumn();
			string filePath = sigMat->get_filePath();
			ImGui::Text("%s", filePath.c_str());
			ImGui::NextColumn();

			ImGui::Text("Size [x, y, t]"); 
			ImGui::NextColumn();
			ImGui::Text("%lu x %lu x %lu", sigMat->get_dim(0), sigMat->get_dim(1), sigMat->get_dim(2)); 
			ImGui::NextColumn();

			ImGui::Text("Memory [Mb]"); 
			ImGui::NextColumn();
			ImGui::Text("%f", (float) sigMat->get_nElements() * sizeof(float) / (1024.0f * 1024.0f));
			ImGui::NextColumn();

			ImGui::Text("Resolution [microm, x, y]"); ImGui::NextColumn();
			ImGui::Text("%.2f x %.2f", 
				sigMat->get_res(0) * 1e6f, 
				sigMat->get_res(1) * 1e6f);
			ImGui::NextColumn();
			
			ImGui::Text("Sampling freqeucny [MHz]"); ImGui::NextColumn();
			ImGui::Text("%.2f",	1.0f / (sigMat->get_res(2) * 1e6f));
			ImGui::NextColumn();

			ImGui::Text("t range [micros]"); ImGui::NextColumn();
			ImGui::Text("%.2f ... %.2f", 
				sigMat->get_minPos(2) * 1e6f, 
				sigMat->get_maxPos(2) * 1e6f);
			ImGui::NextColumn();

			ImGui::Text("x range [mm]"); ImGui::NextColumn();
			ImGui::Text("%.2f ... %.2f", 
				sigMat->get_minPos(0) * 1e3f, 
				sigMat->get_maxPos(0) * 1e3f);
			ImGui::NextColumn();

			ImGui::Text("y range [mm]"); ImGui::NextColumn();
			ImGui::Text("%.2f ... %.2f", 
				sigMat->get_minPos(1) * 1e3f, 
				sigMat->get_maxPos(1) * 1e3f);
			ImGui::NextColumn();

			ImGui::Text("Value range"); 
			ImGui::NextColumn();
			ImGui::Text("%.2f ... %.2f", sigMat->get_minVal(), sigMat->get_maxVal());
			ImGui::NextColumn();

			ImGui::Columns(1);
		}
		
		if (ImGui::CollapsingHeader("Dataset preview"))
		{
			ImGui::SliderInt("tSlice", &currSliceZ, 0, sigMat->get_dim(2) - 1);
			
			ImGui::Columns(2);
			float tSlice = sigMat->get_pos(currSliceZ, 2);
			ImGui::Text("tLayer: %.2f micros", tSlice * 1e6);
			ImGui::NextColumn();
			ImGui::Text("zLayer: %.2f mm", tSlice * sett->get_sos() * 1e3);
			ImGui::NextColumn();

			ImGui::Columns(1);
			ImImagesc(sigMat->get_psliceZ((uint64_t) currSliceZ),
				sigMat->get_dim(0), sigMat->get_dim(1), &inDataTexture, inDataMapper);
			
			int width = 550;
			int height = (float) width / sigMat->get_length(0) * sigMat->get_length(1);
			ImGui::Image((void*)(intptr_t)inDataTexture, ImVec2(width, height));
			
			ImGui::SliderInt("xSlice", &currSliceY, 0, sigMat->get_dim(1) - 1);

			ImImagesc(sigMat->get_psliceY((uint64_t) currSliceY),
			 	sigMat->get_dim(0), sigMat->get_dim(2), &inDataTextureSlice, inDataMapper);
			height = (float) width / sigMat->get_length(0) * sigMat->get_length(2) * sett->get_sos();
			ImGui::Image((void*)(intptr_t)inDataTextureSlice, ImVec2(width, height));
			
			ImGui::SliderFloat("MinVal", inDataMapper.get_pminVal(), sigMat->get_minVal(), sigMat->get_maxVal(), "%.1f");
			ImGui::SliderFloat("MaxVal", inDataMapper.get_pmaxVal(), sigMat->get_minVal(), sigMat->get_maxVal(), "%.1f");
			ImGui::ColorEdit4("Min color", inDataMapper.get_pminCol(), ImGuiColorEditFlags_Float);
			ImGui::ColorEdit4("Max color", inDataMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
		}
	}

	ImGui::End();
	return;
}

// a small button to reload the previously loaded dataset from file
inline void gui::ReloadDataButton()
{
	if (!isSigMatLoaded)
	{
		ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
	}

	if (ImGui::Button("Reload"))
	{
		sigMat->readFromFile();
		sigMat->calcMinMax();

			if (sigMat->get_maxPos(2) < sett->get_tMax())
			{
				sett->set_tMax(sigMat->get_maxPos(2));
			}

			if (sigMat->get_minPos(2) > sett->get_tMin())
			{
				sett->set_tMin(sigMat->get_minPos(2));
			}
	}

	if (!isSigMatLoaded)
	{
		ImGui::PopItemFlag();
  	ImGui::PopStyleVar();
	}
	ImGui::SameLine();
	HelpMarker("Will reload the file from the same location as initially specified");
	return;
}

inline void gui::LoadDataButton()
{
	if (ImGui::Button("Load input data"))
	{
		ImGuiFileDialog::Instance()->OpenDialog("fileDialogDataLoader", 
			"Choose input dataset", ".h5\0", ".");
	}

	if (ImGuiFileDialog::Instance()->FileDialog("fileDialogDataLoader")) 
	{
		if (ImGuiFileDialog::Instance()->IsOk == true)
		{
			std::string inputFilePath = ImGuiFileDialog::Instance()->GetFilepathName();
			sigMat->readFromFile(inputFilePath);
			sigMat->calcMinMax();

			if (sigMat->get_maxPos(2) < sett->get_tMax())
			{
				sett->set_tMax(sigMat->get_maxPos(2));
			}

			if (sigMat->get_minPos(2) > sett->get_tMin())
			{
				sett->set_tMin(sigMat->get_minPos(2));
			}

			isSigMatLoaded = 1;
		}
		ImGuiFileDialog::Instance()->CloseDialog("fileDialogDataLoader");
	}
	return;
}

// helper function to display stuff
void gui::ImImagesc(
	const float* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap)
{
	
	glDeleteTextures(1, out_texture);

	// Create an OpenGL texture identifier
	GLuint image_texture;
	glGenTextures(1, &image_texture);
	glBindTexture(GL_TEXTURE_2D, image_texture);

	// setup filtering parameters for display
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			
	// use color transfer function to convert from float to rgba
	unsigned char* data_conv = new unsigned char[4 * sizex * sizey];
	myCMap.convert_to_map(data, sizex * sizey, data_conv);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, sizex, sizey, 0, GL_RGBA, GL_UNSIGNED_BYTE, data_conv);

	// give pointer back to main program
	*out_texture = image_texture;
	delete[] data_conv; // free memory for temporary array
	return;
}

void gui::ReconPreview()
{
	if (rec.get_iIter() > 1)
	{
		absDataMapper.set_minVal(absMat->get_minVal());
		absDataMapper.set_maxVal(absMat->get_maxVal());
		ImGui::Begin("Reconstruction result");

		// preview of resulting dataset
		if (ImGui::CollapsingHeader("Slicer"))
		{
			ImGui::SliderInt("tSlice", &absSliceZ, 0, absMat->get_dim(2) - 1);
			const float zSlice = absMat->get_pos(absSliceZ, 2);
			ImGui::Text("Current z layer: %f mm", zSlice * 1e3f);
			ImImagesc(absMat->get_psliceZ((uint64_t) absSliceZ),
				absMat->get_dim(0), absMat->get_dim(1), &absDataTexture, absDataMapper);
			
			int width = 550;
			int height = (float) width / absMat->get_length(0) * absMat->get_length(1);
			ImGui::Image((void*)(intptr_t)absDataTexture, ImVec2(width, height));
			
			ImGui::SliderInt("xSlice", &absSliceY, 0, absMat->get_dim(0) - 1);

			ImImagesc(absMat->get_psliceY((uint64_t) absSliceY),
			 	absMat->get_dim(0), absMat->get_dim(2), &absDataTextureSlice, absDataMapper);
			height = (float) width / absMat->get_length(0) * absMat->get_length(2);
			ImGui::Image((void*)(intptr_t)absDataTextureSlice, ImVec2(width, height));
			
			ImGui::SliderFloat("MinVal", absDataMapper.get_pminVal(), 
				absMat->get_minVal(), absMat->get_maxVal(), "%.1f");
			ImGui::SliderFloat("MaxVal", absDataMapper.get_pmaxVal(), 
				absMat->get_minVal(), absMat->get_maxVal(), "%.1f");
			ImGui::ColorEdit4("Min color", absDataMapper.get_pminCol(), 
				ImGuiColorEditFlags_Float);
			ImGui::ColorEdit4("Max color", absDataMapper.get_pmaxCol(), 
				ImGuiColorEditFlags_Float);
		}

		if (ImGui::CollapsingHeader("Export"))
		{

		}

		ImGui::End();
	}
	return;
}

bool gui::canRecon()
{
	bool canRecon = 
		isSigMatLoaded && // data matrix needs to be loaded into file
		mod->get_isLoaded() && // model matrix needs to be loaded
		!isReconRunning && // there should not be any reconstruction running
		(ginfo.get_nGpus() > 0); // we need at least one valid GPU
	return canRecon;
}

// indicator if a bool is true or not. for now just a y n decision
void gui::boolIndicator(const bool status)
{
	if (status)
		ImGui::Text("y");
	else
		ImGui::Text("n");

	return;
}