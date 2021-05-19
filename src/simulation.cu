#include "simulation.cuh"

// class constructor
simulation::simulation()
{
	trans = recon.get_ptransProp(); // pull over transducer properties from reconstruction
	model = recon.get_ptransModel(); // pull over transducer model
	absorberMatrix = recon.get_preconVol();
	signalMatrix = recon.get_ppreprocVol();
	sett = recon.get_psett();

	// set default values for absorber volume
	for (uint8_t iDim = 0; iDim < 3; iDim++)
		absorberMatrix->setRes(iDim, 10e-6);
	
	absorberMatrix->setOrigin(0, 6e-3);
	absorberMatrix->setOrigin(1, -1e-3);
	absorberMatrix->setOrigin(2, -1e-3);

	absorberMatrix->setDim(0, 201);
	absorberMatrix->setDim(1, 201);
	absorberMatrix->setDim(2, 201);

	myMapper.set_mapType(1);

	// get pointers to fluence modeling
	arfiber = mcsim.get_pfiber();
	simprop = mcsim.get_psim();
	tissueprop = mcsim.get_ptissue();
	mcfieldprop = mcsim.get_pfield();
}

simulation::~simulation()
{


}

// initialize ImGui window
void simulation::InitWindow(int *argcp, char**argv)
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
	  throw "Failed to initialize OpenGL loader!";
		return;
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

void simulation::ImImagesc(
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

void simulation::TransducerSettingsWindow()
{
	ImGui::Begin("US Transducer", &show_trans_window);
	
	std::string name = trans->getName();
	ImGui::Text(name.c_str()); ImGui::NextColumn();
	
	static char arrayName[64];
	for (uint8_t iChar = 0; iChar < 64; iChar++)
	{
		if (iChar < name.length())
			arrayName[iChar] = name[iChar];
	}
	ImGui::InputText("Transducer name", arrayName, 64);
	name = arrayName;
	trans->setName(name);
	
	ImGui::Columns(2);
	
	// load and save buttons for transducer properties
	if (ImGui::Button("Load transducer")){
		ImGuiFileDialog::Instance()->OpenDialog("ReadTransDialog", 
			"Choose File", ".h5\0", ".");
	}
	if (ImGuiFileDialog::Instance()->FileDialog("ReadTransDialog")) 
	{
		if (ImGuiFileDialog::Instance()->IsOk == true)
		{
			std::string filePathName = ImGuiFileDialog::Instance()->GetFilepathName();
			trans->readFromFile(filePathName);
		}
		ImGuiFileDialog::Instance()->CloseDialog("ReadTransDialog");
	}
	ImGui::NextColumn();

	if (ImGui::Button("Save transducer")){
		ImGuiFileDialog::Instance()->OpenDialog("SaveTransDialog", 
			"Choose File", ".h5\0", ".");
	}
	if (ImGuiFileDialog::Instance()->FileDialog("SaveTransDialog")) 
	{
		if (ImGuiFileDialog::Instance()->IsOk == true)
		{
			std::string filePathName = ImGuiFileDialog::Instance()->GetFilepathName();
			trans->saveToFile(filePathName);
		}
		ImGuiFileDialog::Instance()->CloseDialog("SaveTransDialog");
	}
	
	// transducer property manipulation
	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Transducer properties"))
	{
		// define focal distance
		float focalDistanceMm = trans->getFocalDistance() * 1000;
		ImGui::InputFloat("fDist [mm]", &focalDistanceMm, 0.1f, 1.0f, "%0.3f");
		trans->setFocalDistance(focalDistanceMm / 1000);

		ImGui::SameLine(); HelpMarker("Focal distance of ultrasound transducer");
		float rApertureMm = trans->getRAperture() * 1000;
		ImGui::InputFloat("rAperture [mm]", &rApertureMm, 0.1f, 1.0f, "%0.3f");
		trans->setRAperture(rApertureMm / 1000);

		ImGui::SameLine(); HelpMarker("Aperture radius of transducer");
		float rHoleMm = trans->getRHole() * 1000;
		ImGui::InputFloat("rHole [m]", &rHoleMm, 0.1f, 1.0f, "%0.3f");
		trans->setRHole(rHoleMm / 1000);

		ImGui::SameLine(); HelpMarker("Hole radius of transducer");
		ImGui::InputFloat("fCentral [Hz]", trans->get_pcentralFrequency());
		ImGui::SameLine(); HelpMarker("Central freqeucny of transducer");
	}

	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Impulse signal"))
	{
		// todo: diable button if not found
		if (ImGui::Button("Load impulse signal"))
		{
			model->loadImpulseResponse();
			isImpulseLoaded = 1;
		}
		// todo: add plot of impulse response signal here
		
		if (isImpulseLoaded)
		{
			ImGui::PlotConfig conf;
			conf.values.xs = model->get_ptImp(); // this line is optional
			conf.values.ys = model->get_pimpResp();
			conf.values.count = model->getNElem();
			conf.scale.min = -1;
			conf.scale.max = 1;
			conf.tooltip.show = true;
			conf.tooltip.format = "x=%.2f, y=%.2f";
			conf.grid_x.show = false;
			conf.grid_y.show = false;
			conf.frame_size = ImVec2(450, 200);
			conf.line_thickness = 2.f;

			ImGui::Plot("plot", conf);
		}

		ImGui::Columns(2);
		if (isImpulseLoaded)
		{
			ImGui::Text("%.2f", model->get_tOffset() * 1e9);			
		}
		ImGui::NextColumn();
		ImGui::Text("Offset time [ns]"); ImGui::NextColumn();
		
		if (isImpulseLoaded)
		{
			float sigLength = 1 / model->get_fSampling() * 
				(float) model->getNElem();
			ImGui::Text("%.2f", sigLength * 1e9);
		}
		ImGui::NextColumn();
		ImGui::Text("Signal length [ns]"); ImGui::NextColumn();
		
		if (isImpulseLoaded)
		{
			ImGui::Text("%.2f", model->get_fSampling() * 1e-6);
		}
		ImGui::NextColumn();
		ImGui::Text("Sampling frequency [MHz]"); ImGui::NextColumn();
	}

	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Model building"))
	{
		int nElements = model->getNElements();
		ImGui::InputInt("Number of elements", &nElements);
		model->setNElements(nElements);


		ImGui::SameLine(); HelpMarker("Number of elements used for discretization");
		if (ImGui::Button("Build model"))
		{
			if (!isImpulseLoaded) // try to load impulse response if not done yet
				model->loadImpulseResponse();
			model->setTransProp(*trans);
			model->buildGPUModel();
			modMatrix = model->getModel();
			modMatrix->calcSensField();
			isModelBuilt = 1;
		}

		if (isModelBuilt)
		{
			// make two sliders for r and z position selection
			
			// plot sensitivity field
			ImImagesc(modMatrix->get_psensField(), modMatrix->getNrModel(), 
				modMatrix->getNzModel(), &mod_texture, modelMapper);
			
			int width = 250; 
			int height = 500;
			ImGui::Image((void*)(intptr_t)mod_texture, ImVec2(width, height));
			ImGui::SameLine();
			
			float* modSlice = modMatrix->get_pmodSlice(zLvlModel);
			modelCutMapper.set_maxVal(modMatrix->get_maxValSlice());
			modelCutMapper.set_minVal(modMatrix->get_minValSlice());
			ImImagesc(modSlice, modMatrix->getNrModel(), 
			 	modMatrix->getNtModel(), &modcut_texture, modelCutMapper);
			ImGui::Image((void*)(intptr_t)modcut_texture, ImVec2(width, height));

			float zDepthMm = zLvlModel * 1000;
			float rDepthMm = rLvlModel * 1000;
			ImGui::SliderFloat("z depth [mm]", &zDepthMm, modMatrix->get_zMin() * 1000, modMatrix->get_zMax() * 1000, "%.01f");
			ImGui::SliderFloat("r depth [mm]", &rDepthMm, 0.0, modMatrix->get_rMax() * 1000, "%.2f");
			zLvlModel = zDepthMm / 1000;
			rLvlModel = rDepthMm / 1000;

			ImGui::SliderFloat("MinVal", modelMapper.get_pminVal(), 0, modMatrix->getMaxVal(), "%.2f");
			ImGui::SliderFloat("MaxVal", modelMapper.get_pmaxVal(), 0, modMatrix->getMaxVal(), "%.2f");
			ImGui::ColorEdit4("Min color", modelMapper.get_pminCol(), ImGuiColorEditFlags_Float);
			ImGui::ColorEdit4("Max color", modelMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
			
			// // plot cut through sensitivity field
			
			// ImGui::SliderFloat("MinVal", modelMapper.get_pminVal(), 0, modMatrix->getMaxVal(), "%.1f");
			// ImGui::SliderFloat("MaxVal", modelMapper.get_pmaxVal(), 0, modMatrix->getMaxVal(), "%.1f");
			// ImGui::ColorEdit4("Min color", modelMapper.get_pminCol(), ImGuiColorEditFlags_Float);
			// ImGui::ColorEdit4("Max color", modelMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);

			// plot single model vector
			ImGui::PlotConfig conf;
			conf.values.ys = modMatrix->get_modelVec(rLvlModel, zLvlModel); // this line is optional
			//conf.values.ys = model->get_pimpResp();
			conf.values.count = modMatrix->getNtModel();
			conf.scale.min = -modMatrix->get_modVecAbsMax();
			conf.scale.max = modMatrix->get_modVecAbsMax();
			conf.tooltip.show = true;
			conf.tooltip.format = "x=%.2f, y=%.2f";
			conf.grid_x.show = false;
			conf.grid_y.show = false;
			conf.frame_size = ImVec2(500, 200);
			conf.line_thickness = 2.f;

			ImGui::Plot("plot", conf);
		}
		else
		{
			ImGui::Text("Model not built yet, nothing to show here");
		}
	}

	ImGui::End();
}

// interface to position differently absorbing spheres at different positions
void simulation::GenerateAbsorbanceVolume()
{
	// first generate empty matrix and fill with zeros
	// --> absorberMatrix
	absorberMatrix->allocMemory();
	absorberMatrix->setValue(0.0);

	// loop through absorbers and assign values
	for (uint8_t iAbsorber = 0; iAbsorber < pAbs.size(); iAbsorber++)
	{
		pointAbsorber currAbs = pAbs[iAbsorber];
		uint64_t posIdx[3];
		for (uint8_t iDim = 0; iDim < 3; iDim++)
			posIdx[iDim] = absorberMatrix->getIdx(currAbs.pos[iDim], iDim);
		absorberMatrix->setValue(posIdx, currAbs.absorbance);
	}

	absorberMatrix->calcMinMax();

	isAbsorberVolDefined = 1;

	return;
}


void simulation::SimulationSettings()
{
	ImGui::Begin("Simulation");
	ImGui::Text("Simulation settings");

	ImGui::InputFloat("Noise level", &noiseLevel);
	ImGui::InputFloat("SOS [m/s]", sett->get_pSOS());

	ImGui::InputFloat("Sampling frequency [MHz]", &fSampl);
	float rRes = sett->get_rRes() * 1e6;
	ImGui::InputFloat("rRes [microm]", &rRes);
	sett->set_rRes(rRes / 1e6);

	ImGui::Separator();
	ImGui::Columns(2);
	if (isAbsorberVolDefined)
		ImGui::Text("Yes");
	else
		ImGui::Text("No");
	ImGui::NextColumn();
	ImGui::Text("Absorbers defined?");

	ImGui::NextColumn();
	if (isModelBuilt)
		ImGui::Text("Yes");
	else
		ImGui::Text("No");
	ImGui::NextColumn();
	ImGui::Text("Model defined?");
	
	ImGui::Separator();
	ImGui::Columns(1);

	// disable button if volume not defined yet (isAbsorberVolDefined)
	if (!isAbsorberVolDefined)
    {
        ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
        ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    }

	if (ImGui::Button("Run simulation"))
	{
		// push all important properties over to reconSett
		for (uint8_t iDim = 0; iDim < 3; iDim++)
			sett->set_drRecon(iDim, absorberMatrix->getRes(iDim));

		sett->set_zLower(absorberMatrix->getOrigin(0));
		sett->set_zUpper(absorberMatrix->getMaxPos(0));
		sett->set_xLower(absorberMatrix->getOrigin(1));
		sett->set_xUpper(absorberMatrix->getMaxPos(1));
		sett->set_yLower(absorberMatrix->getOrigin(2));
		sett->set_yUpper(absorberMatrix->getMaxPos(2));

		// push all important settings over to fluence modeling
		// z positioning and dimensioning depends on 
		mcfieldprop->set_zMax(absorberMatrix->getMaxPos(0));
		mcfieldprop->set_dz(absorberMatrix->getRes(0));
		
		// mcfieldprop->set_rMax();
		// mcfieldprop->set_dr();

		// define sampling frequency for simulation
		signalMatrix->setRes(0, 1.0 / (fSampl * 1e6));

		recon.simulate();
		signalMatrix->calcMinMax();

		// add noise on top of simulation
		noiser.set_noiseLevel(noiseLevel * signalMatrix->get_maxAbsVal());
		noiser.addNoise(signalMatrix->get_pdata(), signalMatrix->get_nElements());

		isSimulationDone = 1;
	}

	if (!isAbsorberVolDefined)
    {
        ImGui::PopItemFlag();
        ImGui::PopStyleVar();	
        ImGui::Text("Absorber volume not defined yet.");
    }

	ImGui::Separator();
	ImGui::Text("Export functions");

	if (!isSimulationDone)
    {
        ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
        ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    }

	if (ImGui::Button("Export as vtk"))
	{
		signalMatrix->printInformation();
		signalMatrix->exportVtk("/home/hofmannu/simulatedData.vtk");
	}	
	if (ImGui::Button("Export as h5"))
	{
		signalMatrix->saveToFile("/home/hofmannu/simulatedData.h5");
	}

	if (!isSimulationDone)
    {
        ImGui::PopItemFlag();
        ImGui::PopStyleVar();
    }

    // print some information about volume if done
    if (isSimulationDone)
    {
    	PrintVolumeInformation(signalMatrix);
    	ImGui::Columns(1);
    	ImGui::SliderInt("zLayer", &currSliceZ, 0, signalMatrix->getDim(0) - 1);
		ImImagesc(signalMatrix->get_psliceZ((uint64_t) currSliceZ),
			signalMatrix->getDim(1), signalMatrix->getDim(2), &in_vol_texture, myMapper);
		
		int width = 550;
		int height = (float) width / signalMatrix->getLength(1) * signalMatrix->getLength(2);
		ImGui::Image((void*)(intptr_t)in_vol_texture, ImVec2(width, height));
		
		ImGui::SliderInt("yLayer", &currSliceY, 0, signalMatrix->getDim(2) - 1);

		ImImagesc(signalMatrix->get_psliceY((uint64_t) currSliceY),
		 	signalMatrix->getDim(1), signalMatrix->getDim(0), &in_vol_texture_slice, myMapper);
		height = (float) width / signalMatrix->getLength(1) * signalMatrix->getLength(0) * 1495.0;
		ImGui::Image((void*)(intptr_t)in_vol_texture_slice, ImVec2(width, height));
		
		ImGui::SliderFloat("MinVal", myMapper.get_pminVal(), signalMatrix->getMinVal(), signalMatrix->getMaxVal(), "%.1f");
		ImGui::SliderFloat("MaxVal", myMapper.get_pmaxVal(), signalMatrix->getMinVal(), signalMatrix->getMaxVal(), "%.1f");
		ImGui::ColorEdit4("Min color", myMapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", myMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
	
    }

	ImGui::End();
	return;
}

// interface to define the position of different absorbers in our field
void simulation::AbsorbanceMatrixWindow()
{
	ImGui::Begin("Absorber volume", &show_absorbance_window);

	float originMm[3];
	int dim[3];
	float resMm[3];
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		originMm[iDim] = absorberMatrix->getOrigin(iDim) * 1000;
		dim[iDim] = absorberMatrix->getDim(iDim);
		resMm[iDim] = absorberMatrix->getRes(iDim) * 1000;
	}

	ImGui::Text("Volume definition"); 
	ImGui::InputFloat3("Origin [mm, z/x/y]", originMm);
	ImGui::InputInt3("Dimension [mm, z/x/y]", dim);
	ImGui::InputFloat3("Resolution [mm, z/x/y]", resMm);
	
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		absorberMatrix->setOrigin(iDim, originMm[iDim] / 1000);
		absorberMatrix->setRes(iDim, resMm[iDim] / 1000);
		absorberMatrix->setDim(iDim, dim[iDim]);
	}

	ImGui::Separator();
	ImGui::Text("Point absorber definition");

	ImGui::Columns(2);
	for (unsigned int iAbsorber = 0; iAbsorber < pAbs.size(); iAbsorber++)
	{
		pointAbsorber currAbs = pAbs[iAbsorber];
		float posMm [3];
		posMm[0] = currAbs.pos[0] * 1000;
		posMm[1] = currAbs.pos[1] * 1000;
		posMm[2] = currAbs.pos[2] * 1000;

		ImGui::PushID(iAbsorber);
		ImGui::InputFloat3("Position [mm]", posMm);
		ImGui::NextColumn();
		ImGui::InputFloat("Absorbance", &currAbs.absorbance);
		ImGui::SameLine();
		currAbs.pos[0] = posMm[0] / 1000;
		currAbs.pos[1] = posMm[1] / 1000;
		currAbs.pos[2] = posMm[2] / 1000;
		pAbs[iAbsorber] = currAbs;
		
		if (ImGui::Button("x"))
		{
			pAbs.erase(pAbs.begin() + iAbsorber);
		}
		ImGui::PopID();
		ImGui::NextColumn();
	}
	ImGui::Columns(1);

	// add another absorber to list above
	if (ImGui::Button("Add absorber"))
	{
		pointAbsorber newAbs;
		newAbs.pos[0] = 0;
		newAbs.pos[1] = 0;
		newAbs.pos[2] = 0;
		newAbs.absorbance = 1;
		pAbs.push_back(newAbs);
	}
	
	ImGui::Separator();
	ImGui::Text("Generate volumetric representation");
	// generate volume out of defined absorbers
	if (ImGui::Button("Generate volume"))
	{
		GenerateAbsorbanceVolume();
	}

	if (isAbsorberVolDefined)
	{
		PrintVolumeInformation(absorberMatrix);
	}

	ImGui::End();
}

void simulation::PrintVolumeInformation(volume* vol)
{
	ImGui::Columns(2);
	ImGui::Text("%d x %d x %d", vol->getDim(0), 
		vol->getDim(1), vol->getDim(2));
	ImGui::NextColumn();
	ImGui::Text("Dimensions");
	ImGui::NextColumn();
	ImGui::Text("%f x %f x %f", vol->getRes(0) * 1e6,
		vol->getRes(1) * 1e6, vol->getRes(2) * 1e6);
	ImGui::NextColumn();
	ImGui::Text("Resolution [microm/micros]");
	ImGui::NextColumn();
	ImGui::Text("%f ... %f", vol->getMinVal(), vol->getMaxVal());
	ImGui::NextColumn();
	ImGui::Text("Data range");
	ImGui::NextColumn();
	ImGui::Text("%f x %f x %f", vol->getOrigin(0) * 1e6,
		vol->getOrigin(1) * 1e6, vol->getOrigin(2) * 1e6);
	ImGui::NextColumn();
	ImGui::Text("Origin [micros/microm]");
	return;
}

// define all the fluence settings for our simulation
void simulation::FluenceWindow()
{

	ImGui::Begin("Fluence settings", &show_fluence_window);
	ImGui::Checkbox("Fluence modeling", &flagFluenceModeling);
	ImGui::SameLine();
	HelpMarker("If disabled, fluence will be assumed constant");
	if (flagFluenceModeling)
	{
		ImGui::InputFloat("z level surface [mm]", &zLevelSurfaceMm);
		mcfieldprop->set_zBound(zLevelSurfaceMm / 1000);

		if (ImGui::CollapsingHeader("Fiber definition"))
		{
			// arfiber
			float dCore = arfiber->get_dCore() * 1e6; // core diameter in microm
			float na = arfiber->get_numAp();
			string name = arfiber->get_name();
			ImGui::Text("Name: %s\n", name);
			ImGui::InputFloat("Numerical aperture", &na);
			ImGui::InputFloat("Core diameter [microm]", &dCore);
			arfiber->set_dCore(dCore / 1e6); // convert back to m
			arfiber->set_numAp(na);
		}

		if (ImGui::CollapsingHeader("Tissue definition"))
		{
			ImGui::InputFloat("Absorption coefficient [1/m]", tissueprop->get_pmua());
			ImGui::InputFloat("Scattering coefficient [1/m]", tissueprop->get_pmus());
			ImGui::InputFloat("Anisotropy coefficient [1]", tissueprop->get_pg());
			ImGui::InputFloat("Refractive index [1]", tissueprop->get_pn());
			tissueprop->calc_albedo();
		}

		if (ImGui::CollapsingHeader("Preview of fluence result"))
		{
			if (!isSimulationDone) // if simulation is not done yet
			{
				ImGui::Text("Please run simulation first");
			}
			else
			{
				ImGui::Text("Not implemented yet");		
			}
		}
	}

	ImGui::End();
	return;
}

void simulation::MainDisplayCode()
{
	AbsorbanceMatrixWindow();
	TransducerSettingsWindow();
	SimulationSettings();
	FluenceWindow();
	return;
}
