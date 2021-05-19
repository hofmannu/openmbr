#include "interface.cuh"

interface::interface()
{
	inputDataVol = recon.get_ppreprocVol();
	trans = recon.get_ptransProp();
	model = recon.get_ptransModel();
	reconVol = recon.get_preconVol();
	gauss = model->get_pgauss();
	fieldProps = model->get_pfieldProp();
	
	// create con
	modelCutMapper.set_mapType(1);

}

interface::~interface()
{

}

void interface::InitWindow(int *argcp, char**argv)
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
}

void interface::ReconSettingWindow()
{
	ImGui::Begin("Reconstruction settings", &show_recon_window);
	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Reconstruction geometry"))
	{
		float zReconMm[2];
		float xReconMm[2];
		float yReconMm[2];
		zReconMm[0] = sett.get_zLower() * 1e3;
		zReconMm[1] = sett.get_zUpper() * 1e3;
		xReconMm[0] = sett.get_xLower() * 1e3;
		xReconMm[1] = sett.get_xUpper() * 1e3;
		yReconMm[0] = sett.get_yLower() * 1e3;
		yReconMm[1] = sett.get_yUpper() * 1e3;
		ImGui::InputFloat2("z crop [mm]", zReconMm);
		ImGui::SameLine(); HelpMarker("Reconstruction range in z direction.");
		ImGui::InputFloat2("x crop [mm]", xReconMm);
		ImGui::SameLine(); HelpMarker("Reconstruction range in x direction.");
		ImGui::InputFloat2("y crop [mm]", yReconMm);
		ImGui::SameLine(); HelpMarker("Reconstruction range in y direction.");
		sett.set_zLower(zReconMm[0] / 1e3);
		sett.set_zUpper(zReconMm[1] / 1e3);
		sett.set_xLower(xReconMm[0] / 1e3);
		sett.set_xUpper(xReconMm[1] / 1e3);
		sett.set_yLower(yReconMm[0] / 1e3);
		sett.set_yUpper(yReconMm[1] / 1e3);

		float resMicroM[3];
		resMicroM[0] = sett.get_dr(0) * 1e6;
		resMicroM[1] = sett.get_dr(1) * 1e6;
		resMicroM[2] = sett.get_dr(2) * 1e6;
		ImGui::InputFloat3("res (z, x, y) [microm]", resMicroM);
		for (uint8_t iDim = 0; iDim < 3; iDim++)
			sett.set_drRecon(iDim, resMicroM[iDim] / 1e6);
		ImGui::SameLine(); HelpMarker("Resolution of output dataset after reconstruction.");
	}
	ImGui::Separator();
	ImGui::InputFloat("SOS [m/s]", sett.get_pSOS());
	ImGui::SameLine(); HelpMarker("Speed of sound assumed in coupling medium");
	ImGui::InputFloat("Threshold [m]", sett.get_pthreshold());
	ImGui::InputFloat("Reg param [m]", sett.get_plambdaReg());
	
	int nIter = sett.get_nIter();
	ImGui::InputInt("nIter [m]", &nIter);
	if (nIter > 0)
		sett.set_nIter(nIter);

	ImGui::InputFloat("rRes [m]", sett.get_prRes());
	ImGui::SameLine(); HelpMarker("rRes");
	ImGui::Checkbox("Regularization", sett.get_pflagRegularization());
	ImGui::End();
	return;
}
// define the most important settings of the transducer
void interface::TransducerSettingsWindow()
{
	ImGui::Begin("Model Building", &show_trans_window);
	
	std::string name = trans->getName();
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
	if (ImGui::Button("Load properties"))
	{
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

	if (ImGui::Button("Save properties"))
	{
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
	if (ImGui::CollapsingHeader("Transducer geometry"))
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

	// here we define the most important field properties
	if (ImGui::CollapsingHeader("Field properties"))
	{

		float drmm = fieldProps->get_dr() * 1e6;
		ImGui::InputFloat("radial resolution [microm]", &drmm);
		fieldProps->set_dr(drmm / 1e6);

		float rMax = fieldProps->get_dr() * ((float) fieldProps->get_nr() - 1);
		rMax = rMax * 1000;
		ImGui::InputFloat("radial extend [mm]", &rMax);
		rMax = rMax / (drmm / 1000) + 1;
		fieldProps->set_nr((unsigned int) rMax);

		float z0mm = fieldProps->get_z0() * 1e3;
		ImGui::InputFloat("z start [mm]", &z0mm);
		fieldProps->set_z0(z0mm / 1000);

		float dzmicrom = fieldProps->get_dz() * 1e6;
		ImGui::InputFloat("axial resolution [microm]", &dzmicrom);
		fieldProps->set_dz(dzmicrom / 1e6);

		float z1mm = z0mm + ((float) fieldProps->get_nz() - 1) * (dzmicrom / 1e3); 
		ImGui::InputFloat("z stop [mm]", &z1mm);
		float deltaZ = (z1mm - z0mm) / (dzmicrom / 1000) + 1;
		fieldProps->set_nz((unsigned int) deltaZ);
	}

	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Impulse signal"))
	{
		// todo: diable button if not found
		if (ImGui::Button("Load impulse signal"))
		{
			ImGuiFileDialog::Instance()->OpenDialog("ReadImpulseDialog", 
			"Choose File", ".h5\0", ".");
		}

		if (ImGuiFileDialog::Instance()->FileDialog("ReadImpulseDialog")) 
		{
			if (ImGuiFileDialog::Instance()->IsOk == true)
			{
				const string filePathName = ImGuiFileDialog::Instance()->GetFilepathName();
				model->loadImpulseResponse(filePathName);
				isImpulseLoaded = 1;
			}
			ImGuiFileDialog::Instance()->CloseDialog("ReadImpulseDialog");
		}

		// plot impulse signal here after loading		
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
			conf.grid_x.show = true;
			conf.grid_y.show = true;
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
	if (ImGui::CollapsingHeader("Optical modeling"))
	{
		ImGui::InputFloat("wavelength [nm]", gauss->get_pwavelength());
		ImGui::InputFloat("NA [1]", gauss->get_pna());
		ImGui::InputFloat("n [1]", gauss->get_pn());
		ImGui::InputFloat("I0 [1]", gauss->get_pI0());
	}

	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Model building"))
	{
		int nElements = model->getNElements();
		ImGui::InputInt("Number of elements", &nElements);
		model->setNElements(nElements);

		ImGui::SameLine(); HelpMarker("Number of elements used for discretization");

		ImGui::Columns(2);
		ImGui::Text("Is impulse signal defined?");
		ImGui::NextColumn();
		if (isImpulseLoaded)
			ImGui::Text("y");
		else
			ImGui::Text("n");
		ImGui::NextColumn();

		ImGui::Text("Is model built?");
		ImGui::NextColumn();
		if (isModelBuilt)
			ImGui::Text("y");
		else
			ImGui::Text("n");
		ImGui::NextColumn();

		ImGui::Columns(1);
		if (ImGui::Button("Build model"))
		{
			model->setTransProp(*trans);
			if (!isImpulseLoaded) // try to load impulse response if not done yet
				model->loadImpulseResponse();
			
			model->buildGPUModel();
			modMatrix = model->getModel();

			// if requested convolve here with optical model
			gauss->define_r(modMatrix->get_dr(), (uint32_t) modMatrix->get_nr());
			gauss->define_z(modMatrix->get_dz(), modMatrix->get_z0(), modMatrix->get_nz());
			gauss->convolve_model(modMatrix->get_data(), modMatrix->get_nt());

			modMatrix->calcSensField();
			isModelBuilt = 1;
		}

	}

	
	if (isModelBuilt && ImGui::CollapsingHeader("Model preview"))
	{
		// make two sliders for r and z position selection
		
		// plot sensitivity field
		ImImagesc(modMatrix->get_psensField(), modMatrix->get_nr(), 
			modMatrix->get_nz(), &mod_texture, modelMapper);
		
		int width = 250; 
		int height = 500;
		ImGui::Image((void*)(intptr_t)mod_texture, ImVec2(width, height));
		ImGui::SameLine();
		
		float* modSlice = modMatrix->get_pmodSlice(zLvlModel);
		modelCutMapper.set_maxVal(modMatrix->get_maxValSlice());
		modelCutMapper.set_minVal(modMatrix->get_minValSlice());
		ImImagesc(modSlice, modMatrix->get_nr(), 
		 	modMatrix->get_nt(), &modcut_texture, modelCutMapper);
		ImGui::Image((void*)(intptr_t)modcut_texture, ImVec2(width, height));

		float zDepthMm = zLvlModel * 1000;
		float rDepthMm = rLvlModel * 1000;
		ImGui::SliderFloat("z depth [mm]", &zDepthMm, modMatrix->get_zMin() * 1000, modMatrix->get_zMax() * 1000, "%.3f");
		ImGui::SliderFloat("r depth [mm]", &rDepthMm, 0.0, modMatrix->get_rMax() * 1000, "%.3f");
		zLvlModel = zDepthMm / 1000;
		rLvlModel = rDepthMm / 1000;

		ImGui::SliderFloat("MinVal", modelMapper.get_pminVal(), 0, modMatrix->getMaxVal(), "%.2f");
		ImGui::SliderFloat("MaxVal", modelMapper.get_pmaxVal(), 0, modMatrix->getMaxVal(), "%.2f");
		ImGui::ColorEdit4("Min color", modelMapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", modelMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
		
		// // plot cut through sensitivity field
		
		// plot single model vector
		ImGui::PlotConfig conf;
		conf.values.ys = modMatrix->get_modelVec(rLvlModel, zLvlModel); 
		conf.values.count = modMatrix->get_nt();
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

	if (isModelBuilt && ImGui::CollapsingHeader("Model export"))
	{
		static char filePath [64];
		ImGui::InputText("Path", filePath, 64);
		string filePathString = filePath;

		ImGui::Columns(2);
		if (ImGui::Button("Export as h5"))
		{
			modMatrix->saveToFile(filePathString);
		}
		ImGui::NextColumn();
		if (ImGui::Button("Export as vtk"))
		{
			modMatrix->exportSensFieldVtk(filePathString);
		}
		ImGui::Columns(1);
	}

	ImGui::End();
}

void interface::InputVolViz()
{
	ImGui::Begin("Data Loader", &show_data_loader);
	ImGui::NextColumn();
	if (ImGui::Button("Load data")){
		ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", 
			"Choose File", ".h5\0", ".");
	}
	if (ImGuiFileDialog::Instance()->FileDialog("ChooseFileDlgKey")) 
	{
		if (ImGuiFileDialog::Instance()->IsOk == true)
		{
			std::string inputFilePath = ImGuiFileDialog::Instance()->GetFilepathName();
			inputDataVol->readFromFile(inputFilePath);
			inputDataVol->calcMinMax();
			isDataSetDefined = 1;
			
		}
		ImGuiFileDialog::Instance()->CloseDialog("ChooseFileDlgKey");
	}
	ImGui::NextColumn();

	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Size and dimension"))
	{
		ImGui::Columns(2);	
		// print information about dataset if defined
		ImGui::Text("Size"); ImGui::NextColumn();
		
		if (isDataSetDefined)
		{
			ImGui::Text("%i x %i x %i", 
				inputDataVol->getDim(0), 
				inputDataVol->getDim(1), 
				inputDataVol->getDim(2)); 
		}
		ImGui::NextColumn();

		ImGui::Text("Resolution [microm]"); ImGui::NextColumn();
		if (isDataSetDefined)
		{
			ImGui::Text("%.2f x %.2f", 
				inputDataVol->getRes(1) * 1e6, 
				inputDataVol->getRes(2) * 1e6);
		} ImGui::NextColumn();

		ImGui::Text("Sampling freqeucny [MHz]"); ImGui::NextColumn();
		if (isDataSetDefined)
		{
			ImGui::Text("%.2f",	1 / inputDataVol->getRes(0) / 1e6);
		} ImGui::NextColumn();

		ImGui::Text("t range [micros]"); ImGui::NextColumn();
		if (isDataSetDefined)
		{
			ImGui::Text("%.2f ... %.2f", 
				inputDataVol->getMinPos(0) * 1e6, 
				inputDataVol->getMaxPos(0) * 1e6);
		}
		ImGui::NextColumn();

		ImGui::Text("x range [mm]"); ImGui::NextColumn();
		if (isDataSetDefined)
		{
			ImGui::Text("%.2f ... %.2f", 
				inputDataVol->getMinPos(1) * 1e3, 
				inputDataVol->getMaxPos(1) * 1e3);
		}
		ImGui::NextColumn();

		ImGui::Text("y range [mm]"); ImGui::NextColumn();
		if (isDataSetDefined)
		{
			ImGui::Text("%.2f ... %.2f", 
				inputDataVol->getMinPos(2) * 1e3, 
				inputDataVol->getMaxPos(2) * 1e3);
		}
		ImGui::NextColumn();
	}

	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Information"))
	{
		ImGui::Columns(2);

		ImGui::Text("Maximum value"); ImGui::NextColumn();
		if (isDataSetDefined)
		{
			ImGui::Text("%f", inputDataVol->getMinVal());

		} ImGui::NextColumn();

		ImGui::Text("Minimum value"); ImGui::NextColumn();
		if (isDataSetDefined)
		{
			ImGui::Text("%f", inputDataVol->getMaxVal());
		} ImGui::NextColumn();

	}
	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Operations"))
	{
		ImGui::Text("Not implemented yet");
		ImGui::Text("- Normalize");
		ImGui::Text("- Remove infs");
		ImGui::Text("- Scale");
	}

	if (ImGui::CollapsingHeader("Export functions"))
	{
		static char filePath [64];
		ImGui::InputText("Path", filePath, 64);
		string filePathString = filePath;

		if (ImGui::Button("Export as vtk"))
		{
			inputDataVol->exportVtk(filePathString);
		}
		if (ImGui::Button("Export as h5"))
		{
			inputDataVol->saveToFile(filePath);	
		}
	}

	// make imagesc here with plot below
	if (isDataSetDefined)
	{
		ImGui::SliderInt("zLayer", &currSliceZ, 0, inputDataVol->getDim(0) - 1);
		float tSlice = inputDataVol->getPos(currSliceZ, 0);
		ImGui::Text("Time of current layer: %f micros", tSlice * 1e6);
		ImGui::Text("Approximated z layer: %f mm", tSlice * 1490e3);
		ImImagesc(inputDataVol->get_psliceZ((uint64_t) currSliceZ),
			inputDataVol->getDim(1), inputDataVol->getDim(2), &in_vol_texture, myMapper);
		
		int width = 550;
		int height = (float) width / inputDataVol->getLength(1) * inputDataVol->getLength(2);
		ImGui::Image((void*)(intptr_t)in_vol_texture, ImVec2(width, height));
		
		ImGui::SliderInt("yLayer", &currSliceY, 0, inputDataVol->getDim(2) - 1);

		ImImagesc(inputDataVol->get_psliceY((uint64_t) currSliceY),
		 	inputDataVol->getDim(1), inputDataVol->getDim(0), &in_vol_texture_slice, myMapper);
		height = (float) width / inputDataVol->getLength(1) * inputDataVol->getLength(0) * 1495.0;
		ImGui::Image((void*)(intptr_t)in_vol_texture_slice, ImVec2(width, height));
		
		ImGui::SliderFloat("MinVal", myMapper.get_pminVal(), inputDataVol->getMinVal(), inputDataVol->getMaxVal(), "%.1f");
		ImGui::SliderFloat("MaxVal", myMapper.get_pmaxVal(), inputDataVol->getMinVal(), inputDataVol->getMaxVal(), "%.1f");
		ImGui::ColorEdit4("Min color", myMapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", myMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
	
			
	}

	ImGui::End();
	return;
}

void interface::ImImagesc(
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

void interface::OutputVolViz()
{
	ImGui::Begin("Data exporter", &show_recon_vol);
	ImGui::NextColumn();
	
	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Size and dimension"))
	{
		ImGui::Columns(2);	
		// print information about dataset if defined
		ImGui::Text("Size"); ImGui::NextColumn();
		
		if (isReconVolDefined)
		{
			ImGui::Text("%i x %i x %i", 
				reconVol->getDim(0), 
				reconVol->getDim(1), 
				reconVol->getDim(2)); 
		}
		ImGui::NextColumn();

		ImGui::Text("Resolution [microm]"); ImGui::NextColumn();
		if (isReconVolDefined)
		{
			ImGui::Text("%.2f x %.2f x %.2f", 
				reconVol->getRes(1) * 1e6, 
				reconVol->getRes(1) * 1e6, 
				reconVol->getRes(2) * 1e6);
		} ImGui::NextColumn();

		ImGui::Text("z range [micros]"); ImGui::NextColumn();
		if (isReconVolDefined)
		{
			ImGui::Text("%.2f ... %.2f", 
				reconVol->getMinPos(0) * 1e3, 
				reconVol->getMaxPos(0) * 1e3);
		}
		ImGui::NextColumn();

		ImGui::Text("x range [mm]"); ImGui::NextColumn();
		if (isReconVolDefined)
		{
			ImGui::Text("%.2f ... %.2f", 
				reconVol->getMinPos(1) * 1e3, 
				reconVol->getMaxPos(1) * 1e3);
		}
		ImGui::NextColumn();

		ImGui::Text("y range [mm]"); ImGui::NextColumn();
		if (isReconVolDefined)
		{
			ImGui::Text("%.2f ... %.2f", 
				reconVol->getMinPos(2) * 1e3, 
				reconVol->getMaxPos(2) * 1e3);
		}
		ImGui::NextColumn();
	}

	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Information"))
	{
		ImGui::Columns(2);

		ImGui::Text("Maximum value"); ImGui::NextColumn();
		if (isReconVolDefined)
		{
			ImGui::Text("%f", reconVol->getMinVal());

		} ImGui::NextColumn();

		ImGui::Text("Minimum value"); ImGui::NextColumn();
		if (isReconVolDefined)
		{
			ImGui::Text("%f", reconVol->getMaxVal());
		} ImGui::NextColumn();

	}
	ImGui::Columns(1);
	if (ImGui::CollapsingHeader("Operations"))
	{
		ImGui::Text("Not implemented yet");
		ImGui::Text("- Normalize");
		ImGui::Text("- Remove infs");
		ImGui::Text("- Scale");
	}

	if (ImGui::CollapsingHeader("Export functions"))
	{
		static char filePath [64];
		ImGui::InputText("Path", filePath, 64);
		string filePathString = filePath;

		ImGui::Columns(2);

		if (ImGui::Button("Export as vtk"))
		{
			reconVol->exportVtk(filePathString);
		}
		ImGui::NextColumn();
		if (ImGui::Button("Export as h5"))
		{
			reconVol->saveToFile(filePath);	
		}
		ImGui::NextColumn();
		ImGui::Columns(1);
	}

	if (isReconVolDefined)
	{
		// make imagesc here with plot below
		ImGui::SliderFloat("z [mm]", &currSliceZRec, 
			reconVol->getMinPos(0) * 1e3, reconVol->getMaxPos(0) * 1e3);
		const float currSliceZRecM = currSliceZRec * 1.0e-3;
		ImImagesc(reconVol->get_psliceZ(currSliceZRecM),
			reconVol->getDim(1), reconVol->getDim(2), &rec_vol_texture, reconVolMapper);
		
		int width = 550;
		int height = (float) width / reconVol->getLength(1) * reconVol->getLength(2);
		ImGui::Image((void*)(intptr_t)rec_vol_texture, ImVec2(width, height));
		
		ImGui::SliderFloat("y [mm]", &currSliceYRec, 
			reconVol->getMinPos(2) * 1e3, reconVol->getMaxPos(2) * 1e3);
		const float currSliceYRecM = currSliceYRec * 1.0e-3;
		ImImagesc(reconVol->get_psliceY(currSliceYRecM),
	 		reconVol->getDim(1), reconVol->getDim(0), &rec_vol_texture_slice, reconVolMapper);
		
		height = (float) width / reconVol->getLength(1) * reconVol->getLength(0);
		ImGui::Image((void*)(intptr_t)rec_vol_texture_slice, ImVec2(width, height));
			
		ImGui::SliderFloat("MinVal", reconVolMapper.get_pminVal(), reconVol->getMinVal(), reconVol->getMaxVal(), "%.9f");
		ImGui::SliderFloat("MaxVal", reconVolMapper.get_pmaxVal(), reconVol->getMinVal(), reconVol->getMaxVal(), "%.9f");
		ImGui::ColorEdit4("Min color", reconVolMapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", reconVolMapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
	}

	ImGui::End();
	return;
}


void interface::Reconstructor()
{
	ImGui::Begin("Reconstructor", &show_reconstructor);
	ImGui::Columns(2);
	if (!isDataSetDefined)
		ImGui::Text("No"); 
	else
		ImGui::Text("Yes");

	ImGui::NextColumn();
	ImGui::Text("Dataset loaded?"); 
	ImGui::SameLine(); HelpMarker("Please load dataset using Data loader first");
	
	ImGui::NextColumn(); 
	if (!isTransducerValid)
		ImGui::Text("No");
	else
		ImGui::Text("Yes");

	ImGui::Text("Current iterat");

	ImGui::NextColumn();
	ImGui::Text("Transducer valid?"); 
	ImGui::SameLine(); HelpMarker("Please define a valid transducer first.");

	ImGui::Columns(1);
	if (ImGui::Button("Reconstruct")){
		// run reconstruction
		recon.set_sett(sett);
		recon.recon();
		reconVol->calcMinMax();
		isReconVolDefined = 1;
		show_recon_window = 1;
		// get model matrix back for gui
		modMatrix = model->getModel();
		modMatrix->calcSensField();
		isModelBuilt = 1;

		// update slices to center of volume to avoid bad range access
		currSliceYRec = reconVol->getCenterPos(2);
		currSliceZRec = reconVol->getCenterPos(1);
	}

	ImGui::End();
	return;
}

void interface::ModelBuilder()
{

	return;
}

void interface::MainDisplayCode()
{
	// if (ImGui::BeginMenuBar())
	// {
	// 	if (ImGui::BeginMenu("Menu"))
 //        {
 //        	ImGui::MenuItem("Transducer settings", NULL, &show_trans_window);
 //        	ImGui::EndMenu();
 //        }
 //        ImGui::EndMenuBar();
	// }


	ReconSettingWindow();
	TransducerSettingsWindow();
	InputVolViz();
	OutputVolViz();
	Reconstructor();
	
	return;
}
