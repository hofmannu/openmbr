#include "gui.h"

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
	SettingsWindow();
	ModelWindow();
	DataLoader();
}


void gui::SettingsWindow()
{
	ImGui::Begin("Reconstruction settings");

	float regStength = sett->get_regStrength();
	ImGui::InputFloat("Reg strength", &regStength);
	sett->set_regStrength(regStength);

	int nIter = sett->get_nIter();
	ImGui::InputInt("Number of iterations", &nIter);
	sett->set_nIter(nIter);

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

	ImGui::Columns(2);
	ImGui::Text("Model loaded?"); ImGui::NextColumn();

	if (mod->get_isLoaded())
	{
		ImGui::Text("y");
	}
	else
	{
		ImGui::Text("n");

	}

	ImGui::NextColumn();
	ImGui::Text("Data loaded?"); ImGui::NextColumn();
	if (isSigMatLoaded)
	{
		ImGui::Text("y");
	}
	else
	{
		ImGui::Text("n");

	}

	ImGui::NextColumn();
	ImGui::Columns(1);
	if (ImGui::Button("Reconstruct"))
	{
		rec.crop();
		rec.lsqr();
	};

	ImGui::End();
	return;
}

void gui::ModelWindow()
{
	ImGui::Begin("Model matrix");
	if (ImGui::Button("Load from file"))
	{
		ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", 
			"Choose File", ".h5\0", ".");
	}

	if (ImGuiFileDialog::Instance()->FileDialog("ChooseFileDlgKey")) 
	{
		if (ImGuiFileDialog::Instance()->IsOk == true)
		{
			std::string inputFilePath = ImGuiFileDialog::Instance()->GetFilepathName();
			mod->load_from_file(inputFilePath);
			mod->calc_sensField();

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
	
	if (mod->get_isLoaded())
	{
		// some general informations about model size
		if (ImGui::CollapsingHeader("Model information"))
		{
			ImGui::Columns(2);
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

void gui::DataLoader()
{
	ImGui::Begin("Data loader");
	if (ImGui::Button("Load input data"))
	{
		ImGuiFileDialog::Instance()->OpenDialog("testmeow", 
			"Choose File", ".h5\0", ".");
	}

	if (ImGuiFileDialog::Instance()->FileDialog("testmeow")) 
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
		ImGuiFileDialog::Instance()->CloseDialog("testmeow");
	}


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

			ImGui::Text("Size [x, y, t]"); ImGui::NextColumn();
					ImGui::Text("%i x %i x %i", 
				sigMat->get_dim(0), 
				sigMat->get_dim(1), 
				sigMat->get_dim(2)); 
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

			ImGui::Text("Maximum value"); ImGui::NextColumn();
			ImGui::Text("%.2f", sigMat->get_minVal());
			ImGui::NextColumn();

			ImGui::Text("Minimum value"); ImGui::NextColumn();
			ImGui::Text("%.2f", sigMat->get_maxVal());
			ImGui::NextColumn();
			ImGui::Columns(1);
		}
		

		if (ImGui::CollapsingHeader("Dataset preview"))
		{
			ImGui::SliderInt("tSlice", &currSliceZ, 0, sigMat->get_dim(2) - 1);
			float tSlice = sigMat->get_pos(currSliceZ, 2);
			ImGui::Text("Time of current layer: %f micros", tSlice * 1e6);
			ImGui::Text("Approximated z layer: %f mm", tSlice * sett->get_sos());
			ImImagesc(sigMat->get_psliceZ((uint64_t) currSliceZ),
				sigMat->get_dim(0), sigMat->get_dim(1), &inDataTexture, inDataMapper);
			
			int width = 550;
			int height = (float) width / sigMat->get_length(0) * sigMat->get_length(1);
			ImGui::Image((void*)(intptr_t)inDataTexture, ImVec2(width, height));
			
			ImGui::SliderInt("xSlice", &currSliceY, 0, sigMat->get_dim(0) - 1);

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
