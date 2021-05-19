#include "interface.cuh"

interface::interface()
{
	arfiber = sim.get_pfiber();
	field = sim.get_pfield();
	simprop = sim.get_psim();
	tissue = sim.get_ptissue();
}

interface::~interface()
{

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

void interface::InitWindow(int *argcp, char**argv){

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


// interfacing for important properties like fiber stuff, field size etc
void interface::Properties()
{

	ImGui::Begin("Simulation properties", &show_properties_window);

	// fiber properties
	float dCore = arfiber->get_dCore() * 1e6; // core diameter in microm
	float na = arfiber->get_numAp();
	string name = arfiber->get_name();
	ImGui::Text("Fiber properties");
	ImGui::Text("Name: %s\n", name);
	ImGui::InputFloat("Numerical aperture", &na);
	ImGui::InputFloat("Core diameter [microm]", &dCore);
	arfiber->set_dCore(dCore / 1e6); // convert back to m
	arfiber->set_numAp(na);
	
	// loading and saving of predefined fiber properties
	ImGui::Columns(2);
	if (ImGui::Button("Save"))
	{
		// nothing here yet
	}
	ImGui::NextColumn();
	if (ImGui::Button("Load"))
	{
		// nothing here yet
	}
	ImGui::NextColumn();

	ImGui::Columns(1);
	
	ImGui::Separator();
	ImGui::Text("Field properties");
	float rRes = field->get_dr() * 1e6; // im microm
	float zRes = field->get_dz() * 1e6; // in microm
	float rMax = field->get_rMax() * 1e3; // in mm
	float zMax = field->get_zMax() * 1e3; // in mm
	float zBound = field->get_zBound() * 1e3; // in mm
	ImGui::InputFloat("Radial resolution [microm]", &rRes);
	ImGui::InputFloat("Axial resolution [microm]", &zRes);
	ImGui::InputFloat("Maximum radius [mm]", &rMax);
	ImGui::InputFloat("Maximum z distance [mm]", &zMax);
	ImGui::InputFloat("Boundary position [mm]", &zBound);
	
	// push values back to field class after converting back to m
	field->set_dz(zRes * 1e-6);
	field->set_dr(rRes * 1e-6);
	field->set_zMax(zMax * 1e-3);
	field->set_zBound(zBound * 1e-3);
	field->set_rMax(rMax * 1e-3);


	ImGui::Separator();
	ImGui::Text("Simulation properties");
	int nPhotons = simprop->get_nPhotons();
	int nPPerThread = simprop->get_nPPerThread();
	ImGui::InputInt("Number of photons", &nPhotons); ImGui::SameLine();
	HelpMarker("Number of photon packages simulated using Monte Carlo.s");
	ImGui::InputInt("Number of photons per thread", &nPPerThread);
	simprop->set_nPPerThread(nPPerThread);
	simprop->set_nPhotons(nPhotons);
	ImGui::Text("True number of simulated photons: %d", simprop->get_nPhotonsTrue());
	ImGui::Text("Number of blocks: %d, threads per block: %d", 
		simprop->get_nBlocks(), simprop->get_threadsPerBlock());

	if (ImGui::Button("Run"))
	{
		// run simulation
		sim.run();
		is_output_defined = 1;
		lastNr = field->get_nr();
		lastNz = field->get_nz();
	}

	ImGui::End();
	return;
}

// c implementation of matlabs imagesc
void interface::ImImagesc(
	const float* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap){
	
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

void interface::TissueProperties()
{
	ImGui::Begin("Optical properties", &show_tissue_properties);

	ImGui::Text("Tissue");
	
	float mua = tissue->get_mua();
	float mus = tissue->get_mus();
	ImGui::InputFloat("Absorption coefficient [1/m]", &mua);
	ImGui::InputFloat("Scattering coefficient [1/m]", &mus);
	ImGui::InputFloat("Anisotropy", tissue->get_pg());
	ImGui::InputFloat("Refractive index", tissue->get_pn());
	tissue->set_mua(mua);
	tissue->set_mus(mus);


	if (ImGui::Button("Load preset"))
	{

	}

	ImGui::End();
	return;
}

void interface::Result()
{
	ImGui::Begin("Result", &show_results);

	// make imagesc here with plot below
	int width = 550;
	int height = (float) width / field->get_rMax() * field->get_zExtend();
	
	ImGui::Checkbox("Logarthmic scale", &flagLogPlot);

	if (flagLogPlot) // use logarithmic plot in this case
	{
		ImImagesc(sim.get_pheat_log(),
	 		lastNr, lastNz, &fluence_texture, fluence_mapper);
		ImGui::Image((void*)(intptr_t)fluence_texture, ImVec2(width, height));
		
		ImGui::SliderFloat("MinVal", fluence_mapper.get_pminVal(), sim.get_minValLog(), sim.get_maxValLog(), "%.9f");
		ImGui::SliderFloat("MaxVal", fluence_mapper.get_pmaxVal(), sim.get_minValLog(), sim.get_maxValLog(), "%.9f");
		ImGui::ColorEdit4("Min color", fluence_mapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", fluence_mapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
	}
	else // normal plot
	{ 
		ImImagesc(sim.get_pheat(),
	 		field->get_nr(), field->get_nz(), &fluence_texture, fluence_mapper);
		ImGui::Image((void*)(intptr_t)fluence_texture, ImVec2(width, height));
		
		ImGui::SliderFloat("MinVal", fluence_mapper.get_pminVal(), sim.get_minVal(), sim.get_maxVal(), "%.9f");
		ImGui::SliderFloat("MaxVal", fluence_mapper.get_pmaxVal(), sim.get_minVal(), sim.get_maxVal(), "%.9f");
		ImGui::ColorEdit4("Min color", fluence_mapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", fluence_mapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
	}


	if (ImGui::Button("Export as vtk"))
	{
		sim.exportVtk("/home/hofmannu/fluenceExport.vtk");
	}
	if (ImGui::Button("Export as h5"))
	{

	}

	ImGui::End();
	return;
}

void interface::MainDisplayCode()
{

	// nothing here yet
	Properties();
	TissueProperties();
	if (is_output_defined)
	{
		Result();
	}
	return;
}