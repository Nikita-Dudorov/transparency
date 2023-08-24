// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------

#define _USE_MATH_DEFINES

#include <glad/glad.h>

#include <cstdlib>
#include <cstdio>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include <exception>
#include <filesystem>

namespace fs = std::filesystem;

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/quaternion.hpp>

#include "Resources.h"
#include "Error.h"
#include "Console.h"
#include "IO.h"
#include "Scene.h"
#include "Image.h"
#include "Rasterizer.h"
#include "RayTracer.h"

using namespace std;

// Window parameters
static GLFWwindow * windowPtr = nullptr;
static std::shared_ptr<Scene> scenePtr;
static std::shared_ptr<Rasterizer> rasterizerPtr;
static std::shared_ptr<RayTracer> rayTracerPtr;

// Camera control variables
static glm::vec3 center = glm::vec3 (0.0); // To update based on the mesh position
static float meshScale = 1.0; // To update based on the mesh size, so that navigation runs at scale
static bool isRotating (false);
static bool isPanning (false);
static bool isZooming (false);
static double baseX (0.0), baseY (0.0);
static glm::vec3 baseTrans (0.0);
static glm::vec3 baseRot (0.0);

// Files
static std::string basePath;
static std::string meshFilename;
static std::string materialDirname;

// Raytraced rendering
static bool isDisplayRaytracing = false;

// Animation
static bool isAnimated = false;

void clear ();

void printHelp () {
	Console::print (std::string ("Help:\n") 
			  + "\tMouse commands:\n" 
			  + "\t* Left button: rotate camera\n" 
			  + "\t* Middle button: zoom\n" 
			  + "\t* Right button: pan camera\n" 
			  + "\tKeyboard commands:\n" 
   			  + "\t* ESC: quit the program\n"
   			  + "\t* H: print this help\n"
   			  + "\t* F12: reload GPU shaders\n"
   			  + "\t* S: save the raytraced image in "+ DEFAULT_RAYTRACED_IMAGE_OUTPUT_FILENAME + "\n"
   			  + "\t* D: save the rasterized image in "+ DEFAULT_RASTERIZED_IMAGE_OUTPUT_FILENAME + "\n"
   			  + "\t* F/Shift+F: increase/decrease field of view\n"
   			  + "\t* TAB: switch between rasterization and ray tracing display\n"
   			  + "\t* SPACE: execute ray tracing\n"
   			  + "\t* P: animate the scene\n"); // raytracing falls down after using this, bug to investigate
}

/// Adjust the ray tracer target resolution and runs it.
void raytrace () {
	int width, height;
	glfwGetWindowSize(windowPtr, &width, &height);
	std::chrono::high_resolution_clock clock;
	Console::print ("Start ray tracing at " + std::to_string (width) + "x" + std::to_string (height) + " resolution...");
	std::chrono::time_point<std::chrono::high_resolution_clock> before = clock.now();
	rayTracerPtr->setResolution (width, height);
	rayTracerPtr->render (scenePtr);
	std::chrono::time_point<std::chrono::high_resolution_clock> after = clock.now();
	double elapsedTime = (double)std::chrono::duration_cast<std::chrono::milliseconds>(after - before).count();
	Console::print ("Ray tracing executed in " + std::to_string(elapsedTime) + "ms");
}

/// Executed each time a key is entered.
void keyCallback (GLFWwindow * windowPtr, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_H) {
			printHelp ();
		} else if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE) {
			glfwSetWindowShouldClose (windowPtr, true); // Closes the application if the escape key is pressed
		} else if (action == GLFW_PRESS && key == GLFW_KEY_F12) {
			rasterizerPtr->loadShaderProgram (basePath);
			Console::print ("Reloading all shaders");
		} else if (action == GLFW_PRESS && key == GLFW_KEY_S) {
			Console::print ("Start writing raytraced rendering to " + DEFAULT_RAYTRACED_IMAGE_OUTPUT_FILENAME);
			rayTracerPtr->image()->save (DEFAULT_RAYTRACED_IMAGE_OUTPUT_FILENAME);
			Console::print ("End writing raytraced rendering");
		} else if (action == GLFW_PRESS && key == GLFW_KEY_D) {
			Console::print ("Start writing rasterized rendering to " + DEFAULT_RASTERIZED_IMAGE_OUTPUT_FILENAME);
			auto imagePtr = rasterizerPtr->generateImage();
			imagePtr->save (DEFAULT_RASTERIZED_IMAGE_OUTPUT_FILENAME);
			Console::print ("End writing rasterized rendering");
		} else if (action == GLFW_PRESS && key == GLFW_KEY_F) {
			scenePtr->camera()->setFoV (glm::clamp (scenePtr->camera()->getFoV () + (mods&GLFW_MOD_SHIFT ? -1.f : 1.f) * 5.f, 5.f, 120.f));
			Console::print ("Camera FoV set to: " + std::to_string (scenePtr->camera()->getFoV ()));
		} else if (action == GLFW_PRESS && key == GLFW_KEY_TAB) {
			isDisplayRaytracing = !isDisplayRaytracing;
			Console::print ("Raytracing display " + isDisplayRaytracing ? "activated" : "desactivated");
		} else if (action == GLFW_PRESS && key == GLFW_KEY_SPACE) {
			raytrace ();
		} else if (action == GLFW_PRESS && key == GLFW_KEY_P) {
			isAnimated = !isAnimated;
			Console::print ("Animation " + isAnimated ? "started" : "stopped");
		} else {
			printHelp ();
		}
	}
}

/// Called each time the mouse cursor moves
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
	int width, height;
	glfwGetWindowSize (windowPtr, &width, &height);
	float normalizer = static_cast<float> ((width + height)/2);
	float dx = static_cast<float> ((baseX - xpos) / normalizer);
	float dy = static_cast<float> ((ypos - baseY) / normalizer);
	if (isRotating) {
		glm::vec3 dRot (-dy * M_PI, dx * M_PI, 0.0);
		scenePtr->camera()->setRotation (baseRot + dRot);
		scenePtr->camera()->rayAt (0,0);
	} else if (isPanning) {
		scenePtr->camera()->setTranslation (baseTrans + meshScale * glm::vec3 (dx, dy, 0.0));
		scenePtr->camera()->rayAt (0,0);
	} else if (isZooming) {
		scenePtr->camera()->setTranslation (baseTrans + meshScale * glm::vec3 (0.0, 0.0, dy));
		scenePtr->camera()->rayAt (0,0);
	}
}

/// Called each time a mouse button is pressed
void mouseButtonCallback (GLFWwindow * window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    	if (!isRotating) {
    		isRotating = true;
    		isDisplayRaytracing = false;
    		glfwGetCursorPos (window, &baseX, &baseY);
    		baseRot = scenePtr->camera()->getRotation ();
        } 
    } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    	isRotating = false;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
    	if (!isPanning) {
    		isPanning = true;
    		isDisplayRaytracing = false;
    		glfwGetCursorPos (window, &baseX, &baseY);
    		baseTrans = scenePtr->camera()->getTranslation ();
        } 
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
    	isPanning = false;
    } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS) {
    	if (!isZooming) {
    		isZooming = true;
    		isDisplayRaytracing = false;
    		glfwGetCursorPos (window, &baseX, &baseY);
    		baseTrans = scenePtr->camera()->getTranslation ();
        } 
    } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_RELEASE) {
    	isZooming = false;
    }
}

/// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window. 
void windowSizeCallback (GLFWwindow * windowPtr, int width, int height) {
	scenePtr->camera()->setAspectRatio (static_cast<float>(width) / static_cast<float>(height));
	rasterizerPtr->setResolution (width, height);
	rayTracerPtr->setResolution (width, height);
}

void initGLFW () {
	// Initialize GLFW, the library responsible for window management
	if (!glfwInit ()) {
		Console::print ("ERROR: Failed to init GLFW");
		std::exit (EXIT_FAILURE);
	}

	// Before creating the window, set some option flags
	glfwWindowHint (GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint (GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint (GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint (GLFW_RESIZABLE, GL_TRUE);
	glfwWindowHint (GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // added to fix "Failed to open window" crash on MacOS

	// Create the window
	windowPtr = glfwCreateWindow (DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT, BASE_WINDOW_TITLE.c_str (), nullptr, nullptr);
	if (!windowPtr) {
		Console::print ("ERROR: Failed to open window");
		glfwTerminate ();
		std::exit (EXIT_FAILURE);
	}

	// Load the OpenGL context in the GLFW window using GLAD OpenGL wrangler
	glfwMakeContextCurrent (windowPtr);

	// Connect the callbacks for interactive control 
	glfwSetWindowSizeCallback (windowPtr, windowSizeCallback);
	glfwSetKeyCallback (windowPtr, keyCallback);
	glfwSetCursorPosCallback (windowPtr, cursorPosCallback);
	glfwSetMouseButtonCallback (windowPtr, mouseButtonCallback);
}

void initScene () {
	scenePtr = std::make_shared<Scene> ();
	//scenePtr->setBackgroundColor (glm::vec3 (0.1f, 0.5f, 0.95f));
	scenePtr->setBackgroundColor (glm::vec3 (0.0f, 0.0f, 0.0f));

	// Main object

	std::shared_ptr<Mesh> mainMeshPtr;
	try {
		mainMeshPtr = IO::loadOFFMesh (meshFilename);
	} catch (std::exception & e) {
		exitOnCriticalError (std::string ("[Error loading mesh]") + e.what ());
	}
	const auto & P = mainMeshPtr->vertexPositions();
	glm::vec3 minP = P[0];
	glm::vec3 maxP = P[0];
	for (const auto & x : P) {
		for (size_t j = 0; j < 3; ++j) {
			if (x[j] < minP[j])
				minP[j] = x[j];
			if (x[j] > maxP[j])
				maxP[j] = x[j];
		}
	}
	center = (minP+maxP)/2.f;
	float width = maxP[0] - minP[0];
	float height = maxP[1] - minP[1];
	float length = maxP[2] - minP[2];
	meshScale = glm::max (width, glm::max (height, length));

	auto mainMaterialPtr = std::make_shared<Material> (glm::vec3 (1.0, 0.85, 0.0f), 0.4, 0.6, 1.33);
	scenePtr->add (mainMeshPtr);
	scenePtr->add (mainMaterialPtr);
	scenePtr->assignMaterial (0, 0);

	// Adding a ground adapted to the loaded model

	std::shared_ptr<Mesh> groundMeshPtr = std::make_shared<Mesh> ();
	float extent = meshScale;
	glm::vec3 startP = center + glm::vec3 (-extent, -height/2.f, -extent);
	groundMeshPtr->vertexPositions().push_back (startP); 
	groundMeshPtr->vertexPositions().push_back (startP + glm::vec3 (0.f, 0.f, 2.f*extent));
	groundMeshPtr->vertexPositions().push_back (startP + glm::vec3 (2.f*extent, 0.f, 2.f*extent));
	groundMeshPtr->vertexPositions().push_back (startP + glm::vec3 (2.f*extent, 0.f, 0.f));
	groundMeshPtr->triangleIndices().push_back (glm::uvec3 (0, 1, 2));
	groundMeshPtr->triangleIndices().push_back (glm::uvec3 (0, 2, 3));
	groundMeshPtr->recomputePerVertexNormals ();
    auto groundMaterialPtr = std::make_shared<Material> (glm::vec3 (0.6, 0.6, 0.6f), 0.6, 0.0);
    scenePtr->add (groundMeshPtr);
    scenePtr->add (groundMaterialPtr);
	scenePtr->assignMaterial (1, 1);

	// Adding a wall adapted to the loaded model

	std::shared_ptr<Mesh> wallMeshPtr = std::make_shared<Mesh> ();
	startP = center + glm::vec3 (-extent, -height/2.f, -extent);
	wallMeshPtr->vertexPositions().push_back (startP); 
	wallMeshPtr->vertexPositions().push_back (startP + glm::vec3 (2.f*extent, 0.f, 0.f));
	wallMeshPtr->vertexPositions().push_back (startP + glm::vec3 (2.f*extent, 2.f*extent, 0.f));
	wallMeshPtr->vertexPositions().push_back (startP + glm::vec3 (0.f, 2.f*extent, 0.f));
	wallMeshPtr->vertexTextCoords().push_back (glm::vec2(0.f, 1.f));
	wallMeshPtr->vertexTextCoords().push_back (glm::vec2(1.f, 1.f));
	wallMeshPtr->vertexTextCoords().push_back (glm::vec2(1.f, 0.f));
	wallMeshPtr->vertexTextCoords().push_back (glm::vec2(0.f, 0.f));
	wallMeshPtr->triangleIndices().push_back (glm::uvec3 (0, 1, 2));
	wallMeshPtr->triangleIndices().push_back (glm::uvec3 (0, 2, 3));
	wallMeshPtr->recomputePerVertexNormals ();
    auto wallMaterialPtr = std::make_shared<Material> (glm::vec3 (0.9, 0.5, 0.3f), 0.3, 0.0);
    scenePtr->add (wallMeshPtr);
    scenePtr->add (wallMaterialPtr);
	scenePtr->assignMaterial (2, 2);

	// Light sources
	scenePtr->add (std::make_shared<LightSource> (normalize (glm::vec3(0.f, -1.f, -1.f)), glm::vec3(1.f, 1.f, 1.f), 0.5f)); // Key light
	scenePtr->add (std::make_shared<LightSource> (normalize (glm::vec3(-2.f, -0.5f, 0.f)), glm::vec3(0.2f, 0.6f, 1.f), 0.25f)); // Fill light
	scenePtr->add (std::make_shared<LightSource> (normalize (glm::vec3(2.f, -0.5f, 0.f)), glm::vec3(1.0f, 0.25f, 0.1f), 0.25f)); // Rim light

	// Camera
	int w, h;
	glfwGetWindowSize (windowPtr, &w, &h);
	auto cameraPtr = std::make_shared<Camera> ();
	cameraPtr->setAspectRatio (static_cast<float>(w) / static_cast<float>(h));
	glm::vec3 eye = center + glm::vec3 (0.0, 0.0, 1.5 * meshScale); // To make the object visible from the camera
	cameraPtr->setTranslation (eye);
	cameraPtr->setNear (0.1f * meshScale);
	cameraPtr->setFar (100.f * meshScale);
	scenePtr->set (cameraPtr);
}

void init () {
	initGLFW (); // Windowing system
	if (!gladLoadGLLoader ((GLADloadproc)glfwGetProcAddress)) // Load extensions for modern OpenGL
		exitOnCriticalError ("[Failed to initialize OpenGL context]");
	initScene (); // Actual scene to render
	rasterizerPtr = make_shared<Rasterizer> ();
	rasterizerPtr->init (basePath, scenePtr); // Mut be called before creating the scene, to generate an OpenGL context and allow mesh VBOs
	rayTracerPtr = make_shared<RayTracer> ();
	rayTracerPtr->init(scenePtr);
}

void clear () {
	glfwDestroyWindow (windowPtr);
	glfwTerminate ();
}

// The main rendering call
void render () {
	if (isDisplayRaytracing)
		rasterizerPtr->display (rayTracerPtr->image ());
	else
		rasterizerPtr->render (scenePtr);
}

// Update any accessible variable based on the current time
void update (float currentTime) {
	// Animate any entity of the program here
	static const float initialTime = currentTime;
	static float lastTime = 0.f;
	static unsigned int frameCount = 0;
	static float fpsTime = currentTime;
	static unsigned int FPS = 0;
	float elapsedTime = currentTime - initialTime;
	float dt = currentTime - lastTime;
	static float rotationAngle = 0.f;
	if (isAnimated) {
		rotationAngle += dt*glm::pi<float>()/(4.0);
	}
	scenePtr->mesh(0)->setRotation (glm::vec3 (0., rotationAngle, 0.));
	if (frameCount == 99) {
		float delai = (currentTime - fpsTime)/100;
		FPS = static_cast<unsigned int> (1.f/delai);
		frameCount = 0;
		fpsTime = currentTime;
	}
	std::string titleWithFPS = BASE_WINDOW_TITLE + " - " + std::to_string (FPS) + "FPS" 
							   + " - " + (isDisplayRaytracing ? "Raytracing" : "Rasterization");
	glfwSetWindowTitle (windowPtr, titleWithFPS.c_str ());
	lastTime = currentTime;
	frameCount++;
}

void usage (const char * command) {
	Console::print ("Usage : " + std::string(command) + " [<meshfile.off> [<material directory>]]");
	std::exit (EXIT_FAILURE);
}

void parseCommandLine (int argc, char ** argv) {
	if (argc > 2)
		usage (argv[0]);
	fs::path appPath = argv[0];
	basePath = appPath.parent_path().string(); 
	meshFilename = basePath + "/" + (argc >= 2 ? argv[1] : DEFAULT_MESH_FILENAME);
}

int main (int argc, char ** argv) {
	parseCommandLine (argc, argv);
	init (); 
	while (!glfwWindowShouldClose (windowPtr)) {
		update (static_cast<float> (glfwGetTime ()));
		render ();
		glfwSwapBuffers (windowPtr);
		glfwPollEvents ();
	}
	clear ();
	Console::print ("Quit");
	return EXIT_SUCCESS;
}

