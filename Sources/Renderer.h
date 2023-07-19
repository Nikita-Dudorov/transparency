#pragma once

#include "Scene.h"

class Renderer {
public:
	virtual void render (std::shared_ptr<Scene> scenePtr) = 0;
};