#pragma once

#include "../../geomodcore.h"
#include "simpletmesh.h"
#include "../grid.h"

namespace geomodcore {

	class TMeshSamples
	{
	public:
		static SimpleTMesh<>* getSampleTMesh();
		static SimpleTMesh<>* getSampleTMesh2();
		//static SimpleTMesh<>* getSampleTMeshMoscow();
	};
}