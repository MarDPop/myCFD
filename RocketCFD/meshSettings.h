#pragma once

enum MESH_TYPE
{
	DEFAULT, // structured
	TRIMMED, // structured
	POLYHEDRAL,
	TRIANGULAR
};

enum MESH_DIMENSION
{
	TWOD,
	AXISYMMETRIC,
	THREED
};

struct MeshSettings
{
	MESH_TYPE meshType;
	MESH_DIMENSION dimension;
	bool prismLayers;
	bool viscous;
	bool structured;
};

