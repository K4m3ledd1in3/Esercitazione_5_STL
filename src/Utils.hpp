
#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolygonalLibrary
{
/// Import the triangular mesh and test if the mesh is correct
bool ImportMesh(PolygonalMesh& mesh);
/// Import the Cell0D properties from Cell0Ds.csv file
bool ImportCell0Ds(PolygonalMesh& mesh);
/// Import the Cell1D properties from Cell1Ds.csv file
bool ImportCell1Ds(PolygonalMesh& mesh);
/// Import the Cell2D properties from Cell2Ds.csv file
bool ImportCell2Ds(PolygonalMesh& mesh);

void NonZeroLength(PolygonalMesh& mesh);
void NonZero_Area(PolygonalMesh& mesh);
}
