#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
namespace PolygonalLibrary
{
bool ImportMesh(PolygonalMesh& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}
// ***************************************************************************
bool ImportCell0Ds(PolygonalMesh& mesh)
{
	char separator=';';
    ifstream file("./Cell0Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);
	bool flag_marker = true;
	double x,y;
	double eps=2.2*10e-16;
    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2d coord;

        converter >>  id >> separator >> marker 
		>> separator
		>> mesh.Cell0DsCoordinates(0, id) >> separator
		>> mesh.Cell0DsCoordinates(1, id);
              /// Memorizza i marker
    	x = mesh.Cell0DsCoordinates(0, id);
    	y= mesh.Cell0DsCoordinates(1, id);
        mesh.Cell0DsId.push_back(id);
        if(marker!=0 &&(eps < x && x < 1-eps ) && (eps < y && y < 1-eps) )
		{
        	cout << "The cell " << id << " is as the cell that hasn't the marker stored correctly." << endl;
        	flag_marker = false;
		}
    }
    if(flag_marker){
    	cout << "All cells have the marker stored correctly" << endl;
	}
    
    /*for(size_t i=0; i<mesh.NumCell0Ds; i++)
 		cout << "(" << mesh.Cell0DsCoordinates(0,i) << ","<< mesh.Cell0DsCoordinates(1,i) << ")"<<endl;*/
 	
    return true;
}
// ***************************************************************************

bool ImportCell1Ds(PolygonalMesh& mesh)
{
	char separator=';';
    ifstream file("./Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2i vertices;

        converter >>  id >> separator 
		>> marker >> separator >> mesh.Cell1DsExtrema(0, id) >>
		separator >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);
        
    }
    /*cout << "1D" << endl;
    for(size_t i =0; i<mesh.NumCell1Ds; i++){
    	cout << mesh.Cell1DsExtrema(0,i) << "," << mesh.Cell1DsExtrema(1,i) << endl;
	}*/
    
    
    return true;
}
// ***************************************************************************
void NonZeroLength(PolygonalMesh& mesh)
{
	bool flag = true;
	for(size_t i = 0; i<mesh.NumCell1Ds; i++){
		unsigned int cell_0 = mesh.Cell1DsExtrema(0, i);
		unsigned int cell_1 = mesh.Cell1DsExtrema(1, i);
	     VectorXd v1(2);
    	 VectorXd v2(2);
    	 v1 << mesh.Cell0DsCoordinates(0,cell_0) , mesh.Cell0DsCoordinates(1,cell_0);
    	 v2 << mesh.Cell0DsCoordinates(0,cell_1) , mesh.Cell0DsCoordinates(1,cell_1);
	     VectorXd v3 = v1 - v2;
    	 double norm = v3.norm();
    	 if(norm < 2.2*10e-16){
    	 	flag = false;
    	 	cout << "Cell " << i << " has zero length" << endl;
		 }
	}
	if(flag)
		cout << "No cell has zero length." << endl;
}

// ***************************************************************************
bool ImportCell2Ds(PolygonalMesh& mesh)
{
	char separator=';';
    ifstream file;
    file.open("./Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();
    // remove header
    listLines.pop_front();
    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);
        unsigned int id;
        unsigned int num_vertices;
        unsigned int num_edges;
        unsigned int marker;
        
        vector<unsigned int> vertices;
        vector<unsigned int> edges;
        
		converter >>  id 
		>> separator >> marker
		>> separator >> num_vertices;
		
		//cout << id << " " << marker << " " << num_vertices << endl;
		vertices.reserve(num_vertices);  
		//cout << "ver: ";      
     for(unsigned int i = 0; i < num_vertices; i++)
			{
            unsigned int v;
            converter >> separator >> v;
    		vertices.push_back(v); 
			//cout << vertices[i] << " ";
			}
			
        converter >> separator >> num_edges;
        
        //cout << "num " <<  num_edges <<" ";
        edges.reserve(num_edges);
         
        for(unsigned int i = 0; i < num_edges; i++)
			{
            unsigned int e;
            converter >> separator >> e;
    		edges.push_back(e); 
    		//cout << edges[i] << " ";
			}
		//	cout << endl;
        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);
    }
    return true;
}
void NonZero_Area(PolygonalMesh& mesh){
	
	//mappa[id_vertex]: output: archi adiacenti.
	for(auto& c : mesh.Cell2DsVertices){
		for(auto& k : c){
			cout << k << " ";
		}cout << endl;
	}
}
}
