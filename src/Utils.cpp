#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unordered_map>
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
    return true;
}
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

    
    
    return true;
}
 
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
 
    listLines.pop_front();
    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
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
		
  		mesh.Cell2DsVertices[id].reserve(num_vertices);
		vertices.reserve(num_vertices);  
   
     for(unsigned int i = 0; i < num_vertices; i++)
			{
            unsigned int v;
            converter >> separator >> v;
    		vertices.push_back(v); 
			}
			
        converter >> separator >> num_edges;
        edges.reserve(num_edges);
        mesh.Cell2DsVertices[id].reserve(num_vertices);
        for(unsigned int i = 0; i < num_edges; i++)
			{
            unsigned int e;
            converter >> separator >> e;
    		edges.push_back(e); 
			}
        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsVertices[id] = vertices;
        mesh.Cell2DsEdges[id] = edges;
    }
    return true;
}

void NonZero_Area(PolygonalMesh& mesh){
	double tau = 2.2*10e-16*2.2*10e-16;
	double TotArea = 0.0;
	double Area = 0.0;
	bool flag = true;
	for(auto& c: mesh.Cell2DsId)
		{
			Area = 0.0;
			vector<double> buff;
			double X,Y;
			buff.reserve(2*mesh.Cell2DsEdges[c].size());
			for(size_t j = 0; j<mesh.Cell2DsEdges[c].size(); j++){
				unsigned int k = mesh.Cell2DsEdges[c][j];
				unsigned int s = mesh.Cell1DsExtrema(0,k);
				double x_s = mesh.Cell0DsCoordinates(0,s);
				double y_s = mesh.Cell0DsCoordinates(1,s);
				buff[2*j] = x_s;
				buff[2*j+1] = y_s;
				X+=(x_s)/(double)(mesh.Cell2DsEdges[c].size());
				Y+=(y_s)/(double)(mesh.Cell2DsEdges[c].size());
				}
		 			
			for(size_t i = 0; i<mesh.Cell2DsEdges[c].size(); i++)
			{
				double x_i = buff[(2*i)%(mesh.Cell2DsEdges[c].size())];
				double y_i1 = buff[(2*i+3)%(mesh.Cell2DsEdges[c].size())];
				double x_i1 = buff[(2*i+2)%(mesh.Cell2DsEdges[c].size())]; 
				double y_i = buff[(2*i+1)%(mesh.Cell2DsEdges[c].size())];
				x_i-=X;
				y_i-=Y;
				x_i1-=X;
				y_i1-=Y;
				Area+=(double) (x_i*y_i1-x_i1*y_i)*(0.5);				
			}
			Area=(double)abs(Area);
			TotArea+=Area;
			if(Area<tau){
				cout << "The polygon at cell " << c << " has zero area" << endl;
				flag = false;
			}
		 
			
		}
		if(flag)
			cout << "None of the cell is such as that has null area."<<endl;
		cout << "Tot_Area " << TotArea << endl; 

}

}
	/*	for(auto& k: mesh.Cell2DsEdges[c]) 
			{
				unsigned int p = mesh.Cell1DsExtrema(1,k);
				unsigned int s = mesh.Cell1DsExtrema(0,k);
				double x_p = mesh.Cell0DsCoordinates(0,p);
				double y_p = mesh.Cell0DsCoordinates(1,p);
				buff.push_back(mesh.Cell0DsCoordinates(0,s));
				buff.push_back(mesh.Cell0DsCoordinates(1,s));
				buff.push_back(x_p);
				buff.push_back(y_p);
				for(size_t i = 0; i<mesh.Cell2DsEdges[c].size(); i++){
						unsigned int j = mesh.Cell2DsEdges[c][i]; //restituisce arco
						unsigned int w = mesh.Cell1DsExtrema(0,j); //estremo iniziale arco
						double x_w = mesh.Cell0DsCoordinates(0,w);
						double y_w = mesh.Cell0DsCoordinates(1,w);
						if(abs(x_w-x_p)<2.2*10e-16 && abs(y_w-y_p)<2.2*10e-16){
							unsigned int h = mesh.Cell1DsExtrema(1,j);
							buff.push_back(x_w);
							buff.push_back(y_w);
							buff.push_back(mesh.Cell0DsCoordinates(0,h));
							buff.push_back(mesh.Cell0DsCoordinates(1,h));
							i = mesh.Cell2DsEdges[c].size();
						}
						
				}
				
				
			}*/