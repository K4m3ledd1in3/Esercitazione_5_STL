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
	cout << "f"<<endl;
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
        mesh.Marker0d.push_back(marker);
        //cout<< mesh.Marker0d[id] << endl;
        //cout << id << endl;
        if( marker!=0 && ( eps < x && x < 1-eps ) && (eps < y && y < 1-eps) )
		{
        	cout << "The cell " << id << " is as the cell that hasn't the marker stored correctly." << endl;
        	flag_marker = false;
		}
    }
    if(flag_marker){
    	cout << "All cells 0d have the marker stored correctly." << endl;
	}else{
		cout << "The markers 0d aren't stored correctly" << endl;
	}
    return true;
}
bool ImportCell1Ds(PolygonalMesh& mesh)
{
	bool flag_marker = true;
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
	double eps = 2.2*10e-16;
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
        unsigned int p0= mesh.Cell1DsExtrema(0, id);
        unsigned int p1= mesh.Cell1DsExtrema(1, id);
        unsigned int m0 = mesh.Marker0d[p0];
        unsigned int m1 = mesh.Marker0d[p1];
        double x0 = mesh.Cell0DsCoordinates(0,p0);
        double y0 = mesh.Cell0DsCoordinates(1,p0);
        double x1 = mesh.Cell0DsCoordinates(0,p1);
        double y1 = mesh.Cell0DsCoordinates(1,p1);
        if(m0!=0 && m1!=0){
        	switch((int) (m0!=m1)){
        		case 1:
        			if((( x0>eps && x0<1-eps)||(y0>eps && y0<1-eps)) ||
					   (( x1>eps && x1<1-eps )&&( y1>eps && y1<1-eps)) ){
        				flag_marker=false;
					}
        		break;
        		default:
        			flag_marker=true;
        			break;
			}
		}
    }
	if(flag_marker){
    	cout << "All cells 1d have the marker stored correctly." << endl;
	} else{
		cout << "The markers 1d aren't stored correctly" << endl;
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
	cout << "\t\tArea" << endl;
	for(auto& c: mesh.Cell2DsId)
		{
			Area = 0.0;
			vector<double> buff;
			double X,Y;
			buff.reserve( 2*mesh.Cell2DsVertices[c].size() );
			for(size_t j = 0; j<mesh.Cell2DsVertices[c].size(); j++){
				unsigned int k = mesh.Cell2DsVertices[c][j];
				double x_s = mesh.Cell0DsCoordinates(0,k); //vertex_x_coordinate.
				double y_s = mesh.Cell0DsCoordinates(1,k); //vertex_y_coordinate.
				buff[2*j] = x_s;
				buff[2*j+1] = y_s;
				//cout << c << " > (x,y)=" <<"(" <<buff[2*j] << "," <<buff[2*j+1] <<")" << endl;
				X+=(x_s) / ( (double)(mesh.Cell2DsVertices[c].size()) ); //barycenter X by arithmetic mean on uniform density. .
				Y+=(y_s) / ( (double)(mesh.Cell2DsVertices[c].size()) ); //barycenter by arithmentic mean on uniform density.
				}
			for(size_t i = 0; i<mesh.Cell2DsVertices[c].size(); i++)
			{
				double x_i = buff[(2*i)%(2*mesh.Cell2DsVertices[c].size())]; //%mod(#buff) needs in order to slide from last vertex up to the first, 
				double y_i1 = buff[(2*i+3)%(2*mesh.Cell2DsVertices[c].size())];
				double x_i1 = buff[(2*i+2)%(2*mesh.Cell2DsVertices[c].size())]; 
				double y_i = buff[(2*i+1)%(2*mesh.Cell2DsVertices[c].size())];
				x_i-=X;
				y_i-=Y;
				x_i1-=X;  //Computation of the area from the barycenter frame, it wouldn't change anything 'cause the area doesn't affected by the reference frame.
				y_i1-=Y;
				Area+=(double) (x_i*y_i1-x_i1*y_i)*(0.5);
			}
			
			Area=(double)abs(Area);
			
			cout << " (ID) c: " << c << "-> " << Area << endl;
			TotArea+=Area;
			if(Area<tau){
				cout << "The polygon at cell " << c << " has zero area" << endl;
				flag = false;
			}
		}
		if(flag)
			cout << "None of the cell is such as that has null area."<<endl;
		cout << "Tot_Area of Mesh: " << TotArea << endl; 
}
}
	