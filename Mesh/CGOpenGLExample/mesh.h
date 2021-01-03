#pragma once 
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace std;

struct HalfEdge;
struct Vertex;
struct Face;
struct Edge;


struct HalfEdge {
	Vertex* vert;
	Face* face;
	Vertex* mid;

	HalfEdge* next;
	HalfEdge* prev;
	HalfEdge* twin;
};

struct Edge {
	HalfEdge *parEdges[2];
	Vertex* mid;

	Edge(HalfEdge* half1) {
		parEdges[0] = half1;
	}
	Edge(HalfEdge* half1, HalfEdge* half2)
	{
		parEdges[0] = half1;
		parEdges[1] = half2;
	}

	Vertex* prom(int id);
	
};

struct Vertex {
	double pos[3];
	int id;

	HalfEdge* halfEdge;

	Vertex(double x, double y, double z, int _id) {
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
		id = _id;
	}

	Vertex * suma(Vertex *v, double val) {
		int a, b, c;
		a = pos[0];
		b = pos[1];
		c = pos[2];

		if(val==0.0)
			return new Vertex(a + v->pos[0], b + v->pos[1], c + v->pos[2], 0);
		else
			return new Vertex((a + v->pos[0])*val, (b + v->pos[1])*val, (c + v->pos[2])*val, 0);

	}
};

struct Face {
	HalfEdge* Edges[3];
	//HalfEdge* edge;
	int id;




	Face(Vertex* v0, Vertex* v1, Vertex* v2, int _id)
	{
		Edges[0]->vert = v0;
		Edges[1]->vert = v1;
		Edges[2]->vert = v2;
		id = _id;
	}

	Face(HalfEdge* he1, HalfEdge* he2, HalfEdge* he3, int _id)
	{
		Edges[0] = he1;
		Edges[1] = he2;
		Edges[2] = he3;
		id = _id;
	}
};

class Mesh {

public:
	vector<Vertex * > vertices;
	vector<Vertex* > verticesAux;
	vector<Face *> faces;
	vector<Edge*>  edges;
	vector<HalfEdge*>  half_edges;

	map<pair<int, int>, HalfEdge*> heMap;


	int n_vertices, n_faces, n_edges;

	void Neighbors(int a);
	void Paint(int id);
	void Draw(bool v);
	void leer(string path, int val);
	vector<Vertex*> division();
	void loop();
	void colapse();
	pair<double,int> numVecinos(Vertex* v);
	Vertex* sumVecinos(Vertex* v,double n);
	vector<Vertex*>Vecinos(Vertex* v);
	double angulo(Face* v, Face* u);
	Vertex* eqPlane(Face* f);
};
