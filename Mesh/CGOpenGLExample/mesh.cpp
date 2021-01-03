#define GLUT_DISABLE_ATEXIT_HACK
#include "mesh.h"
#include <windows.h>
#include <GL/glut.h>
#include "mesh.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#define PI 3.1415926


using namespace std;

void Mesh::leer(string path, int xxx)
{
	ifstream input_file(path);

	const int LINE_LENGTH = 100;

	n_vertices = 1;
	n_edges = 0;
	n_faces = 0;

	while (input_file) {
		string type;

		input_file >> type;

		if (type.length() == 1) {

			switch (type[0]) {

			case 'v': {
				double x, y, z;
				input_file >> x >> y >> z;
				Vertex* vert;
				if(xxx==0)
					vert = new Vertex(x, y, z, n_vertices);
				else 
					vert = new Vertex(x*xxx, y*xxx, z*xxx, n_vertices);

				this->vertices.push_back(vert);

				//cout << n_vertices << " vertice: " << x << "," << y << "," << z << endl;

				n_vertices++;
				
				break;
			}
			case 'f': {
				int a, b, c;
				input_file >> a >> b >> c;

				//cout << "cara: " << a << "," << b << "," << c << endl;

				HalfEdge* he1, * he2, * he3;

				int change = 0;
				map<pair<int, int>, HalfEdge*>::iterator mapIt2;
				mapIt2 = heMap.find(pair<int, int>(a, b));
				if (mapIt2 != heMap.end()) {
					change = a;
					a = b;
					b = change;
				}

				mapIt2 = heMap.find(pair<int, int>(b, c));
				if (mapIt2 != heMap.end()) {
					change = b;
					b = c;
					c = change;
				}

				mapIt2 = heMap.find(pair<int, int>(c, a));
				if (mapIt2 != heMap.end()) {
					change = c;
					c = a;
					a = change;
				}

				he1 = new HalfEdge();
				he2 = new HalfEdge();
				he3 = new HalfEdge();

				he1->vert = this->vertices[b-1];
				he2->vert = this->vertices[c-1];
				he3->vert = this->vertices[a-1];

				he1->next = he2;
				he2->next = he3;
				he3->next = he1;

				he1->prev = he3;
				he2->prev = he1;
				he3->prev = he2;

				vertices[a - 1]->halfEdge = he1;
				vertices[b - 1]->halfEdge = he2;
				vertices[c - 1]->halfEdge = he3;

				half_edges.push_back(he1);
				half_edges.push_back(he2);
				half_edges.push_back(he3);

				map<pair<int, int>, HalfEdge*>::iterator mapIt;
			
				pair<int, int> edgeIndex(a, b);
				heMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he1));

				edgeIndex = pair<int, int>(b, c);
				heMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he2));

				edgeIndex = pair<int, int>(c, a);
				heMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he3));


				mapIt = heMap.find(pair<int, int>(b, a));
				if (mapIt != heMap.end()) {
						heMap[pair<int, int>(a, b)]->twin = mapIt->second;
						he1->twin = mapIt->second;
						mapIt->second->twin = he1;
						Edge* edge = new Edge(he1, he1->twin);
						edges.push_back(edge);
						//cout << "EDGE AÑADIDA" << endl;


				}

				mapIt = heMap.find(pair<int, int>(c, b));
				if (mapIt != heMap.end()) {
						heMap[pair<int, int>(b, c)]->twin = mapIt->second;
						he2->twin = mapIt->second;
						mapIt->second->twin = he2;
						Edge* edge = new Edge(he2, he2->twin);
						edges.push_back(edge);
						//cout << "EDGE AÑADIDA" << endl;
				}

				mapIt = heMap.find(pair<int, int>(a, c));
				if (mapIt != heMap.end()) {
						heMap[pair<int, int>(c, a)]->twin = mapIt->second;
						he3->twin = mapIt->second;
						mapIt->second->twin = he3;
						Edge* edge = new Edge(he3, he3->twin);
						edges.push_back(edge);
						//cout << "EDGE AÑADIDA" << endl;
				}


				Face* f = new Face(he1, he2, he3, n_faces);
				this->faces.push_back(f);


				edgeIndex = pair<int, int>(a, b);
				heMap[edgeIndex]->face = f;
				edgeIndex = pair<int, int>(b, c);
				heMap[edgeIndex]->face = f;
				edgeIndex = pair<int, int>(c, a);
				heMap[edgeIndex]->face = f;

				n_faces++;
			}
			}
		}
	}
}

void Mesh::Draw(bool v) {

	if(v)
	for (int i = 0; i < this->faces.size(); i++) {
		glLineWidth(3);
		glBegin(GL_LINES);
		glColor3d(0.49, 0.44, 0.05);
		glVertex3d(this->faces[i]->Edges[0]->vert->pos[0], this->faces[i]->Edges[0]->vert->pos[1], this->faces[i]->Edges[0]->vert->pos[2]);
		glVertex3d(this->faces[i]->Edges[1]->vert->pos[0], this->faces[i]->Edges[1]->vert->pos[1], this->faces[i]->Edges[1]->vert->pos[2]);

		glVertex3d(this->faces[i]->Edges[1]->vert->pos[0], this->faces[i]->Edges[1]->vert->pos[1], this->faces[i]->Edges[1]->vert->pos[2]);
		glVertex3d(this->faces[i]->Edges[2]->vert->pos[0], this->faces[i]->Edges[2]->vert->pos[1], this->faces[i]->Edges[2]->vert->pos[2]);

		glVertex3d(this->faces[i]->Edges[2]->vert->pos[0], this->faces[i]->Edges[2]->vert->pos[1], this->faces[i]->Edges[2]->vert->pos[2]);
		glVertex3d(this->faces[i]->Edges[0]->vert->pos[0], this->faces[i]->Edges[0]->vert->pos[1], this->faces[i]->Edges[0]->vert->pos[2]);

		glEnd();
	}
	
	else
	for (int i = 0; i < this->faces.size(); i++) {

		glBegin(GL_TRIANGLES);
		glColor3d(1, 0.88, 0.07);
		glVertex3d(this->faces[i]->Edges[0]->vert->pos[0], this->faces[i]->Edges[0]->vert->pos[1], this->faces[i]->Edges[0]->vert->pos[2]);
		glVertex3d(this->faces[i]->Edges[1]->vert->pos[0], this->faces[i]->Edges[1]->vert->pos[1], this->faces[i]->Edges[1]->vert->pos[2]);
		glVertex3d(this->faces[i]->Edges[2]->vert->pos[0], this->faces[i]->Edges[2]->vert->pos[1], this->faces[i]->Edges[2]->vert->pos[2]);
		glEnd();

	}

}



pair<double, int> Mesh::numVecinos(Vertex* v) {
	pair<double, int> res;
	double cont = 0;
	HalfEdge* inicio = v->halfEdge;
	HalfEdge* c = v->halfEdge;

	do {
		c = c->twin->next;
		cont++;
	} while (c != inicio);

	if (cont > 3) {
		res.first=(3/(8 * cont));
		res.second = cont;
		return res;
	}

	else if (cont == 3) {
		res.first = (0.1875);
		res.second = cont;
		return res;
	}
}

Vertex* Mesh::sumVecinos(Vertex* v,double n) {
	HalfEdge* c = v->halfEdge;
	Vertex* res = c->vert;
	Vertex* resultado;

	do {
		c = c->twin->next;
		res = res->suma(c->vert,0);
	} while (c != v->halfEdge);

	resultado = new Vertex(res->pos[0] * n, res->pos[1] * n, res->pos[2] * n,0);
	return resultado;
}


void Mesh::loop() {
	vector<Vertex*> nuevosVertices;
	vector<Vertex*> verticesOriginales;
	vector<Vertex*> nuevosVerticesMedios;

	vector<Face*> nuevasFaces;
	vector<Edge*> edges_aux;
	vector<Face*> faces_aux;
	vector<HalfEdge*> half_edges_aux;

	map<pair<int, int>, HalfEdge*>  hMap;


	Vertex* vi;
	Vertex* vj;
	Vertex* vk;
	Vertex* vl;
	Vertex* xn;
	Vertex* aux1;
	Vertex* aux2;


	for (int i = 0; i < vertices.size(); i++) {
		double a, b, c;
		double a1 = 0, b1 = 0, c1 = 0;
		a = vertices[i]->pos[0];
		b = vertices[i]->pos[1];
		c = vertices[i]->pos[2];

		double N = 0;
		HalfEdge* ini = vertices[i]->halfEdge;
		do {
			a1 += ini->vert->pos[0];
			b1 += ini->vert->pos[1];
			c1 += ini->vert->pos[2];
			ini = ini->twin->next;
			N++;

		} while (ini != vertices[i]->halfEdge);

		double w;
		w = (0.625 - pow((0.375 + (0.25 * cos(2 * PI / N))), 2)) / N;

		a = a * (1 - N * w);
		a1 = a1 * w;

		b = b * (1 - N * w);
		b1 = b1 * w;

		c = c * (1 - N * w);
		c1 = c1 * w;

		a = a + a1;
		b = b + b1;
		c = c + c1;

		Vertex* nV = new Vertex(a, b, c, vertices[i]->id);

		verticesOriginales.push_back(nV);
	}

	for (int i = 0; i < edges.size(); i++) {
		double a, b, c, a1, b1, c1;
		vi = edges[i]->parEdges[0]->vert;
		vj = edges[i]->parEdges[1]->vert;
		vk = edges[i]->parEdges[0]->next->vert;
		vl = edges[i]->parEdges[1]->next->vert;

		a = (vi->pos[0] + vj->pos[0]) * 0.375;
		b = (vi->pos[1] + vj->pos[1]) * 0.375;
		c = (vi->pos[2] + vj->pos[2]) * 0.375;

		a1 = (vk->pos[0] + vl->pos[0]) * 0.125;
		b1 = (vk->pos[1] + vl->pos[1]) * 0.125;
		c1 = (vk->pos[2] + vl->pos[2]) * 0.125;

		a += a1;
		b += b1;
		c += c1;

		Vertex* nV = new Vertex(a, b, c, i + 1 + vertices.size());

		edges[i]->parEdges[0]->mid = nV;
		edges[i]->parEdges[1]->mid = nV;
		edges[i]->mid = nV;

		nuevosVerticesMedios.push_back(nV);
	}

	for (int i = 0; i < vertices.size(); i++) {
		vertices[i]->pos[0] = verticesOriginales[i]->pos[0];
		vertices[i]->pos[1] = verticesOriginales[i]->pos[1];
		vertices[i]->pos[2] = verticesOriginales[i]->pos[2];
		nuevosVertices.push_back(vertices[i]);
	}

	for (int i = 0; i < nuevosVerticesMedios.size(); i++) {
		nuevosVertices.push_back(nuevosVerticesMedios[i]);
	}


	int n_f = 0;
	for (int i = 0; i < faces.size(); i++) {
		for (int e = 0; e < 3; e++) {

			HalfEdge* ed = faces[i]->Edges[e];
			int a = ed->mid->id;
			int b = ed->vert->id;
			int c = ed->next->mid->id;
			int change = 0;

			map<pair<int, int>, HalfEdge*>::iterator mapIt;
			mapIt = hMap.find(pair<int, int>(a, b));
			if (mapIt != hMap.end()) {
				change = a;
				a = b;
				b = change;
			}

			mapIt = hMap.find(pair<int, int>(b, c));
			if (mapIt != hMap.end()) {
				change = b;
				b = c;
				c = change;
			}

			mapIt = hMap.find(pair<int, int>(c, a));
			if (mapIt != hMap.end()) {
				change = c;
				c = a;
				a = change;
			}

			HalfEdge* he1;
			HalfEdge* he2;
			HalfEdge* he3;

			he1 = new HalfEdge();
			he2 = new HalfEdge();
			he3 = new HalfEdge();

			he1->vert = ed->vert;
			he2->vert = ed->next->mid;
			he3->vert = ed->mid;

			he1->next = he2;
			he2->next = he3;
			he3->next = he1;

			he1->prev = he3;
			he2->prev = he1;
			he3->prev = he2;

			ed->mid->halfEdge = he1;
			ed->vert->halfEdge = he2;
			ed->next->mid->halfEdge = he3;

			half_edges_aux.push_back(he1);
			half_edges_aux.push_back(he2);
			half_edges_aux.push_back(he3);

			pair<int, int> edgeIndex(a, b);
			hMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he1));

			edgeIndex = pair<int, int>(b, c);
			hMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he2));

			edgeIndex = pair<int, int>(c, a);
			hMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he3));

			mapIt = hMap.find(pair<int, int>(b, a));
			if (mapIt != hMap.end()) {
				hMap[pair<int, int>(a, b)]->twin = mapIt->second;
				he1->twin = mapIt->second;
				mapIt->second->twin = he1;
				Edge* edge = new Edge(he1, he1->twin);
				edges_aux.push_back(edge);
			}

			mapIt = hMap.find(pair<int, int>(c, b));
			if (mapIt != hMap.end()) {
				hMap[pair<int, int>(b, c)]->twin = mapIt->second;
				he2->twin = mapIt->second;
				mapIt->second->twin = he2;
				Edge* edge = new Edge(he2, he2->twin);
				edges_aux.push_back(edge);
			}

			mapIt = hMap.find(pair<int, int>(a, c));
			if (mapIt != hMap.end()) {
				hMap[pair<int, int>(c, a)]->twin = mapIt->second;
				he3->twin = mapIt->second;
				mapIt->second->twin = he3;
				Edge* edge = new Edge(he3, he3->twin);
				edges_aux.push_back(edge);
			}

			Face* f = new Face(he1, he2, he3, n_f);
			faces_aux.push_back(f);


			edgeIndex = pair<int, int>(a, b);
			hMap[edgeIndex]->face = f;
			edgeIndex = pair<int, int>(b, c);
			hMap[edgeIndex]->face = f;
			edgeIndex = pair<int, int>(c, a);
			hMap[edgeIndex]->face = f;

			n_f++;

		}

	}

	for (int i = 0; i < faces.size(); i++) {
		HalfEdge* ed = faces[i]->Edges[0];
		int a = ed->mid->id;
		int b = ed->next->mid->id;
		int c = ed->next->next->mid->id;
		int change = 0;

		map<pair<int, int>, HalfEdge*>::iterator mapIt;
		mapIt = hMap.find(pair<int, int>(a, b));
		if (mapIt != hMap.end()) {
			change = a;
			a = b;
			b = change;
		}

		mapIt = hMap.find(pair<int, int>(b, c));
		if (mapIt != hMap.end()) {
			change = b;
			b = c;
			c = change;
		}

		mapIt = hMap.find(pair<int, int>(c, a));
		if (mapIt != hMap.end()) {
			change = c;
			c = a;
			a = change;
		}

		HalfEdge* he1;
		HalfEdge* he2;
		HalfEdge* he3;

		he1 = new HalfEdge();
		he2 = new HalfEdge();
		he3 = new HalfEdge();

		he1->vert = ed->next->mid; //b
		he2->vert = ed->next->next->mid; //c
		he3->vert = ed->mid; //a

		he1->next = he2;
		he2->next = he3;
		he3->next = he1;

		he1->prev = he3;
		he2->prev = he1;
		he3->prev = he2;

		ed->mid->halfEdge = he1;
		ed->next->mid->halfEdge = he2;
		ed->next->next->mid->halfEdge = he3;

		half_edges_aux.push_back(he1);
		half_edges_aux.push_back(he2);
		half_edges_aux.push_back(he3);

		pair<int, int> edgeIndex(a, b);
		hMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he1));

		edgeIndex = pair<int, int>(b, c);
		hMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he2));

		edgeIndex = pair<int, int>(c, a);
		hMap.insert(map<pair<int, int>, HalfEdge*>::value_type(edgeIndex, he3));

		mapIt = hMap.find(pair<int, int>(b, a));
		if (mapIt != hMap.end()) {
			hMap[pair<int, int>(a, b)]->twin = mapIt->second;
			he1->twin = mapIt->second;
			mapIt->second->twin = he1;
			Edge* edge = new Edge(he1, he1->twin);
			edges_aux.push_back(edge);
		}

		mapIt = hMap.find(pair<int, int>(c, b));
		if (mapIt != hMap.end()) {
			hMap[pair<int, int>(b, c)]->twin = mapIt->second;
			he2->twin = mapIt->second;
			mapIt->second->twin = he2;
			Edge* edge = new Edge(he2, he2->twin);
			edges_aux.push_back(edge);
		}

		mapIt = hMap.find(pair<int, int>(a, c));
		if (mapIt != hMap.end()) {
			hMap[pair<int, int>(c, a)]->twin = mapIt->second;
			he3->twin = mapIt->second;
			mapIt->second->twin = he3;
			Edge* edge = new Edge(he3, he3->twin);
			edges_aux.push_back(edge);
		}

		Face* f = new Face(he1, he2, he3, n_f);
		faces_aux.push_back(f);

		edgeIndex = pair<int, int>(a, b);
		hMap[edgeIndex]->face = f;
		edgeIndex = pair<int, int>(b, c);
		hMap[edgeIndex]->face = f;
		edgeIndex = pair<int, int>(c, a);
		hMap[edgeIndex]->face = f;

		n_f++;
	}

	vertices = nuevosVertices;
	faces = faces_aux;
	edges = edges_aux;
	half_edges = half_edges_aux;

}


Vertex* Mesh::eqPlane(Face* f) {
	double x1 = f->Edges[0]->vert->pos[0];
	double y1 = f->Edges[0]->vert->pos[1];
	double z1 = f->Edges[0]->vert->pos[2];
	double x2 = f->Edges[1]->vert->pos[0];
	double y2 = f->Edges[1]->vert->pos[1];
	double z2 = f->Edges[1]->vert->pos[2];
	double x3 = f->Edges[2]->vert->pos[0];
	double y3 = f->Edges[2]->vert->pos[1];
	double z3 = f->Edges[2]->vert->pos[2];

	double a1 = x2 - x1;
	double b1 = y2 - y1;
	double c1 = z2 - z1;
	double a2 = x3 - x1;
	double b2 = y3 - y1;
	double c2 = z3 - z1;
	double a = b1 * c2 - b2 * c1;
	double b = a2 * c1 - a1 * c2;
	double c = a1 * b2 - b1 * a2;
	double d = (-a * x1 - b * y1 - c * z1);

	Vertex* p = new Vertex(a, b, c, d);
	return p;
}


double Mesh::angulo(Face* f1, Face* f2) {
	Vertex* v = eqPlane(f1);
	Vertex* u = eqPlane(f2);

	double num = abs((v->pos[0] * u->pos[0]) + (v->pos[1] * u->pos[1]) + (v->pos[2] * u->pos[2]));
	double den = sqrt(pow(v->pos[0], 2) + pow(v->pos[1], 2) + pow(v->pos[2], 2)) * sqrt(pow(u->pos[0], 2) + pow(u->pos[1], 2) + pow(u->pos[2], 2));
	double res = num / den;
	return res;

}

void Mesh::colapse() {
	Face* f1 = faces[0];
	Face* f2 = faces[0];
	double minimo = 1000000;
	double aux;

	for (int i = 0; i < faces.size() - 1; i++) {
		aux = angulo(faces[i], faces[i + 1]);
		if (aux < minimo) {
			minimo = aux;
			Face* f1 = faces[i];
			f2 = faces[i + 1];
		}
	}

	Edge* ed = edges[0];
	for (int i = 0; i < 3; i++) {
		for (int e = 0; e < 3; e++) {
			if (f1->Edges[i]->twin == f2->Edges[e]) {
				ed->parEdges[0] = f2->Edges[e];
				ed->parEdges[0] = f2->Edges[e]->twin;
				i = 1000;
				break;
			}
		}
	}

	Vertex* v1 = ed->parEdges[0]->vert;
	Vertex* v2 = ed->parEdges[0]->twin->vert;

	vector<HalfEdge*> eliminadas;
	eliminadas.push_back(ed->parEdges[0]);
	eliminadas.push_back(ed->parEdges[0]->next);
	eliminadas.push_back(ed->parEdges[0]->next->next);
	eliminadas.push_back(ed->parEdges[0]->twin);
	eliminadas.push_back(ed->parEdges[0]->twin->next);
	eliminadas.push_back(ed->parEdges[0]->twin->next->next);

	vector<Vertex*> eliminados;
	eliminados.push_back(ed->parEdges[0]->vert);
	eliminados.push_back(ed->parEdges[0]->next->vert);
	eliminados.push_back(ed->parEdges[0]->twin->vert);
	eliminados.push_back(ed->parEdges[0]->twin->next->vert);



	ed->parEdges[0]->prev->twin->twin = ed->parEdges[0]->next->twin;	
	ed->parEdges[0]->twin->next->twin->twin = ed->parEdges[0]->twin->prev->twin;

	ed->parEdges[0]->twin->vert->pos[0] = ed->parEdges[0]->vert->pos[0];
	ed->parEdges[0]->twin->vert->pos[1] = ed->parEdges[0]->vert->pos[1];
	ed->parEdges[0]->twin->vert->pos[2] = ed->parEdges[0]->vert->pos[2];



	for (int i = 0; i < edges.size(); i++) {
		if (edges[i] == ed) {
			edges.erase(edges.begin() + i);
			cout << "edge borrada" << endl;
		}
	}

	//for (int i = 0; i < half_edges.size(); i++) {
	//	if (half_edges[i]->vert == v1) {
	//		half_edges[i]->vert = v2;
	//		cout << " vertice igualado " << endl;
	//	}
	//}

	//for (int i = 0; i < vertices.size(); i++) {
	//	if (vertices[i] == v1) {
	//		vertices.erase(vertices.begin() + i);
	//		cout << "vertice borrado" << endl;
	//	}
	//}

	//for (int i = 0; i < half_edges.size(); i++) {
	//	for (int e = 0; e < eliminadas.size(); e++) {
	//		if (half_edges[i] == eliminadas[e]) {
	//			half_edges.erase(half_edges.begin() + i);
	//			cout << "he borrada" << endl;
	//		}

	//	}
	//}

	//for (int i = 0; i < faces.size(); i++) {
	//	if (faces[i] == f1 || faces[i] == f2) {
	//		faces.erase(faces.begin() + i);
	//		cout << "face borrada" << endl;
	//	}
	//}

}
 