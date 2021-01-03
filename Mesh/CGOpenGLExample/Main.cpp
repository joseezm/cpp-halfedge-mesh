#define GLUT_DISABLE_ATEXIT_HACK
#include <windows.h>
#include<GL/glut.h>
#include "mesh.h"
#define KEY_ESC 27
#define KEY_D 100
#define KEY_C 99


Mesh * mesh;
bool check = false;
bool check2 = false;
double rX = 0, rY = 0;

void idle() { // AGREGAR ESTA FUNCION
	glutPostRedisplay();
}

//funcion llamada a cada imagen
void glPaint(void) {

	//El fondo de la escena al color initial
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //CAMBIO
	glLoadIdentity();
	glRotatef(rX, 1.0, 0.0, 0.0);
	glRotatef(rY, 0.0, 1.0, 0.0);

	glOrtho(-15.0f, 15.0f, -15.0f, 15.0f, -15.0f, 15.0f);


	//Vertex* v = mesh->vertices[5];

	//mesh->Neighbors(3);
	//mesh->Paint(5);

	mesh->Draw(1); // ARISTAS
	mesh->Draw(0); // CARAS

	if (check) {
		mesh->loop();
		check = false;
	}

	if (check2) {
		mesh->colapse();
		check2 = false;
	}


	glFlush();
	glutSwapBuffers();

}

//
//inicializacion de OpenGL
//
void init_GL(void) {
	//Color del fondo de la escena
	glClearColor(0.79f, 0.93f, 0.76f, 0.0f); //(R, G, B, transparencia) en este caso un fondo negro

	//modo projeccion
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
}

//en el caso que la ventana cambie de tamaño
GLvoid window_redraw(GLsizei width, GLsizei height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();


}

GLvoid window_key(unsigned char key, int x, int y) {
	switch (key) {
	case KEY_D:
		check = true;
		break;

	case KEY_ESC:
		exit(0);
		break;

	case KEY_C:
		check2 = true;
		break;

	default:
		break;
	}
}
void specialKeys(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_RIGHT:
		rY += 10; break;
	case GLUT_KEY_LEFT:
		rY -= 10; break;
	case GLUT_KEY_UP:
		rX += 10; break;
	case GLUT_KEY_DOWN:
		rX -= 10; break;
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800, 800); //tamaño de la ventana
	glutInitWindowPosition(100, 100); //posicion de la ventana
	glutCreateWindow("Puntos"); //titulo de la ventana
	glEnable(GL_DEPTH_TEST);
	init_GL(); //funcion de inicializacion de OpenGL

	glutDisplayFunc(glPaint);
	glutReshapeFunc(&window_redraw);
	// Callback del teclado
	glutKeyboardFunc(&window_key);
	//glutMouseFunc(&OnMouseClick);d
	//glutMotionFunc(&OnMouseMotion);
	glutSpecialFunc(&specialKeys);
	glutIdleFunc(&idle);

	mesh = new Mesh;


	//mesh->leer("objs/icosaedro.obj",0);
	//mesh->leer("objs/cubo.obj",0);
	//mesh->leer("objs/conejo.obj",20);
	//mesh->leer("objs/piramide.obj",0);
	//mesh->leer("objs/cubin.obj",0);
	mesh->leer("objs/jirafa.obj",0);

	
	
	glutMainLoop(); //bucle de rendering
	return 0;
}