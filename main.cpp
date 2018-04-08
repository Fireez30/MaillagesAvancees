///////////////////////////////////////////////////////////////////////////////
// Imagina
// ----------------------------------------------------------------------------
// IN - Synthèse d'images - Modélisation géométrique
// Auteur : Gilles Gesquière
// ----------------------------------------------------------------------------
// Base du TP 1
// programme permettant de créer des formes de bases.
// La forme représentée ici est un polygone blanc dessiné sur un fond rouge
///////////////////////////////////////////////////////////////////////////////  

#include <stdio.h>      
#include <stdlib.h>     
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;
/* Dans les salles de TP, vous avez généralement accès aux glut dans C:\Dev. Si ce n'est pas le cas, téléchargez les .h .lib ...
Vous pouvez ensuite y faire référence en spécifiant le chemin dans visual. Vous utiliserez alors #include <glut.h>.
Si vous mettez glut dans le répertoire courant, on aura alors #include "glut.h"
*/

#include <GL/glut.h> 
#include "Point.h"
#include "Vector.h"

// Définition de la taille de la fenêtre
#define WIDTH  480

#define HEIGHT 480

// Définition de la couleur de la fenêtre
#define RED   0
#define GREEN 0
#define BLUE  0
#define ALPHA 1


// Touche echap (Esc) permet de sortir du programme
#define KEY_ESC 27

// PI
# define M_PI           3.14159265358979323846  /* pi */

// Entêtes de fonctions
void init_scene();
void render_scene();
GLvoid initGL();
GLvoid window_display();
GLvoid window_reshape(GLsizei width, GLsizei height);
GLvoid window_key(unsigned char key, int x, int y);

double last_meridien = 3, last_parallele = 2;
double radian = 0, meridient = 3, parallele = 2;
bool dessine = true;

int main(int argc, char **argv)
{
	// initialisation  des paramètres de GLUT en fonction
	// des arguments sur la ligne de commande
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA);

	// définition et création de la fenêtre graphique, ainsi que son titre
	glutInitWindowSize(WIDTH, HEIGHT);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Premier exemple : carré");

	// initialisation de OpenGL et de la scène
	initGL();
	init_scene();

	// choix des procédures de callback pour 
	// le tracé graphique
	glutDisplayFunc(&window_display);
	// le redimensionnement de la fenêtre
	glutReshapeFunc(&window_reshape);
	// la gestion des événements clavier
	glutKeyboardFunc(&window_key);

	// la boucle prinicipale de gestion des événements utilisateur
	glutMainLoop();

	return 1;
}

// initialisation du fond de la fenêtre graphique : noir opaque
GLvoid initGL()
{
	glClearColor(RED, GREEN, BLUE, ALPHA);
	glEnable(GL_DEPTH_TEST); 	// Active le test de profondeur
	glEnable(GL_LIGHTING); 	// Active l'éclairage
	glEnable(GL_LIGHT0);

}

// Initialisation de la scene. Peut servir à stocker des variables de votre programme
// à initialiser
void init_scene()
{
}

// fonction de call-back pour l´affichage dans la fenêtre

GLvoid window_display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	// C'est l'endroit où l'on peut dessiner. On peut aussi faire appel
	// à une fonction (render_scene() ici) qui contient les informations 
	// que l'on veut dessiner
	render_scene();

	// trace la scène grapnique qui vient juste d'être définie
	glFlush();
}

// fonction de call-back pour le redimensionnement de la fenêtre


GLvoid window_reshape(GLsizei width, GLsizei height)
{
	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// ici, vous verrez pendant le cours sur les projections qu'en modifiant les valeurs, il est
	// possible de changer la taille de l'objet dans la fenêtre. Augmentez ces valeurs si l'objet est 
	// de trop grosse taille par rapport à la fenêtre.
	glOrtho(-40.0, 40.0, -40.0, 40.0, -40.0, 40.0);
	// toutes les transformations suivantes s´appliquent au modèle de vue 
	glMatrixMode(GL_MODELVIEW);
}

GLvoid window_reshape(GLsizei width, GLsizei height, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max)
{
	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	float x_min_f = x_min - (abs(x_max - x_min));
	float x_max_f = x_max + (abs(x_max - x_min));
	float y_min_f = y_min - (abs(y_max - y_min));
	float y_max_f = y_max + (abs(y_max - y_min));
	float z_min_f = z_min - (abs(z_max - z_min));
	float z_max_f = z_max + (abs(z_max - z_min));

	glOrtho(x_min_f, x_max_f, y_min_f, y_max_f, z_min_f, z_max_f);
	glMatrixMode(GL_MODELVIEW);
}

// fonction de call-back pour la gestion des événements clavier

GLvoid window_key(unsigned char key, int x, int y)
{
	switch (key) {
	case KEY_ESC:
		exit(1);
		break;
	case '+':
		parallele++;
		meridient++;
		dessine = true;
		break;
	case '-':
		if (parallele > 2)
			parallele--;
		if (meridient > 3)
			meridient--;
		dessine = true;
		break;
	default:
		printf("La touche %d n´est pas active.\n", key);
		dessine = false;
		break;
	}
	glutPostRedisplay();
}

int factorielle(int val)
{
	int result = 1;
	for (int i = val; i > 0; --i)
		result *= i;
	return result;
}

void normcrossprod(float v1[3], float v2[3], float out[3])
{
	GLint i, j;
	GLfloat length;

	out[0] = v1[1] * v2[2] - v1[2] * v2[1];
	out[1] = v1[2] * v2[0] - v1[0] * v2[2];
	out[2] = v1[0] * v2[1] - v1[1] * v2[0];
	//normalize(out); 
}

void renderMaillage(GLfloat coordinates[], GLfloat all_normal[], int indices[], int nbTriangles) {
	glEnable(GL_COLOR_MATERIAL);
	int MatSpec[4] = { 0,0,0,0 };
	glMaterialiv(GL_FRONT_AND_BACK, GL_SPECULAR, MatSpec);
	glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 128);

	// Créer un tableau de normal et utiliser "glNormalPointer". Se baser sur le système des vertex.
	// https://www.opengl.org/discussion_boards/showthread.php/130241-using-glNormalPointer

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, coordinates);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, all_normal);
	glDrawElements(GL_TRIANGLES, 3 * nbTriangles, GL_UNSIGNED_INT, indices);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

Point** calculeSphere(double rayon, int nbMeridient, int nbParallele) {
	Point** tab = new Point*[nbParallele];
	nbParallele++;
	for (int i = 0; i < nbParallele; i++) {
		tab[i] = new Point[nbMeridient];
	}
	const double PI = atan(1) * 4;
	double angleBase = 0, angleHauteur = PI / (double)(nbParallele);
	for (int i = 0; i < nbParallele; i++) {
		angleBase = 0;
		for (int j = 0; j < nbMeridient; j++) {
			double x = rayon * sin(angleHauteur) * cos(angleBase);
			double y = rayon * cos(angleHauteur);
			double z = rayon * sin(angleHauteur) * sin(angleBase);
			tab[i][j] = *new Point(x, y, z);
			glVertex3f(x, y, z);
			angleBase += 2.0 * PI / (double)(nbMeridient);
		}
		angleHauteur += PI / (double)nbParallele;
	}
	return tab;
}

void DrawSphere(int nbMeridient, int nbParallele, int rayon) {
	Point ** sphere = calculeSphere(rayon, nbMeridient, nbParallele);

	GLfloat* coordinates;
	GLfloat* all_normal;
	GLubyte* colors;
	int* indices, *adjacents;

	int nbPoints = nbMeridient * nbParallele + 2, nbTriangles = nbMeridient * nbParallele * 2;
	float x_min = -rayon, x_max = rayon, y_min = -rayon, y_max = rayon, z_min = -rayon, z_max = rayon;

	coordinates = new GLfloat[nbPoints * 3];
	all_normal = new GLfloat[nbTriangles * 3 ];
	colors = new GLubyte[nbPoints*3];

	for (int i = 0; i < nbPoints * 3; i++) {
		colors[i] = 255;
	}

	indices = new int[nbTriangles * 3];
	adjacents = new int[nbTriangles * 3];

	for (int i = 0; i < nbTriangles * 3; i++) {
		adjacents[i] = -1;
	}

	int k = 0;
	for (int i = 0; i < nbParallele; i++) {
		for (int j = 0; j < nbMeridient; j++) {
			coordinates[k * 3] = sphere[i][j].getX();
			coordinates[k * 3 + 1] = sphere[i][j].getY();
			coordinates[k * 3 + 2] = sphere[i][j].getZ();
			k++;
		}
	}
	coordinates[k * 3] = 0;
	coordinates[k * 3 + 1] = rayon;
	coordinates[k * 3 + 2] = 0;
	k++;
	coordinates[k * 3] = 0;
	coordinates[k * 3 + 1] = -rayon;
	coordinates[k * 3 + 2] = 0;

	Vector normRef;
	int nbTr = 0;
	for (int i = 0; i < nbParallele - 1; i++) {
		for (int j = 0; j < nbMeridient; j++) {
			indices[nbTr * 3] = i * nbMeridient + j;
			indices[nbTr * 3 + 1] = (i + 1) * nbMeridient + j;
			if (j + 1 < nbMeridient)
				indices[nbTr * 3 + 2] = (i + 1) * nbMeridient + j + 1;
			else
				indices[nbTr * 3 + 2] = (i + 1) * nbMeridient;
			nbTr++;

			indices[nbTr * 3] = i * nbMeridient + j;
			if (j + 1 < nbMeridient) {
				indices[nbTr * 3 + 1] = (i + 1) * nbMeridient + j + 1;
				indices[nbTr * 3 + 2] = i * nbMeridient + j + 1;
			}
			else {
				indices[nbTr * 3 + 1] = (i + 1) * nbMeridient;
				indices[nbTr * 3 + 2] = i * nbMeridient;
			}
			nbTr++;
		}
	}

	for (int i = 0; i < nbMeridient; i++) {
		indices[nbTr * 3] = nbPoints - 2;
		indices[nbTr * 3 + 1] = i;
		if (i + 1 <nbMeridient)
			indices[nbTr * 3 + 2] = i + 1;
		else
			indices[nbTr * 3 + 2] = 0;
		nbTr++;
	}

	for (int i = 0; i < nbMeridient; i++) {
		indices[nbTr * 3 + 1] = (nbParallele - 1) * nbMeridient + i;
		indices[nbTr * 3] = nbPoints - 1;
		if (i + 1 < nbMeridient)
			indices[nbTr * 3 + 2] = (nbParallele - 1) * nbMeridient + i + 1;
		else
			indices[nbTr * 3 + 2] = (nbParallele - 1) * nbMeridient;
		nbTr++;
	}

	for (int i = 0; i < nbTriangles; i++) {
		Point* p1 = new Point(coordinates[indices[3 * i] * 3], coordinates[indices[3 * i] * 3 + 1], coordinates[indices[3 * i] * 3 + 2]);
		Point* p2 = new Point(coordinates[indices[3 * i + 1] * 3], coordinates[indices[3 * i + 1] * 3 + 1], coordinates[indices[3 * i + 1] * 3 + 2]);
		Point* p3 = new Point(coordinates[indices[3 * i + +2] * 3], coordinates[indices[3 * i + 2] * 3 + 1], coordinates[indices[3 * i + 2] * 3 + 2]);
		Vector* v1 = new Vector(*p1, *p2);
		Vector* v2 = new Vector(*p1, *p3);

		v1->Normalize();
		v2->Normalize();

		GLfloat t1[3] = { v1->getX(), v1->getY(), v1->getZ() };
		GLfloat t2[3] = { v2->getX(), v2->getY(), v2->getZ() };
		GLfloat out[3];

		normcrossprod(t1, t2, out);
		if (i == 0) {	
			all_normal[i * 3] = out[0];
			all_normal[i * 3 + 1] = out[1];
			all_normal[i * 3 + 2] = out[2];
			normRef = Vector(out[0], out[1], out[2]);
		}
		else {
			Vector v = Vector(out[0], out[1], out[2]);
			if (v.Scalar(normRef) > 0) {
				all_normal[i * 3] = out[0];
				all_normal[i * 3 + 1] = out[1];
				all_normal[i * 3 + 2] = out[2];
				normRef = v;
			}
			else {
				all_normal[i * 3] = -out[0];
				all_normal[i * 3 + 1] = -out[1];
				all_normal[i * 3 + 2] = -out[2];
				normRef = Vector(-out[0], -out[1], -out[2]);
			}
		}

		all_normal[i * 3] = out[0];
		all_normal[i * 3 + 1] = out[1];
		all_normal[i * 3 + 2] = out[2];

		Point* p20 = new Point(all_normal[i * 3], all_normal[i * 3 + 1], all_normal[i * 3 + 2]);

		glBegin(GL_POINTS);
		p20->DrawPoint();
		glEnd();
	}

	for (int i = 0; i < nbTriangles * 3; i += 3) {
		int nbSimilaire = 0;
		for (int j = 0; j < i; j += 3) {
			int similaire = 0;
			for (int k = 0; k < 3; k++) {
				if (indices[i] == indices[j + k]) {
					similaire++;
				}
				if (indices[i + 1] == indices[j + k]) {
					similaire++;
				}
				if (indices[i + 2] == indices[j + k]) {
					similaire++;
				}
			}
			if (similaire >= 2) {
				adjacents[i + nbSimilaire] = j / 3;
				int n = 0;
				while (n < 3 && adjacents[j + n] != -1) {
					n++;
				}
				if (n < 3) {
					adjacents[j + n] = i/3;
				}
				nbSimilaire++;
			}
		}
	}

	/*for (int i = 0; i < nbTriangles; i++) {
		cout << adjacents[i * 3] << " " << adjacents[i * 3 + 1] << " " << adjacents[i * 3 + 2] << endl;
	}*/

	/*for (int i = 0; i < nbTriangles * 3; i+=3) {
		cout << i/3 << " : " << indices[i] << " " << indices[i+1] << " " << indices[i+2] << endl;
	}*/

	/*for(int i = 0; i < nbPoints * 3; i += 3) {
		cout << i / 3 << " : " << coordinates[i] << " " << coordinates[i+1] << " " << 	coordinates[i+2] << endl;
	}*/


	glColor3f(1.0, 1.0, 1.0);
	renderMaillage(coordinates, all_normal, indices, nbTriangles);


	for (int i = 0; i < nbTriangles; i++) {
		Vector* v1 = new Vector(all_normal[i], all_normal[i + 1], all_normal[i + 2]);
		for (int j = 0; j < 3; j++) {
			int indice = adjacents[i*3 + j];
			Vector* v2 = new Vector(all_normal[indice], all_normal[indice + 1], all_normal[indice + 2]);
			double angle = v2->Angle(*v1);
			angle = (angle * 180) / M_PI;
			if (angle > 50.0 && i < indice) {
				cout << i << " " << indice << " " << angle << endl;
				Point* p1 = NULL;
				Point* p2 = NULL;
				for (int k = 0; k < 3; k++) {
					if (indices[i * 3] == indices[indice * 3 + k]) {
						if (p1 == NULL) {
							p1 = new Point(coordinates[indices[indice * 3 + k]*3], coordinates[indices[indice * 3 + k]*3+1], coordinates[indices[indice * 3 + k] * 3 + 2]);
						}
						else {
							p2 = new Point(coordinates[indices[indice * 3 + k] * 3], coordinates[indices[indice * 3 + k] * 3 + 1], coordinates[indices[indice * 3 + k] * 3 + 2]);
						}
					}
					if (indices[i * 3 + 1] == indices[indice * 3 + k]) {
						if (p1 == NULL) {
							p1 = new Point(coordinates[indices[i * 3 + 1]*3], coordinates[indices[i * 3 + 1] * 3 + 1], coordinates[indices[i * 3 + 1] * 3 + 2]);
						}
						else {
							p2 = new Point(coordinates[indices[i * 3 + 1] * 3], coordinates[indices[i * 3 + 1] * 3 + 1], coordinates[indices[i * 3 + 1] * 3 + 2]);
						}
					}
					if (indices[i * 3 + 2] == indices[indice * 3 + k]) {
						if (p1 == NULL) {
							p1 = new Point(coordinates[indices[i * 3 + 2] * 3], coordinates[indices[i * 3 + 2] * 3 + 1], coordinates[indices[i * 3 + 2] * 3 + 2]);
						}
						else {
							p2 = new Point(coordinates[indices[i * 3 + 2] * 3], coordinates[indices[i * 3 + 2] * 3 + 1], coordinates[indices[i * 3 + 2] * 3 + 2]);
						}
					}
				}
				Vector* v3 = new Vector(*p1, *p2);
				v3->DrawLine(*p1);
			}
		}
	}
}

void Cylindre_Tri(int rayon, int hauteur, int nbMeridien) {
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	int nbPoints = nbMeridien * 2, nbTriangles = nbMeridien * 2;
	GLfloat* coordinates = new GLfloat[nbPoints * 3];
	GLfloat* all_normal = new GLfloat[nbTriangles * 3];
	int* indices = new int[nbTriangles * 3];

	int k = 0;

	for (int i = 0; i < nbMeridien; i++) {
		double omega = 2 * M_PI * (double)i / (double)nbMeridien;
		double x = rayon * cos(omega);
		double y = rayon * sin(omega);
		double zH = (double)hauteur / 2.0;
		double zB = -(double)hauteur / 2.0;

		coordinates[3 * k] = x;
		coordinates[3 * k + 1] = zB;
		coordinates[3 * k + 2] = y;

		coordinates[3 * k + 3] = x;
		coordinates[3 * k + 4] = zH;
		coordinates[3 * k + 5] = y;

		k += 2;
	}

	for (int i = 0; i < nbTriangles; i++) {
		indices[3 * i] = i % nbTriangles;
		indices[3 * i + 1] = (i + 1) % nbTriangles;
		indices[3 * i + 2] = (i + 2) % nbTriangles;

		Point* p1 = new Point(coordinates[indices[3 * i] * 3], coordinates[indices[3 * i] * 3 + 1], coordinates[indices[3 * i] * 3 + 2]);
		Point* p2 = new Point(coordinates[indices[3 * i + 1] * 3], coordinates[indices[3 * i + 1] * 3 + 1], coordinates[indices[3 * i + 1] * 3 + 2]);
		Point* p3 = new Point(coordinates[indices[3 * i + 2] * 3], coordinates[indices[3 * i + 2] * 3 + 1], coordinates[indices[3 * i + 2] * 3 + 2]);
		Vector* v1 = new Vector(*p1, *p2);
		Vector* v2 = new Vector(*p1, *p3);

		v1->Normalize();
		v2->Normalize();

		GLfloat t1[3] = { v1->getX(), v1->getY(), v1->getZ() };
		GLfloat t2[3] = { v2->getX(), v2->getY(), v2->getZ() };
		GLfloat out[3];

		normcrossprod(t1, t2, out);

		if (i > 0) {
			Vector* v3 = new Vector(out[0], out[1], out[2]);
			Vector* v4 = new Vector(all_normal[(i - 1) * 3], all_normal[(i - 1) * 3 + 1], all_normal[(i - 1) * 3 + 2]);
			//Vector* v4 = new Vector(all_normal[(i-1) * 3 +1], all_normal[(i-1) * 3 + 2], all_normal[(i-1) * 3 + 3]);

			double k = v4->Scalar(*v3);
			if (k < 0) {
				out[0] *= -1;
				out[1] *= -1;
				out[2] *= -1;
			}
		}

		cout << out[0] << " " << out[1] << " " << out[2] << endl;

	   Point* p20 = new Point(out[0], out[1], out[2]);

	   glBegin(GL_POINTS);
	   p20->DrawPoint();
	   glEnd();
		
		all_normal[i * 3] = out[0];
		all_normal[i * 3 + 1] = out[1];
		all_normal[i * 3 + 2] = out[2];
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, coordinates);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, all_normal);
	glDrawElements(GL_TRIANGLES, 3 * nbTriangles, GL_UNSIGNED_INT, indices);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

//////////////////////////////////////////////////////////////////////////////////////////
// Fonction que vous allez modifier afin de dessiner
/////////////////////////////////////////////////////////////////////////////////////////
void render_scene()
{
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Définition de la couleur

	//readMaillage("bunny.off");
	int rayon = 10, hauteur = 20;
	gluLookAt(rayon + 5, hauteur/2, rayon + 5, 0, 0, 0, 0, 1, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//Cylindre_Tri(rayon, hauteur, meridient);
	DrawSphere(meridient, parallele, rayon);
}
