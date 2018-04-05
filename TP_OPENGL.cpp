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
#include <vector>
#include <string>
//#include "GL/freeglut.h"
#include "glad/glad.h"
#include <fstream>

/* Dans les salles de TP, vous avez généralement accès aux glut dans C:\Dev. Si ce n'est pas le cas, téléchargez les .h .lib ...
Vous pouvez ensuite y faire référence en spécifiant le chemin dans visual. Vous utiliserez alors #include <glut.h>. 
Si vous mettez glut dans le répertoire courant, on aura alors #include "glut.h" 
*/


#include "CVector.h"
#include "CPoint.h"
#include "Voxel.h"

using namespace std;

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
#define PRECISION 0.1
double nbM = 10;
double nbMc = 10;
double nbN = 10;
double ortho1 = -51.0;
double ortho2 = 51.0;
double ortho3 = -51.0;
double ortho4 = 51.0;
double ortho5 = -51.0;
double ortho6 = 51.0;

// Entêtes de fonctions
void init_scene();
void render_scene();
GLvoid initGL();
GLvoid window_display();
GLvoid window_reshape(GLsizei width, GLsizei height); 
GLvoid window_key(unsigned char key, int x, int y); 

struct mesh
{
    float * coordinates;
    int * indices;
    int nbPoints;
    int nbTriangles;
};

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
    gladLoadGL();
    glClearColor(RED, GREEN, BLUE, ALPHA);
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
    glOrtho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);


    // toutes les transformations suivantes s´appliquent au modèle de vue
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
        nbM++;
        nbN++;
        nbMc;
        glutPostRedisplay();
        break;

    case '-':
        if (nbM > 4){
            nbM--;
            nbMc--;
        }
        if (nbN > 3) nbN--;
        glutPostRedisplay();
        break;

    case 'z':;
        ortho4++;
        ortho5++;
        ortho6++;
        glutPostRedisplay();
        break;

    case 's':
        ortho4--;
        ortho5--;
        ortho6--;
        glutPostRedisplay();
        break;

    default:
        printf ("La touche %d n´est pas active.\n", key);
        break;
    }
}
CPoint somme(vector<double> asum, vector<CPoint> point){
    double resultx = 0;
    double resulty = 0;
    double resultz = 0;

    CPoint final;

    for (int i = 0;i < asum.size(); i++){
        resultx+=asum[i]*point[i].getX();
        cout << resultx << endl;
        resulty+=asum[i]*point[i].getY();
        resultz+=asum[i]*point[i].getZ();
    }
    final.setX(resultx);
    final.setY(resulty);
    final.setZ(resultz);

    return final;
}

void DrawCurve(vector<CPoint> TabPoint){
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < TabPoint.size(); i++){
        glVertex3f(TabPoint[i].getX(),TabPoint[i].getY(),TabPoint[i].getZ());
    }
    glEnd();
}

void DrawPoly(vector<CPoint> TabPoint){
    glBegin(GL_POLYGON);
    for (int i = 0; i < TabPoint.size(); i++){
        glVertex3f(TabPoint[i].getX(),TabPoint[i].getY(),TabPoint[i].getZ());
    }
    glEnd();
}

unsigned fact(unsigned x){
    if (x==0) return 1;
    return x*fact(x-1);
}
vector<CPoint> HermiteCubicCurve(CPoint p0,CPoint p1, CVector v0, CVector v1, long nbU){

    vector<CPoint> Pointst;

    for (int i = 0; i < nbU; i++){
        CPoint tmp;

        double j = (double) i/nbU;

        float F1 = 2*pow(j,3)-3*pow(j,2)+1;
        float F2 = -2*pow(j,3)+3*pow(j,2);
        float F3 = pow(j,3)-2*pow(j,2)+j;
        float F4 = pow(j,3)-pow(j,2);

        tmp.setX(F1*p0.getX()+F2*p1.getX()+F3*v0.getX()+F4*v1.getX());
        tmp.setY(F1*p0.getY()+F2*p1.getY()+F3*v0.getY()+F4*v1.getY());
        tmp.setZ(F1*p0.getZ()+F2*p1.getZ()+F3*v0.getZ()+F4*v1.getZ());

        cout << tmp.getX() << tmp.getY() << tmp.getZ() << endl;

        Pointst.push_back(tmp);
    }

    return Pointst;

}
vector<CPoint> BezierCurveByBernstein(vector<CPoint> tabControl,long nbU){
    vector<CPoint> final;
    for (int i = 0; i < nbU; i++){
        double k = (double) i/nbU;
        vector<double> polynomes;
        for (int j = 0; j < tabControl.size(); j++){
            double f1 = (double) fact(tabControl.size()-1)/(fact(j)*fact(tabControl.size()-1 -j));
            double f2 = (double) pow(1-k,tabControl.size()-1 -j);

            double tmp = f1*pow(k,j)*f2;
            cout << tmp << endl;
            polynomes.push_back(tmp);
        }

        CPoint x;
        x = somme(polynomes,tabControl);
        final.push_back(x);
    }

    return final;
}

vector<CPoint> BezierCurveByCasteljau(vector<CPoint> points,long nbu)
{
    vector<CPoint> bezierPoints;

    for(double t = 0 ; t < nbu ; t++)
    {
        double k = (double) t/nbu;

        vector<CPoint> tmp1 = points;

        while(tmp1.size()>1)
        {
            vector<CPoint> tmp2;

            for(int i = 0 ; i<tmp1.size()-1 ; i++)
            {
                CPoint Ptmp1 = tmp1[i];
                CPoint Ptmp2 = tmp1[i+1];

                CPoint tmp = CPoint(Ptmp1.getX() * (1-k) + Ptmp2.getX() * k,Ptmp1.getY() * (1-k) + Ptmp2.getY() * k,0);
                tmp.drawPoint();
                tmp2.push_back(tmp);
            }

            tmp1 = tmp2;
        }

        bezierPoints.push_back(tmp1[0]);
    }

    return bezierPoints;
}
vector<vector<CPoint> > traceSurfaceCylindrique(vector<CPoint> points, CVector v1, long nbv,long nbu){
    vector<vector<CPoint> > result;
    vector<CPoint> cptmp;

    vector<CPoint> bezier = BezierCurveByBernstein(points,nbu);

    for (int j = 0; j < bezier.size(); j++){
        CPoint tmp1 = CPoint(bezier[j].getX()+v1.getX(),bezier[j].getY()+v1.getY(),0);
        cptmp.push_back(bezier[j]);
        cptmp.push_back(tmp1);
        result.push_back(cptmp);
        cptmp.clear();
    }

    vector<CPoint> copy = points;

    double range = (double) 1/nbv;
    for (int i = 0; i < nbv; i++){
        bezier = BezierCurveByBernstein(copy,nbu);
        for (int l = 0; l < copy.size(); l++){
            copy[l].setX(copy[l].getX()+range);
            copy[l].setY(copy[l].getX()+range);
        }

        result.push_back(bezier);
        bezier.clear();
    }

    return result;
}

/* 
void DisplayVoxel(double length,CPoint Centre){
    vector <CPoint> VoxelP1;

    VoxelP1.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() - length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() - length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));

    DrawPoly(VoxelP1);

    vector <CPoint> VoxelP2;

    VoxelP2.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() - length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() - length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));

    DrawPoly(VoxelP2);

    for (unsigned i (0); i < VoxelP1.size() - 1; ++i) {
        glBegin(GL_POLYGON);
        glVertex3f(VoxelP1[i].getX(), VoxelP1[i].getY(), VoxelP1[i].getZ());
        glVertex3f(VoxelP2[i].getX(), VoxelP2[i].getY(), VoxelP2[i].getZ());
        glVertex3f(VoxelP2[i+1].getX(), VoxelP2[i+1].getY(), VoxelP2[i+1].getZ());
        glVertex3f(VoxelP1[i+1].getX(), VoxelP1[i+1].getY(), VoxelP1[i+1].getZ());
        glEnd();
    }

}

/* */

void displayVoxelRecSphere(Voxel* v, CPoint* ct, double ray,double reso, double i) {
  if(i<reso)
  {
    for(int k=0; k<8; k++)
    {
      Voxel* vtemp = new Voxel(&v->getSubCenterPoints()->at(k),(double)v->getLength()/2.0);
      int test = vtemp->isInsideSphere(ct, ray);
      //std::cout<<"test ="<<test<<std::endl;
      if(test!=2)
      {
        cout << "" << endl;
      }
      else // test == 2
      {
        if(i == reso-1)
          vtemp->displayV();
        else 
          displayVoxelRecSphere(vtemp, ct,  ray, reso, i+1);
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////
// Fonction que vous allez modifier afin de dessiner
/////////////////////////////////////////////////////////////////////////////////////////


double rand_float(double a, double b) {
    return ((double)rand() / RAND_MAX) * (b - a) + a;
}

vector<vector<CPoint>> Cylindre(int nbMeridiens){
    double h = 30.0;
    double r = 10.0;

    vector<CPoint> tmp1;
    vector<CPoint> tmp2;
    vector<vector<CPoint> > result;

    for(double i = 0.0; i < nbMeridiens; i++){
        double omega = 2.0 * (atan(1) * 4) *  i / nbMeridiens;

        double x = r * sin(omega);
        double y = r * cos(omega);
        double z = 0.0;
        double z2 = h;

        CPoint p = CPoint(x,y,z);
        CPoint p2 = CPoint(x,y,z2);

        tmp1.push_back(p);
        tmp2.push_back(p2);


    }
    result.push_back(tmp1);

    result.push_back(tmp2);


    double* Points;
    unsigned cptPoints = 0;
    unsigned* Indices;
    unsigned cptIndices = 0;

    for (unsigned i (0); i < result[0].size(); i++) {
        Points[cptPoints++] = result[0][i].getX();
        Points[cptPoints++] = result[0][i].getY();
        Points[cptPoints++] = result[0][i].getZ();

        Points[cptPoints++] = result[1][i].getX();
        Points[cptPoints++] = result[1][i].getY();
        Points[cptPoints++] = result[1][i].getZ();

        if (i >= 1) {
            Indices[cptIndices++] = cptPoints - 12;
            Indices[cptIndices++] = cptPoints - 9;
            Indices[cptIndices++] = cptPoints - 6;

            Indices[cptIndices++] = cptPoints - 3;
            Indices[cptIndices++] = cptPoints - 6;
            Indices[cptIndices++] = cptPoints - 9;
        }
    }

    return result;

}

vector<vector<CPoint>> Sphere(int rayon, double nbMeridiens, double nbParallelles){
    vector<vector<CPoint> > result;

    CPoint nord = CPoint(0,rayon,0);
    CPoint sud = CPoint(0,-rayon,0);
    glBegin(GL_POINTS);
    glColor3f(rand_float(0.0,1.0), rand_float(0.0,1.0), rand_float(0.0,1.0));
    glColor3f(rand_float(0.0,1.0), rand_float(0.0,1.0), rand_float(0.0,1.0));
    sud.drawPoint();
    glColor3f(1.0,1.0,1.0);
    glEnd();


    for (int i = 0; i < nbParallelles; i++){
        double phi = (atan(1) * 4) * i / nbParallelles;
        vector<CPoint> tmp1;
        for (int j = 0; j < nbMeridiens; j++){
            double the = 2 * (atan(1) * 4) * j / nbMeridiens;
            double x = rayon * sin(phi) * cos(the);
            double z = rayon * sin(phi) * sin(the);
            double y = rayon * cos(phi);
            CPoint t = CPoint(x,y,z);
            glBegin(GL_POINTS);
            t.drawPoint();
            glEnd();
            tmp1.push_back(t);
        }
        result.push_back(tmp1);
        tmp1.clear();
    }

    return result;


}

void render_scene()
{


    // Exercice 1

    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // glEnable(GL_LIGHTING);
    // glEnable(GL_LIGHT0);

    glColor3f(1.0, 1.0, 1.0);
    //glMaterial(100);
    glOrtho(ortho1,ortho2,ortho3,ortho4,ortho5,ortho6);
    gluLookAt(15.0,20.0,15.0,0,10,0,0,1,0);
    cout << "Avant remplissage de la sphere ! " << endl;
    vector<vector<CPoint>> result = Sphere(20,nbM,nbN);
    cout << "Apres le remplissage de la sphere ! " << endl;
    //Cylindre(nbMc);
    //Sphere(20,nbM,nbN);

    double* Points;
    unsigned cptPoints = 0;
    unsigned* Indices;
    unsigned cptIndices = 0;
    
    for (unsigned i (0); i < result.size(); ++i) {
        for (unsigned j (0); j < result[i].size(); ++j) {
            if (i == 0) {
                Points[cptPoints++] = result[i][j].getX();
                Points[cptPoints++] = result[i][j].getY();
                Points[cptPoints++] = result[i][j].getZ();
            }
            else {
                 if (j == 0) {
                    Points[cptPoints++] = result[i][j].getX();
                    Points[cptPoints++] = result[i][j].getY();
                    Points[cptPoints++] = result[i][j].getZ();
                    
                    Points[cptPoints++] = result[i][result[i].size() - 1].getX();
                    Points[cptPoints++] = result[i][result[i].size() - 1].getY();
                    Points[cptPoints++] = result[i][result[i].size() - 1].getZ();
                    
                    Points[cptPoints++] = result[i-1][j].getX();
                    Points[cptPoints++] = result[i-1][j].getY();
                    Points[cptPoints++] = result[i-1][j].getZ();
                    
                    Points[cptPoints++] = result[i-1][result[i].size() - 1].getX();
                    Points[cptPoints++] = result[i-1][result[i].size() - 1].getY();
                    Points[cptPoints++] = result[i-1][result[i].size() - 1].getZ(); 
                }

                 else {
                    Points[cptPoints++] = result[i][j].getX();
                    Points[cptPoints++] = result[i][j].getY();
                    Points[cptPoints++] = result[i][j].getZ();
                    
                    Points[cptPoints++] = result[i][j-1].getX();
                    Points[cptPoints++] = result[i][j-1].getY();
                    Points[cptPoints++] = result[i][j-1].getZ();
                    
                    Points[cptPoints++] = result[i-1][j].getX();
                    Points[cptPoints++] = result[i-1][j].getY();
                    Points[cptPoints++] = result[i-1][j].getZ();
                    
                    Points[cptPoints++] = result[i-1][j-1].getX();
                    Points[cptPoints++] = result[i-1][j-1].getY();
                    Points[cptPoints++] = result[i-1][j-1].getZ();     
                }
                Indices[cptIndices++] = cptPoints - 12;
                Indices[cptIndices++] = cptPoints - 9;
                Indices[cptIndices++] = cptPoints - 6;
                
                Indices[cptIndices++] = cptPoints - 3;
                Indices[cptIndices++] = cptPoints - 6;
                Indices[cptIndices++] = cptPoints - 9;
            }
        }
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    //glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, Points);
    //glDrawElements(GL_TRIANGLES, 3 * NbTriangles, GL_UNSIGNED_INT, Indices);
    for (unsigned i (0) ; i < cptIndices; i++) {
        glColor3f(rand_float(0.0,1.0), rand_float(0.0,1.0), rand_float(0.0,1.0));
        glDrawElements(GL_TRIANGLES, 3 * i, GL_UNSIGNED_INT, Indices);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    /* */


}
