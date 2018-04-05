#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "CPoint.h"
#include "CVector.h"

using namespace std;

CPoint::CPoint(){}
CPoint::CPoint(double x2,double y2,double z2) : x (x2), y (y2), z (z2) {}
CPoint::CPoint(const CPoint &p){
    x = p.x;
    y = p.y;
    z = p.z;
}
double CPoint::getX(){
    return x;
}
double CPoint::getY(){
    return y;
}
double CPoint::getZ(){
    return z;
}
void CPoint::setX(double x2){
    x=x2;
}
void CPoint::setY(double y2){
    y=y2;
}
void CPoint::setZ(double z2){
    z=z2;
}

int CPoint::isInsideSphere(CPoint* c,double r){
    cout << "isinside spehere point" << endl;
    double a = c->getX() - this->x;
    cout << "a" << endl;
    double b = c->getY() - this->y;
    cout << "b" << endl;
    double c2 = c->getZ() - this->z;
    double dist =  sqrt(pow(a,2)+pow(b,2)+pow(c2,2));
    cout << "fin calcul distance point " << endl;
    if (dist <= r)
        return 1;
    else return 0;

}

int CPoint::isInsideCylinder(CPoint* c,CVector v,double r){
    CPoint proj = ProjectOnLine(v,c);

    double a = proj.getX() - this->x;
    cout << "a" << endl;
    double b = proj.getY() - this->y;
    cout << "b" << endl;
    double c2 = proj.getZ() - this->z;
    double dist =  sqrt(pow(a,2)+pow(b,2)+pow(c2,2));
    cout << "fin calcul distance point " << endl;
    if (dist <= r)
        return 1;
    else return 0;

}

CPoint CPoint::ProjectOnLine(CPoint pl1,CPoint pl2){
    CVector v = CVector(pl1.getX(), pl1.getY(), pl1.getZ(), pl2.getX(), pl2.getY(), pl2.getZ());
    CPoint* p = new CPoint(pl1.getX(),pl1.getY(),pl1.getZ());
    return ProjectOnLine(v,p);
}

CPoint CPoint::ProjectOnLine(CVector v,CPoint* pl){

    CVector v3= CVector(x-pl->getX(),y-pl->getY(),z-pl->getZ());//vecteur point a projeter
    
    double n = v.Norme();
    double dist = v.Scalar(v3)/n;
    v.Normalize();

    double xtmp = pl->getX()+v.getX()*dist;
    double ytmp = pl->getY()+v.getY()*dist;
    double ztmp = pl->getZ()+v.getZ()*dist;

    CPoint r = CPoint(xtmp,ytmp,ztmp);
    return r;
}

CPoint CPoint::ProjectOnPlan(CPoint pplan,CVector normalplan){
    CVector v2 = CVector(x, y, z, pplan.getX(), pplan.getY(), pplan.getZ());
    double dist = v2.Scalar(normalplan)/normalplan.Norme();

    double xtmp = x-normalplan.getX()*dist;
    double ytmp = y-normalplan.getY()*dist;
    double ztmp = z-normalplan.getZ()*dist;

    CPoint r = CPoint(xtmp,ytmp,ztmp);
    return r;
}

void CPoint::drawPoint(){
    //glBegin(GL_POINTS);
    glVertex3f(x,y,z);
    //glEnd();
}

