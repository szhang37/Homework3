#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rayTriangle.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////
void RayTriangle::initialize(void){
}
double RayTriangle::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
	return -1;
}
BoundingBox3D RayTriangle::setBoundingBox(void){
	Point3D pList[3];
	pList[0]=v[0]->position;
	pList[1]=v[1]->position;
	pList[2]=v[2]->position;
	bBox=BoundingBox3D(pList,3);
	for(int i=0;i<3;i++){
		bBox.p[0][i]-=RAYEPS;
		bBox.p[1][i]+=RAYEPS;
	}
	return bBox;
}

//////////////////
// OpenGL stuff //
//////////////////
int RayTriangle::drawOpenGL(int materialIndex){
	
	Point3D v1 = v[0]->position;
	Point3D v2 = v[1]->position;
	Point3D v3 = v[2]->position;
	Point3D n1 = v[0]->normal;
	Point3D n2 = v[1]->normal;
	Point3D n3 = v[2]->normal;
	Point2D tex1 = v[0]->texCoordinate;
	Point2D tex2 = v[1]->texCoordinate;
	Point2D tex3 = v[2]->texCoordinate;

	if (material->index != materialIndex) {
		
	 	material->drawOpenGL();
	 }

	glBegin(GL_TRIANGLES);
    printf("%f, %f, %f, %f, %f, %f\n", tex1[0],tex1[1],tex2[0],tex2[1],tex3[0],tex3[0]);
 	glTexCoord2f(tex1.p[0], tex1.p[1]);       
	glNormal3d(n1[0], n1[1], n1[2]);
	glVertex3d(v1[0], v1[1], v1[2]);
	

	glTexCoord2f(tex2.p[0], tex2.p[1]);
	glNormal3d(n2[0], n2[1], n2[2]);
	glVertex3d(v2[0], v2[1], v2[2]);

	glTexCoord2f(tex3.p[0], tex3.p[1]);
	glNormal3d(n3[0], n3[1], n3[2]);
	glVertex3d(v3[0], v3[1], v3[2]);

	glEnd();

	return materialIndex;
}
