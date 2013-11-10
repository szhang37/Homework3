#include <stdlib.h>
#include <GL/glut.h>
#include "rayGroup.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////
double RayGroup::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
	return -1;
}

BoundingBox3D RayGroup::setBoundingBox(void){
	Point3D* pList;
	BoundingBox3D tBBox;
	pList=new Point3D[sNum*2];
	for(int i=0;i<sNum;i++){
		tBBox=shapes[i]->setBoundingBox();
		pList[2*i  ]=tBBox.p[0];
		pList[2*i+1]=tBBox.p[1];
	}
	tBBox=BoundingBox3D(pList,sNum*2);

	delete[] pList;
	bBox=tBBox.transform(getMatrix());
	return bBox;
}

int StaticRayGroup::set(void){
	return 1;
}
//////////////////
// OpenGL stuff //
//////////////////
int RayGroup::getOpenGLCallList(void){
        //generate the new call list
        GLuint newList= glGenLists (1);
	
        //get new list
	glNewList(newList, GL_COMPILE);

        // Material Index is set to -1
        drawOpenGL(-1);

        glEndList();

        return newList;

}

int RayGroup::drawOpenGL(int materialIndex){

	 if (openGLCallListID == 0) {
                int i = 0;
                for (i = 0; i < sNum; i++) {

                        shapes[i]->drawOpenGL(materialIndex);
                }

        } else {
                int i = 0;
                for (i = 0; i < sNum; i++) {
                        glCallList(openGLCallListID);
                }

        }

        return -1;
}

/////////////////////
// Animation Stuff //
/////////////////////
Matrix4D ParametrizedRayGroup::getInverseMatrix(void){
	return Matrix4D::IdentityMatrix();
}
Matrix4D ParametrizedRayGroup::getNormalMatrix(void){
	return Matrix4D::IdentityMatrix();
}

Matrix4D ParametrizedEulerAnglesAndTranslation::getMatrix(void){
	return Matrix4D::IdentityMatrix();
}
Matrix4D ParametrizedClosestRotationAndTranslation::getMatrix(void){
	return Matrix4D::IdentityMatrix();
}
Matrix4D ParametrizedRotationLogarithmAndTranslation::getMatrix(void){
	return Matrix4D::IdentityMatrix();
}
Matrix4D ParametrizedQuaternionAndTranslation::getMatrix(void){
	return Matrix4D::IdentityMatrix();
}
