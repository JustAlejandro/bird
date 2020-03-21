#include "RigidBodyTemplate.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

namespace bird1 {

RigidBodyTemplate::RigidBodyTemplate(Eigen::Ref<Eigen::MatrixX3d> V,
                                     Eigen::Ref<Eigen::MatrixX3i> F,
                                     double scale) : volume_(0), radius_(0)
{
    inertiaTensor_.setZero();
    this->V = V * scale;
    this->F = F;
    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate()
{    
}

void RigidBodyTemplate::initialize()
{
    volume_ = computeVolume();
    com_ = computeCenterOfMass();
    translateCOM();
    computeRadius();
    // TODO: Translate center of mass to origin
    inertiaTensor_ = computeInertiaTensor();    
}

void RigidBodyTemplate::translateCOM(){
    for (int i = 0; i < V.rows(); i++)
    {
        V.row(i) = V.row(i) - com_.transpose();
    }
    com_ = Vector3d(0,0,0);
}

double RigidBodyTemplate::computeVolume()
{
    double volume = 0;
    int faceCount = F.rows();
    for (int i = 0; i < faceCount; i++)
    {
        Eigen::Vector3i indices = F.row(i); 
        Vector3d T0 = V.row(indices[0]);
        Vector3d T1 = V.row(indices[1]);
        Vector3d T2 = V.row(indices[2]);
        volume += 1.0/6.0 * (T0[0] + T1[0] + T2[0]) * (T1 - T0).cross(T2 - T0)[0];
    }
    return volume;
}

Vector3d RigidBodyTemplate::computeCenterOfMass()
{
    Vector3d cm(0, 0, 0);
    int faceCount = F.rows();

    for (int i = 0; i < faceCount; i++)
    {
        Eigen::Vector3i indices = F.row(i); 
        Vector3d T0 = V.row(indices[0]);
        Vector3d T1 = V.row(indices[1]);
        Vector3d T2 = V.row(indices[2]);

        Vector3d cur(0,0,0);
        cur += 1.0/2.0 * T0.cwiseProduct(T0);
        cur += 1.0/3.0 * (T1 - T0).cwiseProduct(T0);
        cur += 1.0/3.0 * (T2 - T0).cwiseProduct(T0);
        cur += 1.0/12.0 * (T1 - T0).cwiseProduct(T2 - T0);
        cur += 1.0/12.0 * (T1 - T0).cwiseProduct(T1 - T0);
        cur += 1.0/12.0 * (T2 - T0).cwiseProduct(T2 - T0);
        cur = cur.cwiseProduct((T1 - T0).cross(T2 - T0));
        cm += cur;
    }

    return (cm / 2.0) / volume_;
}

void RigidBodyTemplate::computeRadius(){
    radius_ = 0.0000001;
    for (int i = 0; i < V.rows(); i++)
    {
        radius_ = max(radius_, V.row(i).norm());
    }
}

// a^2 + b (2 a u + b u^2 + 2 c u v) + c (2 a v + c v^2)
double RigidBodyTemplate::computeInertiaQuadraticTerm(const double& a, const double& b, const double& c){
    return  (1.0/3.0) * ((1.0/2.0)*a*a*a + (1.0/2.0)*a*a*b + (1.0/2.0)*a*a*c 
        + (1.0/4.0)*a*b*b + (1.0/4.0)*a*b*c + (1.0/4.0)*a*c*c
        + (1.0/20.0)*b*b*b + (1.0/20.0)*b*b*c + (1.0/20.0)*b*c*c + (1.0/20.0)*c*c*c);
}



double computeXSquaredYSquaredTerm(Vector3d a, Vector3d b, Vector3d c, int type){
    int i0,i1,i2;
    switch(type){
        //0,0
        case 0:
            i0 = 0;
            i1 = 1;
            i2 = 2; 
            break;
        //1,1
        case 1:
            i0 = 1;
            i1 = 0;
            i2 = 2;
            break;
        case 2: 
            i0 = 2;
            i1 = 0;
            i2 = 1;
            break;
        default:
            throw "Invalid Call";
            break;
    }

    double x0 = a(i0), x1 = b(i0), x2 = c(i0);
    double y0 = a(i1), y1 = b(i1), y2 = c(i1);
    double z0 = a(i2), z1 = b(i2), z2 = c(i2);
    double toRet = (x0*y0*y0)/20 + (x0*y1*y1)/60 + (x1*y0*y0)/60 + (x0*y2*y2)/60 + (x1*y1*y1)/20 + (x2*y0*y0)/60 + (x1*y2*y2)/60 + (x2*y1*y1)/60 + (x2*y2*y2)/20 + (x0*z0*z0)/20 + (x0*z1*z1)/60 + (x1*z0*z0)/60 + (x0*z2*z2)/60 + (x1*z1*z1)/20 + (x2*z0*z0)/60 + (x1*z2*z2)/60 + (x2*z1*z1)/60 + (x2*z2*z2)/20 + (x0*y0*y1)/30 + (x0*y0*y2)/30 + (x1*y0*y1)/30 + (x0*y1*y2)/60 + (x1*y0*y2)/60 + (x2*y0*y1)/60 + (x1*y1*y2)/30 + (x2*y0*y2)/30 + (x2*y1*y2)/30 + (x0*z0*z1)/30 + (x0*z0*z2)/30 + (x1*z0*z1)/30 + (x0*z1*z2)/60 + (x1*z0*z2)/60 + (x2*z0*z1)/60 + (x1*z1*z2)/30 + (x2*z0*z2)/30 + (x2*z1*z2)/30;

    return toRet * (b-a).cross(c-a)(type);
}

double computeXYZTerm(Vector3d a, Vector3d b, Vector3d c, int type){
    int i0,i1,i2;
    i0 = 0;
    i1 = 1;
    i2 = 2; 

    double x0 = a(i0), x1 = b(i0), x2 = c(i0);
    double y0 = a(i1), y1 = b(i1), y2 = c(i1);
    double z0 = a(i2), z1 = b(i2), z2 = c(i2);

    double toRet = (x0*y0*z0)/20 + (x0*y0*z1)/60 + (x0*y1*z0)/60 + (x1*y0*z0)/60 + (x0*y0*z2)/60 + (x0*y1*z1)/60 + (x0*y2*z0)/60 + (x1*y0*z1)/60 + (x1*y1*z0)/60 + (x2*y0*z0)/60 + (x0*y1*z2)/120 + (x0*y2*z1)/120 + (x1*y0*z2)/120 + (x1*y1*z1)/20 + (x1*y2*z0)/120 + (x2*y0*z1)/120 + (x2*y1*z0)/120 + (x0*y2*z2)/60 + (x1*y1*z2)/60 + (x1*y2*z1)/60 + (x2*y0*z2)/60 + (x2*y1*z1)/60 + (x2*y2*z0)/60 + (x1*y2*z2)/60 + (x2*y1*z2)/60 + (x2*y2*z1)/60 + (x2*y2*z2)/20;
    return -1.0 * toRet * (b-a).cross(c-a)(type);
}

Eigen::Matrix3d
RigidBodyTemplate::computeInertiaTensor()
{
    Eigen::Matrix3d inertiaTensor;
    inertiaTensor.setZero();
    int faceCount = F.rows();
    // TODO: Computer Inertia Tensor
    for (int i = 0; i < faceCount; i++)
    {
        Eigen::Vector3i indices = F.row(i); 
        Vector3d T0 = V.row(indices[0]);
        Vector3d T1 = V.row(indices[1]);
        Vector3d T2 = V.row(indices[2]);

        //Compute Diagonal
        for(int i = 0; i < 3; i++)
            inertiaTensor(i,i) += computeXSquaredYSquaredTerm(T0, T1, T2, i);
        double integral = computeXYZTerm(T0, T1, T2, 2);
        inertiaTensor(0,1) += integral;
        inertiaTensor(1,0) += integral;

        integral = computeXYZTerm(T0, T1, T2, 1);
        inertiaTensor(0,2) += integral;
        inertiaTensor(2,0) += integral;

        integral = computeXYZTerm(T0, T1, T2, 0);
        inertiaTensor(1,2) += integral;
        inertiaTensor(2,1) += integral;

    }
    
    cout << inertiaTensor <<endl;
    return inertiaTensor;
}

}
