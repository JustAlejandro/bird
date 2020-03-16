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
    // //(z2 + y2)x
    // //First Term
    // double firstTerm = 5.0 * a(i0) * c(i1) * c(i1) + 10.0 * a(i1) * c(i0) * c(i1);
    // firstTerm += 10.0 * a(i1) * c(i0) * c(i1);
    // firstTerm += 5.0 * a(i0) * c(i2) * c(i2);
    // firstTerm += 10.0 * a(i2) * c(i0) * c(i2);
    // firstTerm += b(i0) * c(i1) * c(i1);
    // firstTerm += 2.0 * b(i1) * c(i0) * c(i1);
    // firstTerm += 2.0 * b(i2) * c(i0) * c(i2);
    // firstTerm /= 60.0;

    // //Second term;
    // double secondTerm = a(i1) * (c(i1) * (4.0 * a(i0) + b(i0))) + b(i1) * c(i0);
    // secondTerm += a(i2) * (c(i2) * (4.0 * a(i0) + b(i0)) + b(i2) * c(0));
    // secondTerm /= 12.0;

    // //Third Term
    // double thirdTerm = 5.0 * a(i0) * b(i1) * c(i1);
    // thirdTerm += 5.0 * a(i0) * b(i2) * c(i2);
    // thirdTerm += 10.0 * a(i1) * a(i1) * c(i0);
    // thirdTerm += 10.0 * a(i2) * a(i2) * c(i0);
    // thirdTerm += b(i1) * b(i1) * c(i0);
    // thirdTerm += b(i2) * b(i2) * c(i0);
    // thirdTerm += 2.0 * b(i0) * b(i1) * c(i1);
    // thirdTerm += 2.0 * b(i0) * b(i2) * c(i2);
    // thirdTerm /= 60.0;

    // //Fourth Term
    // double fourthTerm = 4.0 * (a(i1) * b(i1) + a(i2) * b(i2));
    // fourthTerm += 6.0 * (a(i1) * a(i1) + a(i2) * a(i2));
    // fourthTerm += b(i1) * b(i1) + b(i2) * b(i2);
    // fourthTerm *= 5.0 * a(i0);

    // double fourthTerm2 = 10.0 * (a(i1) * b(i1) + a(i2) * b(i2) + a(i1) * a(i1) + a(i2) * a(i2));
    // fourthTerm2 += 3.0 * (b(i1) * b(i1) + b(i2) * b(i2));
    // fourthTerm2 *= b(i0);

    // fourthTerm += fourthTerm2;
    // fourthTerm /= 60.0;

    // //Fifth Term
    // double fifthTerm = (1.0/20.0) * c(i0) * (c(i1) * c(i1) + c(i2) * c(i2));

    double firstTerm = 1.0/12.0*a(i1)*(2.0*a(i1)*(b(i0) + c(i0)) + a(i0)*(6.0*a(i1) + 4.0*(b(i1) + c(i1))) + b(i1)*(2.0*b(i0) + c(i0)) + c(i1)*(b(i0) + 2.0*c(i0)));
    double secondTerm = 1.0/60.0*(5.0*a(i0)*(b(i1)*c(i1) + b(i1)*b(i1) + c(i1)*c(i1)) + b(i0)*(2.0*b(i1)*c(i1) + 3.0*b(i1)*b(i1) + c(i1)*c(i1)) + c(i0)*(2.0*b(i1)*c(i1) + b(i1)*b(i1) + 3.0*c(i1)*c(i1)));
    double thirdTerm = 1.0/12.0*a(i2)*(2.0*a(i2)*(b(i0) + c(i0)) + a(i0)*(6.0*a(i2) + 4.0*(b(i2) + c(i2))) + b(i2)*(2.0*b(i0) + c(i0)) + c(i2)*(b(i0) + 2.0*c(i0)));
    double fourthTerm = 1.0/60.0*(5.0*a(i0)*(b(i2)*c(i2) + b(i2)*b(i2) + c(i2)*c(i2)) + b(i0)*(2.0*b(i2)*c(i2) + 3.0*b(i2)*b(i2) + c(i2)*c(i2)) + c(i0)*(2.0*b(i2)*c(i2) + b(i2)*b(i2) + 3.0*c(i2)*c(i2)));
    
    return (firstTerm + secondTerm + thirdTerm + fourthTerm) * (b-a).cross(c-a)(type);
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

        inertiaTensor(0,0) += computeXSquaredYSquaredTerm(T0, T1, T2, 0);
        
    }
    
    cout << inertiaTensor <<endl;
    return inertiaTensor;
}

}
