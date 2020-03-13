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
        cur = cur.cwiseProduct(T1 - T0).cross(T2 - T0);
        cm += cur;
    }

    return (cm / 2.0) / volume_;
}

Eigen::Matrix3d
RigidBodyTemplate::computeInertiaTensor()
{
    Eigen::Matrix3d inertiaTensor;    
    inertiaTensor.setZero();
    // TODO: Computer Inertia Tensor
    return inertiaTensor;
}

}
