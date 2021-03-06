#ifndef PSYM_CORE_GOO2_SCENEOBJECTS_H
#define PSYM_CORE_GOO2_SCENEOBJECTS_H

#include "SimParameters.h"
#include <Eigen/Core>
#include <set>

namespace goo2 {

struct Particle {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Particle(Eigen::Vector2d pos, double mass, bool isFixed, bool isInert)
        : pos(pos)
        , mass(mass)
        , fixed(isFixed)
        , inert(isInert)
    {
        vel.setZero();
    }

    Eigen::Vector2d pos;
    Eigen::Vector2d vel;
    double mass;
    bool fixed;
    bool inert;

    int32_t uid = -1;
};

struct Connector {
public:
    Connector(int p1, int p2, double mass)
        : p1(p1)
        , p2(p2)
        , mass(mass)
    {
    }
    virtual ~Connector() {}

    virtual SimParameters::ConnectorType getType() = 0;

    int p1;
    int p2;
    double mass;

    std::set<int> associatedBendingStencils;
};

struct Spring : public Connector {
public:
    Spring(int p1, int p2, double mass, double stiffness, double restlen, bool canSnap)
        : Connector(p1, p2, mass)
        , stiffness(stiffness)
        , restlen(restlen)
        , canSnap(canSnap)
    {
    }

    virtual SimParameters::ConnectorType getType()
    {
        return SimParameters::CT_SPRING;
    }

    double stiffness;
    double restlen;
    bool canSnap;
};

struct RigidRod : public Connector {
public:
    RigidRod(int p1, int p2, double mass, double length)
        : Connector(p1, p2, mass)
        , lambda(0)
        , length(length)
    {
    }

    virtual SimParameters::ConnectorType getType()
    {
        return SimParameters::CT_RIGIDROD;
    }

    double lambda;
    double length;
};

struct Saw {
public:
    Saw(Eigen::Vector2d pos, double radius)
        : pos(pos)
        , radius(radius)
    {
    }

    Eigen::Vector2d pos;
    double radius;
};

struct BendingStencil {
public:
    BendingStencil(int p1, int p2, int p3, double kb)
        : p1(p1)
        , p2(p2)
        , p3(p3)
        , kb(kb)
    {
    }

    int p1, p2, p3;
    double kb;
};

}

#endif
