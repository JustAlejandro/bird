#include "BirdsCore.h"
#include "helper.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include <iostream>
#include <Eigen/LU> // Required for .inverse()

using namespace Eigen;
using namespace std;

namespace bird1 {

BirdsCore::BirdsCore()
{
    params_ = std::make_shared<SimParameters>();
    initSimulation();
}

BirdsCore::~BirdsCore()
{
}

void BirdsCore::setM(int nbodies){
    // Have the mass matrix filled for all Tran and Rot fields of an obj
    if (M.rows() == nbodies * 3) return;
    MatrixXd MTemp = MatrixXd::Zero(bodies_.size() * 3, bodies_.size() * 3);
    MatrixXd MInvTemp = MatrixXd::Zero(bodies_.size() * 3, bodies_.size() * 3);
    for(int i = 0; i < bodies_.size(); i++){
        for (int j = 0; j < 3; j++)
        {
            MTemp(i*3+j, i*3+j) = bodies_[i]->density * bodies_[i]->getTemplate().getVolume();
            MInvTemp(i*3+j, i*3+j) =  1.0 / MTemp(i*3+j, i*3+j);
        }
    }
    M = MTemp.sparseView();
    MInv = MInvTemp.sparseView();
}

void BirdsCore::setInertiaTensor(int nbodies){
    if(MInertia.rows() == nbodies * 3) return;
    MatrixXd MInertiaTemp = MatrixXd::Zero(bodies_.size() * 3, bodies_.size() * 3);

    for (int i = 0; i < bodies_.size(); i++)
    {
        MInertiaTemp.block(i * 3, i * 3, 3, 3) = bodies_[i]->density * bodies_[i]->getTemplate().getInertiaTensor();  
    }

    MatrixXd MInertiaInvTemp = MInertiaTemp.inverse();

    MInertia = MInertiaTemp.sparseView();
    MInertiaInv = MInertiaInvTemp.sparseView();
}

bool checkCollision(const shared_ptr<RigidBodyInstance>& c1, const shared_ptr<RigidBodyInstance>& c2, const int& i, const int& j, const bool& elastic){
    bool hit = false;

    //Inelastic never gets unstuck, so no need to check every step
    if(!elastic && c1->collided.find(j) != c1->collided.end()){
        hit = true;
    }

    //Discover collisions
    else if((c1->c - c2->c).norm() < (c1->getTemplate().getBoundingRadius() + c2->getTemplate().getBoundingRadius())){
        c1->collided.insert(j);
        for(int k : c2->collided){
            if(k != i){
                c1->collided.insert(k);
            }
        }
        
        c2->collided.insert(i);
        for(int k : c1->collided){
            if(k != j){
                c2->collided.insert(k);
            }
        }
        hit = true;
    }

    //Erase collisions if inelastic
    else if(elastic){
        c1->collided.erase(j);
        c2->collided.erase(i);
        hit = false;
    }

    c1->collisionCalculated = false;
    c2->collisionCalculated = false;

    return hit;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
BirdsCore::getCurrentMesh() const
{
    int totverts = 0;
    int totfaces = 0;
    for (const auto& rbi : bodies_)
    {
        totverts += rbi->getTemplate().getVerts().rows();
        totfaces += rbi->getTemplate().getFaces().rows();
    }

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    renderQ.resize(totverts, 3);
    renderF.resize(totfaces, 3);
    int voffset = 0;
    int foffset = 0;
    for (const auto& rbi : bodies_)
    {
        int nverts = rbi->getTemplate().getVerts().rows();
        for (int i = 0; i < nverts; i++)
            renderQ.row(voffset + i) = (rbi->c + VectorMath::rotationMatrix(rbi->theta)*rbi->getTemplate().getVerts().row(i).transpose()).transpose();
        int nfaces = rbi->getTemplate().getFaces().rows();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                renderF(foffset + i, j) = rbi->getTemplate().getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }
    // std::cerr << __func__ << " Nbodies " << bodies_.size() << std::endl;
    return std::make_tuple(renderQ, renderF, renderC);
}


void BirdsCore::initSimulation()
{
    rigid_body_id_ = 0;
    time_ = 0;
    setM();
    setInertiaTensor();
}

void BirdsCore::computeForces(VectorXd &Fc, VectorXd &Ftheta)
{
    Fc = VectorXd::Zero(bodies_.size() * 3);
    Ftheta = VectorXd::Zero(bodies_.size() * 3);
    //Only Gravity for now
    if(params_->gravityEnabled || params_->inelasticEnabled || params_->elasticEnabled){
        for (int i = 0; i < bodies_.size(); i++)
        {
            for (int j = i+1; j < bodies_.size(); j++)
            {
                shared_ptr<RigidBodyInstance> c1 = bodies_[i];
                shared_ptr<RigidBodyInstance> c2 = bodies_[j];

                bool hit = checkCollision(c1, c2, i, j, params_->elasticEnabled);

                if(params_->gravityEnabled && !(hit && (params_->elasticEnabled || params_->inelasticEnabled))){
                    Vector3d diff = q.segment(i * 6, 3) - q.segment(j * 6, 3);
                    
                    Vector3d grav = params_->gravityG * c1->mass * c2->mass
                        *(1.0 / diff.squaredNorm()) * diff.normalized();
                    Fc.segment(i * 3, 3) += grav;
                    Fc.segment(j * 3, 3) -= grav;
                }
            }
        }
    }
}

Vector3d BirdsCore::FNewton(const Vector3d& wGuess, const Vector3d& oldW, const SparseMatrix<double>& Inertia, const SparseMatrix<double>& InertiaInv) {
    return wGuess - (oldW.transpose() * Inertia * VectorMath::TMatrix(params_->timeStep * oldW).inverse()
        * VectorMath::TMatrix(-params_->timeStep * wGuess) * InertiaInv).transpose();
}

Matrix3d BirdsCore::dFNewton(const Vector3d& wGuess, const Vector3d& oldW, const SparseMatrix<double>& Inertia, const SparseMatrix<double>& InertiaInv){
    return Matrix3d::Identity(3,3);
}

//Computes Newton's for a single object
Vector3d BirdsCore::newtonsMethod(const Vector3d& oldW, const SparseMatrix<double>& Inertia, const SparseMatrix<double>& InertiaInv){
    Vector3d wGuess = oldW;

    int iter = 0;
    while (FNewton(wGuess, oldW, Inertia, InertiaInv).squaredNorm() > (params_->NewtonTolerance*params_->NewtonTolerance) && iter < params_->NewtonMaxIters){
        SparseMatrix<double> df = dFNewton(wGuess, oldW, Inertia, InertiaInv).sparseView();
        //Eigen::SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver(df);

        //solver.analyzePattern(df);
        //solver.factorize(df);

        //wGuess += solver.solve(-FNewton(wGuess, oldW, Inertia, InertiaInv));
        wGuess += -FNewton(wGuess, oldW, Inertia, InertiaInv);
        
        iter++;
    }
    return wGuess;
}
void BirdsCore::setupConfigVector(int nbodies){
    qPrev = q;
    q = VectorXd(6 * nbodies);

    for (int i = 0; i < nbodies; i++)
    {
        q.segment(i * 6, 3) = bodies_[i]->c;
        q.segment(i * 6 + 3, 3) = bodies_[i]->theta;
    }

    if(q.size() != qPrev.size()) {
        qPrev = q;
    }
}

void BirdsCore::setupConfigVelVector(int nbodies){
    qDotPrev = qDot;
    qDot = VectorXd(6 * nbodies);

    for (int i = 0; i < nbodies; i++)
    {
        qDot.segment(i * 6, 3) = bodies_[i]->cvel;
        qDot.segment(i * 6 + 3, 3) = bodies_[i]->w;
    }

    if(qDot.size() != qDotPrev.size()) {
        qDotPrev = qDot;
    }
}

//Updates the config and velocity of all RigidBodyInstance(s) in scene
void BirdsCore::updateInstances(int nbodies){
    for (int i = 0; i < nbodies; i++)
    {
        shared_ptr<RigidBodyInstance> b = bodies_[i];
        b->c = q.segment(i * 6, 3);
        b->theta = q.segment(i * 6 + 3, 3);
        b->cvel = qDot.segment(i * 6, 3);
        b->w = qDot.segment(i*6+3, 3);
    }
    
}

void BirdsCore::simpleTimeIntegrator(int nbodies){

    for (int i = 0; i < nbodies; i++)
    {
        shared_ptr<RigidBodyInstance> b = bodies_[i];
        q.segment(i * 6, 3) = b->c +  params_->timeStep * b->cvel;
        q.segment(i * 6 + 3, 3) = VectorMath::axisAngle(
            VectorMath::rotationMatrix(b->theta)
            * VectorMath::rotationMatrix(params_->timeStep * b->w));
        
        //No Potential here since gravity not in yet
        qDot.segment(i * 6, 3) = qDotPrev.segment(i * 6, 3);

        qDot.segment(i*6+3, 3) = newtonsMethod(qDot.segment(i*6+3, 3), MInertia.block(i*3,i*3,3,3), MInertiaInv.block(i*3,i*3,3,3));
    }

    VectorXd transForce, rotForce;
    computeForces(transForce, rotForce);

    //Precalc Q with Gravity if no collisions occur
    VectorXd qDotGrav = qDot;
    for(int i; params_->gravityEnabled && i < nbodies; i++){
        qDotGrav.segment(i*6, 3) += -params_->timeStep * MInv.block(i*3, i*3, 3, 3) * transForce.segment(i * 3, 3);
    }

    for (int i = 0; i < nbodies; i++)
    {
        if(!params_->inelasticEnabled || !params_->elasticEnabled || !bodies_[i]->collisionCalculated){
            if(params_->inelasticEnabled && bodies_[i]->collided.size() != 0 && !params_->elasticEnabled){
                //Momentum calc
                Vector3d momentumSum = bodies_[i]->mass * qDot.segment(i*6, 3);
                double massSum = bodies_[i]->mass;

                Vector3d totalGravityForce = transForce.segment(i * 3, 3);
                Matrix3d MTotal = M.block(i*3, i*3, 3, 3);

                for(int j : bodies_[i]->collided) {
                    momentumSum += bodies_[j]->mass * qDot.segment(j*6, 3);
                    massSum += bodies_[j]->mass;
                    totalGravityForce += transForce.segment(j * 3, 3);
                    MTotal += M.block(j*3, j*3, 3, 3);
                }
                MTotal = MTotal.inverse();
                qDot.segment(i*6, 3) = momentumSum / massSum - params_->timeStep * MTotal * totalGravityForce;
                bodies_[i]->collisionCalculated = true;
                
                for(int j : bodies_[i]->collided){
                    qDot.segment(j*6, 3) = momentumSum / massSum - params_->timeStep * MTotal * totalGravityForce;
                    bodies_[j]->collisionCalculated = true;
                }
            }
            else if(params_->elasticEnabled && bodies_[i]->collided.size() != 0){
                elasticCollision(i, qDotGrav);
            }
            else{
                qDot.segment(i*6, 3) = qDotGrav.segment(i*6, 3);
            }
        }
    }
}

void BirdsCore::elasticCollision(int i, const VectorXd& qDotGrav){
    Vector3d totalVelocity = qDotGrav.segment(i*6, 3);
    Vector3d iVel = totalVelocity;

    for(int j : bodies_[i]->collided){
        //Normal to the collision plane
        Vector3d norm = q.segment(j * 6, 3) - q.segment(i * 6, 3);

        Vector3d iVel = totalVelocity;
        Vector3d jVel = qDotGrav.segment(j*6, 3);
        //Project the forces
        Vector3d projIVel = (iVel.dot(norm) / norm.squaredNorm()) * norm;
        Vector3d projJVel = (jVel.dot(norm) / norm.squaredNorm()) * norm;
        
        Vector3d tanIVel = iVel - projIVel;

        double totalMass = bodies_[i]->mass + bodies_[j]->mass;
        double diffMass = bodies_[i]->mass - bodies_[j]->mass;
        totalVelocity = (diffMass/totalMass)*projIVel
            + ((2*bodies_[j]->mass)/totalMass) * projJVel;
        
        totalVelocity += tanIVel;

        if((q.segment(i * 6, 3) - q.segment(j * 6, 3)).norm() < (bodies_[i]->getTemplate().getBoundingRadius() + bodies_[j]->getTemplate().getBoundingRadius())){
            q.segment(j * 6, 3) = (bodies_[i]->getTemplate().getBoundingRadius() + bodies_[j]->getTemplate().getBoundingRadius()) * norm.normalized()
                + q.segment(i*6,3);
        }
    }
    qDot.segment(i*6, 3) = totalVelocity;
}

bool BirdsCore::simulateOneStep()
{
    time_ += params_->timeStep;
    int nbodies = (int)bodies_.size();

    // TODO: Implement Time Integrator here
    setM(nbodies);
    setInertiaTensor(nbodies);
    setupConfigVector(nbodies);
    setupConfigVelVector(nbodies);

    simpleTimeIntegrator(nbodies);

    updateInstances(nbodies);

    return false;
}

void
BirdsCore::clearScene()
{
    bodies_.clear();
    templates_.clear();
    init_configurations_.clear();
}

Eigen::VectorXi
BirdsCore::addMesh(const std::string& file_name,
                   double scale,
                   const Eigen::MatrixXd& Qs)
{
    auto tup = loadOBJ(file_name);
    templates_.emplace_back(new RigidBodyTemplate(std::get<0>(tup), std::get<1>(tup), scale));

    auto rbt = templates_.back();
    init_configurations_.emplace_back(Qs);
    return addInstances(rbt, Qs);
}

std::shared_ptr<RigidBodyInstance>
BirdsCore::queryRigidBodyInstance(int32_t bid)
{
    for (auto& b : bodies_)
        if (b->bid == bid)
            return b;
    return std::shared_ptr<RigidBodyInstance>(nullptr);
}

int32_t
BirdsCore::addSingleInstance(std::shared_ptr<RigidBodyTemplate> rbt,
                             double density,
                             const Eigen::Vector3d &c,
                             const Eigen::Vector3d &theta,
                             const Eigen::Vector3d &cvel,
                             const Eigen::Vector3d &w)
{
    bodies_.emplace_back(new RigidBodyInstance(*rbt, c, theta, cvel, w, density));
    bodies_.back()->bid = rigid_body_id_++;
    return bodies_.back()->bid;
}

Eigen::VectorXi
BirdsCore::addInstances(std::shared_ptr<RigidBodyTemplate> rbt,
                        const Eigen::MatrixXd& Qs)
{
    Eigen::VectorXi ret;
    ret.resize(Qs.rows());
    for (int i = 0; i < Qs.rows(); i++) {
        double density;
        Eigen::Vector3d c, cvel, theta, w;
        density = Qs(i, 0);
        int base = 1;
        c << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        cvel << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        theta << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        w << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        ret(i) = addSingleInstance(rbt, density, c, theta, cvel, w);
    }
    return ret;
}

}