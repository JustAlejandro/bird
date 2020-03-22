#ifndef PSIM_CORE_BIRD1_BIRDSCORE_H
#define PSIM_CORE_BIRD1_BIRDSCORE_H

#include "../PhysicsCore.h"
#include "SimParameters.h"
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <memory>

namespace bird1 {

class RigidBodyTemplate;
class RigidBodyInstance;

class BirdsCore : public PhysicsCore
{
public:
    // BirdsCore() : PhysicsCore(), sceneFile_("box.scn") {}
    BirdsCore();
    ~BirdsCore();

    virtual void initSimulation() override;

    virtual
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
    getCurrentMesh() const override;

    virtual bool simulateOneStep() override;

    std::shared_ptr<SimParameters> getPointerToSimParameters()
    {
            return params_;
    }

    void clearScene();

    /*
     * The facade of adding RigidBodyTemplate and its instances.
     */
    Eigen::VectorXi
    addMesh(const std::string& file_name,
            double scale,
            const Eigen::MatrixXd& Qs);

    /*
     * add an RigidBodyInstance directly
     */
    int32_t
    addSingleInstance(std::shared_ptr<RigidBodyTemplate> rbt,
		      double density,
		      const Eigen::Vector3d &c,
		      const Eigen::Vector3d &theta,
		      const Eigen::Vector3d &cvel,
		      const Eigen::Vector3d &w);

    /*
     * addInstances:
     *    add rigid body instance from Qs.
     *    Qs should be a (N, 13) matrix and each row
     */
    Eigen::VectorXi
    addInstances(std::shared_ptr<RigidBodyTemplate> rbt,
                 const Eigen::MatrixXd& Qs);


    //
    // Similar to Goo::query*, but we are returning a modifiable object.
    // 
    // This interface change allows testing code to modify the phase directly.
    // Therefore the accumulation of numerical error during testing can be
    // eliminated by constantly synchronizing with the phase data from the reference code.
    // 
    std::shared_ptr<RigidBodyInstance>
    queryRigidBodyInstance(int32_t bid);

private:
    int32_t rigid_body_id_;
    void computeForces(Eigen::VectorXd &Fc, Eigen::VectorXd &Ftheta);

    void setM(int nbodies = 0);
    void setInertiaTensor(int nbodies = 0);
    void setupConfigVector(int nbodies);
    void setupConfigVelVector(int nbodies);

    void updateInstances(int nbodies);

    void simpleTimeIntegrator(int nbodies);

    Eigen::Vector3d newtonsMethod(const Eigen::Vector3d& oldW, const Eigen::SparseMatrix<double>& Inertia, const Eigen::SparseMatrix<double>& InertiaInv);
    Eigen::Vector3d FNewton(const Eigen::Vector3d& wGuess, const Eigen::Vector3d& oldW, const Eigen::SparseMatrix<double>& Inertia, const Eigen::SparseMatrix<double>& InertiaInv);
    Eigen::Matrix3d dFNewton(const Eigen::Vector3d& wGuess, const Eigen::Vector3d& oldW, const Eigen::SparseMatrix<double>& Inertia, const Eigen::SparseMatrix<double>& InertiaInv);

    void elasticCollision(int i, const Eigen::VectorXd& qDotGrav);

    double time_;
    std::shared_ptr<SimParameters> params_;

    std::vector<std::shared_ptr<RigidBodyTemplate>> templates_;
    std::vector<std::shared_ptr<RigidBodyInstance>> bodies_;

    std::vector<Eigen::MatrixXd> init_configurations_;
    Eigen::SparseMatrix<double> M;
    Eigen::SparseMatrix<double> MInv;
    
    Eigen::SparseMatrix<double> MInertia;
    Eigen::SparseMatrix<double> MInertiaInv;

    Eigen::VectorXd q;
    Eigen::VectorXd qPrev;

    Eigen::VectorXd qDot;
    Eigen::VectorXd qDotPrev;
};

}

#endif
