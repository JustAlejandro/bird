#ifndef PSIM_CORE_BIRD1_RIGIDBODYINSTANCE_H
#define PSIM_CORE_BIRD1_RIGIDBODYINSTANCE_H

#include <Eigen/Core>
#include <list>
#include <vector>
#include <set>

namespace bird1 {

class RigidBodyTemplate;

class RigidBodyInstance
{
public:
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate,
                      const Eigen::Vector3d &c,
		      const Eigen::Vector3d &theta,
		      const Eigen::Vector3d &cvel,
		      const Eigen::Vector3d &w,
		      double density);
    ~RigidBodyInstance();

    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;

    //Index of the objects colliding with this object
    std::set<int> collided;

    bool collisionCalculated;

    double density;

    double mass;
    
    const RigidBodyTemplate &getTemplate() const {return rbtemplate_;}
    
    int32_t bid;
private:
    const RigidBodyTemplate &rbtemplate_;
};

}

#endif // RIGIDBODYINSTANCE_H
