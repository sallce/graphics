
#ifndef MASS_SPRING_MESH_H
#define MASS_SPRING_MESH_H

#include <vector>
#include <inttypes.h>
#include "Vector3.h"

enum FaceType { MSS_NONE, MSS_TRIANGLE = 3, MSS_QUAD = 4};

struct Particle {
    Vector3f position, velocity, force;
    Vector3f restPosition, restVelocity, restForce;
    float mass;
};

struct Spring {
    uint32_t p;
    uint32_t q;
    float ks, kd, rl;
};

struct Edge {
    bool operator < ( const Edge& e) const { return (this->p * this->p + this->q * this->q) < (e.p * e.p + e.q * e.q); }
    uint32_t p, q;
};

class MassSpringMesh {
public:
    MassSpringMesh();
    ~MassSpringMesh();

    bool load(const std::string& objFilename, float ks, float kd, Vector3f gravity = Vector3f::Zero());
    void update(float dt);
    void resolveCollisions();

    void translate(float dx, float dy, float dz);
    void translate(const Vector3f& disp);

    bool setVelocity(std::size_t index, const Vector3f& velocity);
    bool setForce(std::size_t index, const Vector3f& force);

    bool setRestPositionsFromCurrent();
    bool setRestVelocitiesFromCurrent();
    bool setRestPosition(std::size_t index, const Vector3f& position, bool reset = true);
    bool setRestVelocity(std::size_t index, const Vector3f& velocity, bool reset = true);
    bool setRestForce(std::size_t index, const Vector3f& force, bool reset = true);

    void resetRestPositions();
    void resetRestVelocities();
    void resetRestForces();

    std::vector<Particle>& getNodes();
    const std::vector<Particle>& getNodes() const;
    const std::vector<uint32_t>& getIndices() const;
    const std::vector<Edge>& getEdges() const;

    static uint32_t OptimalSpringMatrixDimension(uint32_t n);

protected:
    std::vector<Particle> nodes;
    std::vector<uint32_t> indices;
    std::vector<Spring> springMatrix;
    std::vector<Edge> edges;

    uint32_t n;

    FaceType faceType;
};

#endif
