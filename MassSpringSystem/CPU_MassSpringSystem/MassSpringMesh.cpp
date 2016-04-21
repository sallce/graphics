
#include "MassSpringMesh.h"
#include "ObjMesh.h"
#include <set>

MassSpringMesh::MassSpringMesh() {
    this->n = 0;
    this->faceType = MSS_NONE;
}

MassSpringMesh::~MassSpringMesh() {}

bool ExtractEdges(std::vector<uint32_t> indices, FaceType face, std::vector<Edge>& edges) {
    if ( indices.size() == 0 ) return false;
    if ( indices.size() % face != 0 ) return false;
    if ( edges.size() != 0 ) edges.erase(edges.begin(), edges.end());

    std::set<Edge> unique_edges;

    Edge e;
    std::size_t index = 0;
    for ( std::size_t i = 0; i < indices.size() / face; i++ ) {
        for ( std::size_t x = 0; x < face; x++ ) {
            uint32_t p, q;

            if ( x == face - 1 ) {
                p = indices[index + x];
                q = indices[index];
            }
            else {
                p = indices[index + x];
                q = indices[index + x + 1];
            }

            e.p = p;
            e.q = q;
            unique_edges.insert(e);
        }

        index += face;
    }

    for ( auto i = unique_edges.begin(); i != unique_edges.end(); i++ ) {
        edges.push_back(*i);
    }

    return true;
}

bool ConstructSprings(const std::vector<Particle>& nodes, const std::vector<Edge>& edges, float ks, float kd, std::vector<Spring>& springMatrix) {
    if ( edges.size() == 0 ) return false;
    if ( ks <= 0.0f ) return false;
    if ( kd <= 0.0f ) return false;

    Spring s;
    for ( std::size_t i = 0; i < edges.size(); i++ ) {
        s.p = edges[i].p;
        s.q = edges[i].q;
        s.ks = ks;
        s.kd = kd;

        const Vector3f& pPosition = nodes[s.p].position;
        const Vector3f& qPosition = nodes[s.q].position;
        s.rl = static_cast<float>(Vector3f::Distance(pPosition, qPosition));

        springMatrix.push_back(s);
    }

    return true;
}

bool MassSpringMesh::load(const std::string& objFilename, float ks, float kd, Vector3f gravity) {
    std::shared_ptr<ObjMesh> mesh = nullptr;

    //------------------------------------------------------------------------------
    // Load OBJ Mesh.
    //------------------------------------------------------------------------------
    if ( !LoadObjMesh(objFilename, mesh) ) return false;

    if ( mesh->faces.size() == 0 ) return false;
    if ( mesh->faces[0].type != TRIANGLE && mesh->faces[0].type != QUAD ) return false;
    if ( mesh->faces[0].type == TRIANGLE ) this->faceType = MSS_TRIANGLE;
    if ( mesh->faces[0].type == QUAD ) this->faceType = MSS_QUAD;

    //------------------------------------------------------------------------------
    // Node Positions from OBJ Vertices
    //------------------------------------------------------------------------------
    std::vector<Vector3f>& vertices = mesh->vertices;
    this->nodes.resize(vertices.size());
    for ( std::size_t i = 0; i < vertices.size(); i++ ) {
        nodes[i].position = vertices[i];
        nodes[i].velocity = Vector3f::Zero();
        nodes[i].force = Vector3f::Zero();
        nodes[i].restPosition = Vector3f::Zero();
        nodes[i].restVelocity = Vector3f::Zero();
        nodes[i].restForce = gravity;
        nodes[i].mass = 1.0f;
    }

    //------------------------------------------------------------------------------
    // Face indices (also for edge/spring genereration)
    //------------------------------------------------------------------------------
    for ( std::size_t i = 0; i < mesh->faces.size(); i++ ) {
        for ( std::size_t x = 0; x < mesh->faces[i].vertexIndices.size(); x++ ) {
            indices.push_back((uint32_t)mesh->faces[i].vertexIndices[x]);
        }
    }

    //------------------------------------------------------------------------------
    // Convert unique edges to springs.
    //------------------------------------------------------------------------------
    ExtractEdges(this->indices, MSS_TRIANGLE, this->edges);
    ConstructSprings(this->nodes, this->edges, ks, kd, this->springMatrix);
    this->n = OptimalSpringMatrixDimension((uint32_t)this->springMatrix.size());

    return true;
}

void CalculateSpringForce(const Spring& spring, std::vector<Particle>& nodes) {
    uint32_t p = spring.p;
    uint32_t q = spring.q;

    Vector3f& pPos = nodes[p].position;
    Vector3f& pVel = nodes[p].velocity;
    Vector3f& pFor = nodes[p].force;
    Vector3f& qPos = nodes[q].position;
    Vector3f& qVel = nodes[q].velocity;
    Vector3f& qFor = nodes[q].force;

    float distance = static_cast<float>(Vector3f::Distance(pPos, qPos));
    float rest = spring.rl;
    float ks = spring.ks;
    float kd = spring.kd;

    float x = -(rest - distance);
    float f = ks * x;

    Vector3f pDir = (qPos - pPos).normalized();
    Vector3f qDir = (pPos - qPos).normalized();

    // Damping
    float pdvel = static_cast<float>(pDir.dot(pVel));
    float qdvel = static_cast<float>(qDir.dot(qVel));
    Vector3f pDampDir = pdvel * pDir;
    Vector3f qDampDir = qdvel * qDir;

    pFor += (f * pDir) - kd * pDampDir;
    qFor += (f * qDir) - kd * qDampDir;
}

void MassSpringMesh::update(float dt) {
    this->resetRestForces();

    //------------------------------------------------------------------------------
    // Iterate through each spring. Each spring force is calculated and contributes
    // to the internal force acting on the nodes of the system.
    //------------------------------------------------------------------------------
    for ( std::size_t i = 0; i < n; i++ ) {
        for ( std::size_t j = 0; j < n; j++ ) {
            std::size_t index = i * n + j;

            if ( index < this->springMatrix.size() ) {
                Spring& s = this->springMatrix[index];
                CalculateSpringForce(s, this->nodes);
            }
        }
    }

    //------------------------------------------------------------------------------
    // For each node, update the new velocity and position (Exp. Euler Integration)
    //------------------------------------------------------------------------------
    for ( std::size_t i = 0; i < this->nodes.size(); i++ ) {
        this->nodes[i].velocity = this->nodes[i].velocity + (dt * this->nodes[i].force * (1.0f / this->nodes[i].mass));
        this->nodes[i].position = this->nodes[i].position + (dt * this->nodes[i].velocity);
    }
}

void MassSpringMesh::resolveCollisions() {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ ) {
        if ( this->nodes[i].position.y() < 0.0f ) {
            this->nodes[i].position.y() = 0.0f;
            this->nodes[i].velocity.y() *= -1.0f;
            this->nodes[i].velocity.y() *= 0.8f;
        }
    }
}

void MassSpringMesh::translate(float dx, float dy, float dz) {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ ) {
        this->nodes[i].position.x() += dx;
        this->nodes[i].position.y() += dy;
        this->nodes[i].position.z() += dz;
    }
}

void MassSpringMesh::translate(const Vector3f& disp) {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ )
        this->nodes[i].position += disp;
}

bool MassSpringMesh::setVelocity(std::size_t index, const Vector3f& velocity) {
    if ( index >= this->nodes.size() ) return false;
    this->nodes[index].velocity = velocity;
    return true;
}

bool MassSpringMesh::setForce(std::size_t index, const Vector3f& force) {
    if ( index >= this->nodes.size() ) return false;
    this->nodes[index].force = force;
    return true;
}

bool MassSpringMesh::setRestPositionsFromCurrent() {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ )
        this->nodes[i].restPosition = this->nodes[i].position;
    return true;
}

bool MassSpringMesh::setRestVelocitiesFromCurrent() {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ )
        this->nodes[i].restVelocity = this->nodes[i].velocity;
    return true;
}

bool MassSpringMesh::setRestPosition(std::size_t index, const Vector3f& position, bool reset) {
    if ( index >= this->nodes.size() ) return false;
    this->nodes[index].restPosition = position;
    if ( reset ) this->resetRestPositions();
    return true;
}

bool MassSpringMesh::setRestVelocity(std::size_t index, const Vector3f& velocity, bool reset) {
    if ( index >= this->nodes.size() ) return false;
    this->nodes[index].restVelocity = velocity;
    if ( reset ) this->resetRestVelocities();
    return true;
}

bool MassSpringMesh::setRestForce(std::size_t index, const Vector3f& force, bool reset) {
    if ( index >= this->nodes.size() ) return false;
    this->nodes[index].restForce = force;
    if ( reset ) this->resetRestForces();
    return true;
}

void MassSpringMesh::resetRestPositions() {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ )
        this->nodes[i].position = this->nodes[i].restPosition;
}

void MassSpringMesh::resetRestVelocities() {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ )
        this->nodes[i].velocity = this->nodes[i].restVelocity;
}

void MassSpringMesh::resetRestForces() {
    for ( std::size_t i = 0; i < this->nodes.size(); i++ )
        this->nodes[i].force = this->nodes[i].restForce;
}

std::vector<Particle>& MassSpringMesh::getNodes() {
    return this->nodes;
}

const std::vector<Particle>& MassSpringMesh::getNodes() const {
    return this->nodes;
}

const std::vector<uint32_t>& MassSpringMesh::getIndices() const {
    return this->indices;
}

const std::vector<Edge>& MassSpringMesh::getEdges() const {
    return this->edges;
}

uint32_t MassSpringMesh::OptimalSpringMatrixDimension(uint32_t n) {
    return static_cast<uint32_t>(std::ceil(std::sqrt(n)));
}

