#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <cfloat>
#include <tuple>
#include <omp.h>
#include <list>
#include <chrono>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"
#include "stb_image.h"

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0.0, 1.0);

// -----------------------------------------------------------------------------------------------------------------------------------------------
// VECTOR CLASS ----------------------------------------------------------------------------------------------------------------------------------

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    Vector normalize() {
        double n = norm();
        Vector temp;
        temp[0] = data[0]/n;
        temp[1] = data[1]/n;
        temp[2] = data[2]/n;
        return temp;
    }
    // We implement this operator to make code cleaner later on...

    void operator+=(const Vector& b) {
    data[0] += b[0];
    data[1] += b[1];
    data[2] += b[2];
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

// Random Cosine for Indirect Lighting.
Vector generateRandomDirection(const Vector &normal) {

    // Generate two random numbers
    double rand1 = uniform(engine);
    double rand2 = uniform(engine);

    // Convert uniform random numbers to spherical coordinates
    double azimuth = 2. * M_PI * rand1;
    double z = sqrt(rand2);
    double x = sqrt(1 - rand2) * cos(azimuth);
    double y = sqrt(1 - rand2) * sin(azimuth);

    // Identify the axis with the smallest component to minimize floating-point errors
    double least = std::numeric_limits<double>::max();
    int smallestAxis = 0;
    for (int i = 0; i < 3; i++) {
        if (fabs(normal[i]) < least) {
            least = fabs(normal[i]);
            smallestAxis = i;
        }
    }

    // Create orthogonal vectors
    Vector helperVector;
    switch (smallestAxis) {
        case 0:
            helperVector = Vector(0., normal[2], -normal[1]).normalize();
            break;
        case 1:
            helperVector = Vector(normal[2], 0., -normal[0]).normalize();
            break;
        case 2:
            helperVector = Vector(normal[1], -normal[0], 0.).normalize();
            break;
    }

    Vector orthogonalVector = cross(normal, helperVector);

    // Return the random direction vector
    return helperVector * x + orthogonalVector * y + normal * z;
}



// -----------------------------------------------------------------------------------------------------------------------------------------------
// RAYTRACER CLASS ----------------------------------------------------------------------------------------------------------------------------------

class Ray {
    public:
        Vector O;
        Vector u;
        explicit Ray(Vector origin, Vector direction) {
            O = origin;
            u = direction;
        }
};


// -----------------------------------------------------------------------------------------------------------------------------------------------
// INTERSECTION ----------------------------------------------------------------------------------------------------------------------------------
// this is the general structure for intersections in the raytracer

struct Intersection {
    
    Vector P;
    Vector N;
    Vector color;
    double t = 0.;
    double refractiveIndex;
    bool isMirror;
    bool intersects = false;

    // Constructor using member initializer list
    Intersection(const Vector& color = Vector(0., 0., 0.), double refractiveIndex = 1.0, bool isMirror = false)
    : color(color), refractiveIndex(refractiveIndex), isMirror(isMirror) {}

};


// -----------------------------------------------------------------------------------------------------------------------------------------------
// GEOMETRY CLASS ----------------------------------------------------------------------------------------------------------------------------------

class Geometry {
public:
    virtual Intersection intersect(const Ray& ray) = 0;
};


// -----------------------------------------------------------------------------------------------------------------------------------------------
// SPHERE CLASS ----------------------------------------------------------------------------------------------------------------------------------

class Sphere : public Geometry {
    private:
        Vector centerPosition_;  // Center position of the sphere
        Vector surfaceColor_;    // Color of the sphere's surface
        double sphereRadius_;    // Radius of the sphere
        double opticalDensity_;  // Refractive index used for optical effects
        bool reflective_;        // Reflectivity flag
        bool hollow_;            // Hollow flag indicating if the sphere is hollow or solid

    public:
        // Constructor initializes sphere properties
        Sphere(const Vector& center, double radius, const Vector& color, bool reflective = false, double refractiveIndex = 1.0, bool hollow = false)
        : centerPosition_(center), sphereRadius_(radius), surfaceColor_(color), opticalDensity_(refractiveIndex), reflective_(reflective), hollow_(hollow) {}

        // Method to compute intersection of a ray with the sphere
        Intersection intersect(const Ray& ray) override {
            Intersection result(surfaceColor_, opticalDensity_, reflective_);

            // Calculate vector from ray origin to sphere center
            Vector rayToCenter = ray.O - this->centerPosition_;
            // Compute discriminant to check for intersection
            double discriminant = pow(dot(ray.u, rayToCenter), 2) - (dot(rayToCenter, rayToCenter) - pow(sphereRadius_, 2));

            // If the discriminant is non-negative, intersections are possible
            if (discriminant >= 0.) {
                double entryPoint = -dot(ray.u, rayToCenter) - sqrt(discriminant);
                double exitPoint = -dot(ray.u, rayToCenter) + sqrt(discriminant);
                result.t = (entryPoint > 0) ? entryPoint : ((exitPoint > 0) ? exitPoint : 0.0);
                result.intersects = (exitPoint < 0.) ? false : true;
            }

            // Compute the intersection point and normal at the intersection
            result.P = ray.O + (result.t * ray.u);
            result.N = (result.P - this->centerPosition_).normalize();
            // Adjust normal for hollow spheres
            result.N = (this->hollow_) ? - 1. * result.N : result.N;
            return result;
        }
};



// -----------------------------------------------------------------------------------------------------------------------------------------------
// BOUNDINGBOX CLASS ----------------------------------------------------------------------------------------------------------------------------------

class BoundingBox {
public:
    Vector B_min;
    Vector B_max;

    explicit BoundingBox(Vector min = Vector(), Vector max = Vector()) {
        B_min = min;
        B_max = max;
    }
};

// -----------------------------------------------------------------------------------------------------------------------------------------------
// TRIANGLEINDICES CLASS ----------------------------------------------------------------------------------------------------------------------------------

// from the lecture
class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

// -----------------------------------------------------------------------------------------------------------------------------------------------
// NODE STRUCTURE ----------------------------------------------------------------------------------------------------------------------------------

struct Node {
    BoundingBox bbox;
    int startingTriangle;
    int endingTriangle;
    Node* leftChild;
    Node* rightChild;
};

// -----------------------------------------------------------------------------------------------------------------------------------------------
// TRIANGLE MESH ----------------------------------------------------------------------------------------------------------------------------------

class TriangleMesh : public Geometry {
    double scalingFactor;
    Vector translation;
    Vector color;
    double refractiveIndex;
    bool isMirror;

public:

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    Node* root;

     ~TriangleMesh() {
        delete root; // Ensure proper cleanup
    }
    TriangleMesh(double scaling_factor, const Vector& translation, const Vector& color = Vector(0., 0., 0.), double refractiveIndex = 1.0, bool isMirror = false)
    : scalingFactor(scaling_factor), translation(translation), color(color), refractiveIndex(refractiveIndex), isMirror(isMirror), root(new Node) {}
    

    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);

        this->buildBVH(root, 0, indices.size());
    }


  
    BoundingBox calculateBoundingBox(int firstTriangle, int lastTriangle) {
        // Initialize bounds to extreme values to capture the correct min and max through comparisons
        double minimumX = DBL_MAX;
        double maximumX = DBL_MIN;

        double minimumY = DBL_MAX;
        double maximumY = DBL_MIN;

        double minimumZ = DBL_MAX;
        double maximumZ = DBL_MIN;

        // Loop over the specified range of triangle indices
        for (int index = firstTriangle; index < lastTriangle; index++) {
            // Get vertices of the current triangle
            auto triangleVertices = {
                this->vertices[this->indices[index].vtxi],
                this->vertices[this->indices[index].vtxj],
                this->vertices[this->indices[index].vtxk]
            };
            // Calculate the bounding box by checking each vertex
            for (const auto& vertex : triangleVertices) {
                // Apply transformation to the vertex
                Vector transformedVertex = scalingFactor * vertex + translation;
                // Update minimum and maximum coordinates
                minimumX = std::min(minimumX, transformedVertex[0]);
                maximumX = std::max(maximumX, transformedVertex[0]);
                minimumY = std::min(minimumY, transformedVertex[1]);
                maximumY = std::max(maximumY, transformedVertex[1]);
                minimumZ = std::min(minimumZ, transformedVertex[2]);
                maximumZ = std::max(maximumZ, transformedVertex[2]);
            }
        }

        // Return the bounding box constructed from the min and max coordinates
        return BoundingBox(Vector(minimumX, minimumY, minimumZ), Vector(maximumX, maximumY, maximumZ));
    }


    // Calculate minimal and maximal t values from bounding box intersections
    std::tuple<double, double> calculateIntersectionParams(Vector normal, Vector boxMin, Vector boxMax, const Ray& ray) {
        double entryTime = dot(boxMin - ray.O, normal) / dot(ray.u, normal);
        double exitTime = dot(boxMax - ray.O, normal) / dot(ray.u, normal);
        double minTime = std::min(entryTime, exitTime);
        double maxTime = std::max(entryTime, exitTime);
        return std::make_tuple(minTime, maxTime);
    }

    // Determines if the ray intersects with the bounding box
    bool intersectsBoundingBox(const Ray& ray, BoundingBox box, double* t) {
        // Intersection tests for each axis
        auto [minX, maxX] = calculateIntersectionParams(Vector(1, 0, 0), box.B_min, box.B_max, ray);
        auto [minY, maxY] = calculateIntersectionParams(Vector(0, 1, 0), box.B_min, box.B_max, ray);
        auto [minZ, maxZ] = calculateIntersectionParams(Vector(0, 0, 1), box.B_min, box.B_max, ray);

        double minIntersection = std::min(maxX, std::min(maxY, maxZ));
        double maxIntersection = std::max(minX, std::max(minY, minZ));

        *t = (minIntersection > maxIntersection) ? maxIntersection : *t;

        return (minIntersection > maxIntersection);
    }

    // Calculate the center of gravity for the triangle specified by index i
    Vector computeCentroid(int index) {
        Vector pointA = (scalingFactor * this->vertices[this->indices[index].vtxi]) + translation;
        Vector pointB = (scalingFactor * this->vertices[this->indices[index].vtxj]) + translation;
        Vector pointC = (scalingFactor * this->vertices[this->indices[index].vtxk]) + translation;
        return (pointA + pointB + pointC) / 3.0;
    }

    // Build the bounding volume hierarchy
    void buildBVH(Node *node, int start, int end) {
        node->bbox = calculateBoundingBox(start, end);
        node->startingTriangle = start;
        node->endingTriangle = end;

        Vector diagonal = node->bbox.B_max - node->bbox.B_min;
        Vector midpoint = node->bbox.B_min + diagonal * 0.5;
        int axis = (fabs(diagonal[0]) > fabs(diagonal[1]) && fabs(diagonal[0]) > fabs(diagonal[2])) ? 0 :
                (fabs(diagonal[1]) > fabs(diagonal[2])) ? 1 : 2;

        int partitionIndex = start;
        for (int i = start; i < end; i++) {
            Vector center = computeCentroid(i);
            if (center[axis] < midpoint[axis]) {
                std::swap(indices[i], indices[partitionIndex]);
                partitionIndex++;
            }
        }

        // Check if we can divide further or terminate recursion
        if (partitionIndex == start || partitionIndex == end) return;
        node->leftChild = new Node();
        node->rightChild = new Node();
        buildBVH(node->leftChild, start, partitionIndex);
        buildBVH(node->rightChild, partitionIndex, end);
    }
    // Determines if the ray intersects with any of the objects in the scene
    Intersection intersect(const Ray &ray) override {

        Intersection collisionResult(this->color, this->refractiveIndex, this->isMirror);
        double nearestDistance = DBL_MAX;
        double intersectionDistance;

        // Early exit if no intersection with the bounding box of the root node
        if (!intersectsBoundingBox(ray, this->root->bbox, &intersectionDistance))
            return Intersection();  // Return an empty collision result if no intersection

        std::list<Node*> nodesToCheck;
        nodesToCheck.push_front(root);
        while (!nodesToCheck.empty()) {
            Node* currentNode = nodesToCheck.back();
            nodesToCheck.pop_back();

            // Recursively check child nodes if they intersect with the ray
            if (currentNode->leftChild) {
                if (intersectsBoundingBox(ray, currentNode->leftChild->bbox, &intersectionDistance))
                    if (intersectionDistance < nearestDistance)
                        nodesToCheck.push_back(currentNode->leftChild);

                if (intersectsBoundingBox(ray, currentNode->rightChild->bbox, &intersectionDistance))
                    if (intersectionDistance < nearestDistance)
                        nodesToCheck.push_back(currentNode->rightChild);
            } else {
                // Leaf node: check all triangles within this node
                Vector vertexA, vertexB, vertexC, edge1, edge2, normal;
                for (int i = currentNode->startingTriangle; i < currentNode->endingTriangle; i++) {
                    TriangleIndices indices = this->indices[i];
                    vertexA = scalingFactor * vertices[indices.vtxi] + translation;
                    vertexB = scalingFactor * vertices[indices.vtxj] + translation;
                    vertexC = scalingFactor * vertices[indices.vtxk] + translation;
                    edge1 = vertexB - vertexA;
                    edge2 = vertexC - vertexA;
                    normal = cross(edge1, edge2);

                    double beta = dot(cross(vertexA - ray.O, ray.u), edge2) / dot(ray.u, normal);
                    double gamma = -dot(cross(vertexA - ray.O, ray.u), edge1) / dot(ray.u, normal);
                    double alpha = 1.0 - beta - gamma;
                    double t = dot(vertexA - ray.O, normal) / dot(ray.u, normal);

                    if (alpha >= 0 && beta >= 0 && gamma >= 0 && t > 0 && t < nearestDistance) {
                        nearestDistance = t;
                        collisionResult.intersects = true;
                        collisionResult.t = t;
                        collisionResult.P = vertexA + beta * edge1 + gamma * edge2;
                        collisionResult.N = normal.normalize();
                    }
                }
            }
        }
        return collisionResult;
    }

};

// -----------------------------------------------------------------------------------------------------------------------------------------------
// SCENE CLASS ----------------------------------------------------------------------------------------------------------------------------------

class Scene {
private:
    std::vector<Geometry*> geometries_;  // Stores all geometrical objects in the scene
    Vector S_;  // Light source position
    double intensity_ = 100000;  // Light source intensity

public:
    // Constructor initializes the light source
    explicit Scene(Vector light) { S_ = light; }

    // Adds a geometry object to the scene
    void addGeometry(Geometry* geometry) { geometries_.push_back(geometry); }

    // Traces a ray through the scene and finds the closest intersection
    Intersection traceRay(const Ray& ray) {
        Intersection nearestIntersection, currentIntersection;
        double closestDistance = std::numeric_limits<double>::infinity();
        for (auto& geometry : geometries_) {
            currentIntersection = geometry->intersect(ray);
            if (currentIntersection.intersects && currentIntersection.t < closestDistance) {
                closestDistance = currentIntersection.t;
                nearestIntersection = currentIntersection;
            }
        }
        return nearestIntersection;
    }

    // Calculates the color at the ray's point of intersection using recursive ray tracing
    Vector getColor(const Ray& ray, int depth) {
        if (depth < 0) return Vector(0., 0., 0.);  // Base case for recursion

        Intersection intersection = traceRay(ray);
        Vector colorOutput(0., 0., 0.);

        if (intersection.intersects) {
            Vector adjustedPosition = intersection.P + (1e-10 * intersection.N);
            Vector normal = intersection.N;

            // Handling reflections
            if (intersection.isMirror) {
                Vector reflectedDirection = ray.u - (2 * dot(ray.u, normal) * normal);
                Ray reflectedRay(adjustedPosition, reflectedDirection);
                return getColor(reflectedRay, depth - 1);
            }

            // Handling refractions
            if (intersection.refractiveIndex != 1.) {
                double viewDotNormal = dot(ray.u, normal);
                double refractiveRatio = (viewDotNormal > 0.) ? intersection.refractiveIndex : 1.;
                double inverseRefractiveRatio = (viewDotNormal > 0.) ? 1. : intersection.refractiveIndex;
                normal = (viewDotNormal > 0.) ? -1.*normal : normal;

                adjustedPosition = intersection.P - (1e-10 * normal);
                viewDotNormal = dot(ray.u, normal);
                double refractionCheck = 1. - pow(refractiveRatio / inverseRefractiveRatio, 2) * (1. - pow(viewDotNormal, 2));

                if (refractionCheck > 0.) {
                    Vector refractedPartT = (refractiveRatio / inverseRefractiveRatio) * (ray.u - viewDotNormal * normal);
                    Vector refractedPartN = -1.*normal * sqrt(refractionCheck);
                    Vector refractedDirection = refractedPartT + refractedPartN;
                    
                    // Implementing Frensel Reflection 
                    double fresnelReflectance = pow((refractiveRatio - inverseRefractiveRatio) / (refractiveRatio + inverseRefractiveRatio), 2);
                    double reflectance = fresnelReflectance + (1 - fresnelReflectance) * pow(1 - fabs(dot(normal, refractedDirection)), 5);
                    if (uniform(engine) < reflectance) {
                        Vector totalReflectionDirection = ray.u - (2 * dot(ray.u, intersection.N) * intersection.N);
                        Ray totalReflectedRay(adjustedPosition, totalReflectionDirection);
                        return getColor(totalReflectedRay, depth - 1);
                    } else {
                        Ray refractedRay(adjustedPosition, refractedDirection);
                        return getColor(refractedRay, depth - 1);
                    }
                } else {
                    // Handle total internal reflection
                    Vector totalInternalReflectionDirection = ray.u - (2 * dot(intersection.N, ray.u) * intersection.N);
                    Ray totalInternalReflectedRay(adjustedPosition, totalInternalReflectionDirection);
                    return getColor(totalInternalReflectedRay, depth - 1);
                }
            }

            // Calculate direct lighting using the point light model
            double lightDistance = (S_ - adjustedPosition).norm();
            Vector lightDirection = (S_ - adjustedPosition).normalize();
            Intersection lightIntersection = traceRay(Ray(S_, lightDirection * (-1)));
            double visibility = (!lightIntersection.intersects || lightIntersection.t > lightDistance) ? 1. : 0.;
            colorOutput = intensity_ / (4 * M_PI * lightDistance * lightDistance) * intersection.color / M_PI * visibility * std::max(0., dot(lightDirection, normal));

            // Add indirect lighting contribution via Monte Carlo ray tracing
            Ray randomRay(adjustedPosition, generateRandomDirection(normal));
            colorOutput += intersection.color * getColor(randomRay, depth - 1);
        }

        return colorOutput;
    }
};


void BoxMuller(double standardDev, double& x, double& y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = standardDev * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
    y = standardDev * sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
}




// -----------------------------------------------------------------------------------------------------------------------------------------------
// MAIN  ----------------------------------------------------------------------------------------------------------------------------------

int main() {

    // BEGIN timer
    auto start = std::chrono::high_resolution_clock::now();

    // Initializations

    // Width and Height
    int W = 512;
    int H = 512;

    std::vector<unsigned char> image(W*H*3, 0);
    Vector Camera = Vector(0, 0, 55);
    
    // Angle in radians
    double angle = 1.047; 
    double gamma = 2.2;
    
    // Adjust these based off Computer Performance.
    int max_bounces = 5;
    int rayPerPixel = 1;

    bool sphere = false;

    Scene scene = Scene(Vector(-10, 20, 40));

    // these are the 3 spheres to check for refraction, frensel, reflection, etc.

    Sphere S1           (Vector(20, 0, 0), 10, Vector(1., 1., 1.), false, 1, false);
    Sphere S2           (Vector(-20, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    Sphere S3           (Vector(0, 0, 0), 10, Vector(1., 1., 1.), true, 1);

    // build the walls from spheres as per lecture notes
    Sphere ceilingWall      (Vector(0, 1000, 0),    940,    Vector(0.5, 0.5, 0));
    Sphere floorWall        (Vector(0, -1000, 0),   990,    Vector(1, 1, 0));
    Sphere frontWall        (Vector(0, 0, -1000),   940,    Vector(0, 0.5, 0.5));
    Sphere backWall         (Vector(0, 0, 1000),    940,    Vector(0.5, 0, 0.9));
    Sphere leftWall         (Vector(1000, 0, 0),    940,    Vector(0.5, 0.5, 1));
    Sphere rightWall        (Vector(-1000, 0, 0),   940,    Vector(1, 0.5, 0.5));

    // we add everything to the scene 
    scene.addGeometry(&S1);
    scene.addGeometry(&S2);
    scene.addGeometry(&S3);

    scene.addGeometry(&ceilingWall);
    scene.addGeometry(&floorWall);
    scene.addGeometry(&frontWall);
    scene.addGeometry(&backWall);
    scene.addGeometry(&leftWall);
    scene.addGeometry(&rightWall);

    // add cat to the scene
    // TriangleMesh cat = TriangleMesh(0.6, Vector(0, -10, 0), Vector(1., 1., 1.));
    // cat.readOBJ("cat.obj");
    // scene.addGeometry(&cat);


    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector pixelColor = Vector(0., 0., 0.);
            double x, y;

            for (int k = 0; k < rayPerPixel; k++) {
                double u = (j + uniform(engine)) / W;
                double v = (i + uniform(engine)) / H;
                Vector Pixel = Vector(Camera[0] + (u - 0.5) * 2 * tan(angle/2) * W,
                                    Camera[1] - (v - 0.5) * 2 * tan(angle/2) * H,
                                    Camera[2] - W/(2*tan(angle/2)));
                Ray ray = Ray(Camera, (Pixel-Camera).normalize());
                pixelColor += scene.getColor(ray, max_bounces);
            }

            image[(i * W + j) * 3 + 0] = std::min(255., pow(pixelColor[0]/rayPerPixel, 1./gamma) * 255);
            image[(i * W + j) * 3 + 1] = std::min(255., pow(pixelColor[1]/rayPerPixel, 1./gamma) * 255);
            image[(i * W + j) * 3 + 2] = std::min(255., pow(pixelColor[2]/rayPerPixel, 1./gamma) * 255);
        }
    }
    // END timer
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "Rendering time: " << duration.count() << " ms." << std::endl;

    stbi_write_png("imag5.png", W, H, 3, &image[0], 0);
    return 0;
}
