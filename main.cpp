#define _CRT_SECURE_NO_WARNINGS 1

#include <cmath>
#include <vector>
#include <limits>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <list>
// #include <omp.h>

// static std::default_random_engine engine[8];
static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);

double sqr(double x) {return x*x;}

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
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }

    void operator+=(const Vector& b) {
        data[0] += b[0];
        data[1] += b[1];
        data[2] += b[2];
    }

    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

class Ray {
public:
    Ray(const Vector &O, const Vector &u) : origin(O), direction(u) {};

    Vector origin;
    Vector direction;
};

class Object{
    public:
        Object(Vector rho, bool mirror) : rho(rho), mirror(mirror) {};
        virtual bool intersect(const Ray &r, Intersection &I) const = 0;
        Vector rho;
        bool mirror;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class BVHNode {
public:
    BoundingBox bounding_box;
    BVHNode* left;
    BVHNode* right;
    size_t start;
    size_t end;

    BVHNode() : left(nullptr), right(nullptr), start(-1), end(-1) {}

    ~BVHNode() {
        delete left;
        delete right;
    }
};

class Intersection {
public:
    Intersection(Vector P = Vector(0, 0, 0), Vector N = Vector(0, 0, 0), double t = std::numeric_limits<double>::max(), Vector rho = Vector(0, 0, 0))
        : P(P), N(N), t(t), rho(rho) {};

    Vector P;
    Vector N;
    double t;
    Vector rho;
    bool isMirror;
};

class BoundingBox {
    public:

        BoundingBox() {
            m = Vector(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
            M = Vector(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
        }

        bool intersect(const Ray& r) const {
            double tx1 = (m[0] - r.origin[0]) / r.direction[0];
            double tx2 = (M[0] - r.origin[0]) / r.direction[0];
            double txmin = std::min(tx1, tx2);
            double txmax = std::max(tx1, tx2);

            double ty1 = (m[1] - r.origin[1]) / r.direction[1];
            double ty2 = (M[1] - r.origin[1]) / r.direction[1];
            double tymin = std::min(ty1, ty2);
            double tymax = std::max(ty1, ty2);

            double tz1 = (m[2] - r.origin[2]) / r.direction[2];
            double tz2 = (M[2] - r.origin[2]) / r.direction[2];
            double tzmin = std::min(tz1, tz2);
            double tzmax = std::max(tz1, tz2);

            double actualT = std::max(txmin, std::max(tymin, tzmin));
            double furthestT = std::min(txmax, std::min(tymax, tzmax));
            if (furthestT > 0 && furthestT > actualT) return true;
            return false;
        }

        void expand(const Vector& v) {
            for (int i = 0; i < 3; i++) {
                m[i] = std::min(m[i], v[i]);
                M[i] = std::max(M[i], v[i]);
            }
        }

        int longestAxis() const {
            Vector d = M - m; // Extent of the bounding box along each axis
            if (d[0] > d[1] && d[0] > d[2])
                return 0;
            else if (d[1] > d[2])
                return 1;
            else
                return 2;
        }

        Vector m, M;
};


Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector random_cos(const Vector& N) {
    // int current_thread = omp_get_thread_num();

    // double r1 = uniform(engine[current_thread]);
    // double r2 = uniform(engine[current_thread]);
    double r1 = uniform(engine);
    double r2 = uniform(engine);

    double x = sqrt(1 - r1) * cos(2. * M_PI * r2);
    double y = sqrt(1 - r1) * sin(2. * M_PI * r2);
    double z = sqrt(r1);

    Vector T1;

    if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])) {
        T1 = Vector(0, N[2], -N[1]);
    }
    else {
        if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2])) {
            T1 = Vector(N[2], 0, -N[0]);
        }
        else {
            T1 = Vector(N[1], -N[0], 0);
        }
    }
    T1.normalize();
    Vector T2 = cross(N, T1);

    return x * T1 + y * T2 + z * N;
}

class Ray {
public:
    Ray(const Vector& O, const Vector& u) : O(O), u(u) {};
    Vector O;
    Vector u;
};

class BoundingBox {
public:
    bool intersect(const Ray& r, double& inter_distance) const {
        double tx1 = (m[0]-r.O[0]) / r.u[0];
        double tx2 = (M[0]-r.O[0]) / r.u[0];
        double txmin = std::min(tx1, tx2);
        double txmax = std::max(tx1, tx2);

        double ty1 = (m[1]-r.O[1]) / r.u[1];
        double ty2 = (M[1]-r.O[1]) / r.u[1];
        double tymin = std::min(ty1, ty2);
        double tymax = std::max(ty1, ty2);

        double tz1 = (m[2]-r.O[2]) / r.u[2];
        double tz2 = (M[2]-r.O[2]) / r.u[2];
        double tzmin = std::min(tz1, tz2);
        double tzmax = std::max(tz1, tz2);

        double actualT = std::max(txmin, std::max(tymin, tzmin));
        double furthestT = std::min(txmax, std::min(tymax, tzmax));

        if (furthestT > 0 && furthestT > actualT) {
            return true;
        }

        inter_distance = furthestT;
        return false;
    }
    Vector m, M;
};

class Object {
public:
    Object(const Vector& rho, bool is_mirror = false, bool is_transparent = false) : rho(rho), is_mirror(is_mirror), is_transparent(is_transparent) {};
    virtual bool intersect(const Ray& r, Vector& P, Vector& N,double& t) const = 0;
    Vector rho;
    bool is_mirror, is_transparent;
};

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

bool intersect_triangle(const Ray& r, const Vector& A, const Vector& B, const Vector& C, Vector& P, Vector& N, double& t) {
    Vector e1 = B - A;
    Vector e2 = C - A;
    N = cross(e1, e2);

    double inverse_uN = 1./dot(r.u, N);
    t = dot(A-r.O, N) * inverse_uN;

    if (t<0) return false;

    Vector OAcrossu = cross(A - r.O, r.u);

    double beta = dot(e2, OAcrossu) * inverse_uN;
    if (beta<0 || beta > 1) return false;

    double gamma = -dot(e1, OAcrossu) * inverse_uN;
    if (gamma < 0 || gamma > 1) return false;

    double alpha = 1 - beta - gamma;
    if (alpha < 0) return false;

    P = r.O + t*r.u;
    N.normalize();
    return true;
}

class Node{
public:
    Node(){};
    int starting_triangle;
    int ending_triangle;
    BoundingBox Bbox;
    Node* child_left;
    Node* child_right;
    bool is_leaf;
};
 
 
class TriangleMesh : public Object {
public:
    ~TriangleMesh() {}
    TriangleMesh(const Vector& rho, bool is_mirror=false, bool is_transparent=false): ::Object(rho, is_mirror, is_transparent) {
        root = new Node();
    }
    
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
 
    }

    void scale_and_translate(double scaling, const Vector& translation) {
        for (int i=0; i<vertices.size(); i++) {
            vertices[i] = vertices[i]*scaling + translation;
        }
    }

    BoundingBox compute_bounding_box(int& starting_triangle, int& ending_triangle) {
        BoundingBox bbox;
		bbox.m = Vector(1E9, 1E9, 1E9);
        bbox.M = Vector(-1E9, -1E9, -1E9);

		for (int i = starting_triangle; i < ending_triangle; i++) 
		{
			TriangleIndices triangle_idx = indices[i];
			for (int j : {triangle_idx.vtxi,triangle_idx.vtxj,triangle_idx.vtxk}) {
				Vector point = vertices[j];

				bbox.M[0] = std::max(bbox.M[0], point[0]);
				bbox.m[0]= std::min(bbox.m[0], point[0]);
				bbox.M[1] = std::max(bbox.M[1], point[1]);
				bbox.m[1] = std::min(bbox.m[1], point[1]);
				bbox.M[2] = std::max(bbox.M[2], point[2]);
				bbox.m[2] = std::min(bbox.m[2], point[2]);
			}
		}
		return bbox;
    }

    // BVH function inspired from page 49 of lecture notes
    void bvh(Node *node, int starting_triangle, int ending_triangle) {
        node->Bbox = compute_bounding_box(starting_triangle, ending_triangle);
        BoundingBox bbox = node->Bbox;
        node->starting_triangle = starting_triangle; 
        node->ending_triangle = ending_triangle;
        Vector diag = bbox.M - bbox.m;
        Vector middle_diag = bbox.m + diag*0.5; 
		int index = 0;
		double max = diag[0];
        if (diag[1] > max){ 
            max = diag[1];
            index = 1;
        }
        if (diag[2] > max){ 
            max = diag[2];
            index = 2; 
        }
        int longest_axis = index;
        int pivot_index = starting_triangle;
        for (int i=starting_triangle ; i<ending_triangle ; i++) {
			TriangleIndices T = indices[i];
			Vector barycenter = Vector(vertices[T.vtxi] + vertices[T.vtxj] + vertices[T.vtxk]) * (1./3);
            if (barycenter[longest_axis] < middle_diag[longest_axis]) { 
                std::swap(indices[i], indices[pivot_index]); 
                pivot_index++;
            }
        }
        if (pivot_index<=starting_triangle || pivot_index>=ending_triangle-1 || ending_triangle-starting_triangle<5) {
            node->child_left = nullptr;
            node->child_right = nullptr;
            node->is_leaf = true;
            return ;
        }

        node->is_leaf = false;
        node->child_left = new Node();
        node->child_right = new Node();
        bvh(node->child_left, starting_triangle, pivot_index); 
        bvh(node->child_right, pivot_index, ending_triangle);
    }

    bool intersect(const Ray& r, Vector &P, Vector &N, double &t) const {
        double inter_distance;
        BoundingBox bbox = root->Bbox;

        if (!bbox.intersect(r, inter_distance)) return false;

        t = std::numeric_limits<double>::max();
        bool has_inter = false;
        std::list<Node*> nodes_to_visit;
        nodes_to_visit.push_front(root);

        while(!nodes_to_visit.empty()) {
            Node* curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();
            // if there is one child, then it is not a leaf, so test the bounding box
            if (curNode->child_left) {
                if (curNode->child_left->Bbox.intersect(r, inter_distance)) {
                    // if (inter_distance < t) {
                    nodes_to_visit.push_back(curNode->child_left);
                    // }
                }
                if (curNode->child_right->Bbox.intersect(r, inter_distance)) {
                    // if (inter_distance < t) {
                    nodes_to_visit.push_back(curNode->child_right);
                    // }
                }
            } else {
                for (int i=curNode->starting_triangle; i<curNode->ending_triangle; i++) {

                    const Vector& A = vertices[indices[i].vtxi];
                    const Vector& B = vertices[indices[i].vtxj];
                    const Vector& C = vertices[indices[i].vtxk];

                    Vector local_P, local_N;
                    double local_t;
                    bool local_has_inter = intersect_triangle(r, A, B, C, local_P, local_N, local_t);
                    if (local_has_inter && local_t < t) {
                        t = local_t;
                        has_inter = true;
                        N = local_N;
                        P = local_P;
                    }
                }
            }
        }
        return has_inter;
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    
    BoundingBox box;
    Node *root;
};

class Sphere : public Object {
public:

    Sphere(const Vector& C, double R, const Vector& rho, bool is_mirror=false, bool is_transparent=false) : C(C), R(R), ::Object(rho, is_mirror, is_transparent) {};

    bool intersect(const Ray& r, Vector &P, Vector &N, double &t) const {

        double delta = sqr(dot(r.u, r.O - C)) - ((r.O - C).norm2() - R*R);

        if (delta >= 0) {

            double t1 = dot(r.u, C - r.O) - sqrt(delta);
            double t2 = dot(r.u, C - r.O) + sqrt(delta);
            if (t2 < 0) return false;
            if (t1 > 0) {
                t = t1;
            }
            else {
                t = t2;
            }
            P = r.O + t * r.u;
            N = P - C;
            N.normalize();
            return true;
        }
        return false;
    }

    Vector C;
    double R;
};

class Scene {
public:
    void addSphere(const Sphere* s) {objects.push_back(s);}

    void addMesh(const TriangleMesh* m) {objects.push_back(m);}

    bool intersect(const Ray& r, Vector& P, Vector& N, double& t, int& sphere_id) {
        bool has_inter = false;
        t = std::numeric_limits<double>::max();

        for (int i=0; i<objects.size(); i++) {

            Vector localP, localN;
            double localt;

            if (objects[i]->intersect(r, localP, localN, localt)) {
                if (localt<t) {
                    t = localt;
                    P = localP;
                    N = localN;
                    has_inter = true;
                    sphere_id = i;
                }
            }
        }
        return has_inter;
    }

    Vector getColor(const Ray& r, int bounce) {
        if (bounce <= 0) return Vector(0, 0, 0);

        Vector P, N;
        double t;
        Vector color(0, 0, 0);
        int sphere_id;
        if (intersect(r, P, N, t, sphere_id)) {

            if (objects[sphere_id]->is_mirror) {
                Vector mirrorDirection = r.u - 2 * dot(r.u, N) * N;

                Ray mirrorR(P + 0.001 * N, mirrorDirection);
                return getColor(mirrorR, bounce - 1);
            }
            if (objects[sphere_id]->is_transparent) {
                double n1 = 1.; // Scene's refraction index
                double n2 = 1.4; // Object's refraction index
                Vector Ntransp = N;
                if (dot(r.u, N) > 0) {
                    std::swap(n1, n2);
                    Ntransp = -Ntransp;
                }

                // Implement Fresnel Light
                double k0 = (n1-n2) * (n1-n2) / ((n1+n2) * (n1+n2));
                double u = (double) rand()/RAND_MAX;
                double K = k0 + (1 - k0) * pow(1 - std::abs(dot(N, r.u)), 5.);
                if (u < K){
                    Ray reflected_ray = Ray(P + 0.0001 * N, r.u - ( 2 * dot(r.u, N) * N ));
                    return getColor(reflected_ray, bounce - 1);
                }

                Vector tTangent, tNormal;
                tTangent = n1 / n2 * (r.u - dot(r.u, Ntransp) * Ntransp);
                double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, Ntransp)));

                // Reflection
                if (rad < 0) {
                    Vector mirrorDirection = r.u - 2 * dot(r.u, N) * N;
                    Ray mirrorR(P - 0.001*N, mirrorDirection);
                    return getColor(mirrorR, bounce - 1);
                }

                // Refraction
                tNormal = -sqrt(rad) * Ntransp;

                Ray refractedRay(P - Ntransp*0.001, tTangent + tNormal);

                return getColor(refractedRay, bounce - 1);
            }

            //else, diffuse sphere (Lambertian Model):

            // DIRECT COMPONENT
            double d2 = (L - P).norm2();
            Vector lightDir = (L - P); lightDir.normalize();
            Ray shadowRay(P + 0.001 * N, lightDir);
            Vector shadowP, shadowN;
            double shadowt;
            int shadow_id;
            bool in_shadow = false;
            if (intersect(shadowRay, shadowP, shadowN, shadowt, shadow_id)) {
                if (shadowt * shadowt < d2) {
                    in_shadow = true;
                }
            }
            if (!in_shadow)
                color = I/(4.*M_PI*d2) * objects[sphere_id]->rho/ M_PI * std::max(0., dot(N,lightDir));
            
            // ADD INDIRECT COMPONENT

            Ray indirect_ray(P + 0.0001 * N, random_cos(N));
            color += objects[sphere_id]->rho * getColor(indirect_ray, bounce - 1);

            return color;
        }
        return Vector(0, 0, 0);
    }
        

    std::vector<const Object*> objects;

    double I;
    Vector L; // light
};

int main() {
    // start timer
    auto start_time = std::chrono::high_resolution_clock::now();

    bool sphere = false; // true if show three spheres, else shows cat
    bool antialiasing = true; // sets antialiasing or not

    int W = 512;
    int H = 512;

    Scene scene;
    scene.I = 3E10;
    scene.L = Vector(-10, 20, 40);

    Sphere left_wall(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
    Sphere right_wall(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
    Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
    Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));
    Sphere front_wall(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
    Sphere behind_wall(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));

    scene.addSphere(&left_wall);
    scene.addSphere(&right_wall);
    scene.addSphere(&ceiling);
    scene.addSphere(&floor);
    scene.addSphere(&front_wall);
    scene.addSphere(&behind_wall);

    Vector camera_center(0,0,55);
    double alpha = 60. * M_PI / 180.; //field of view, rad

    int max_bounces = 5;
    int nb_paths = 64;
    std::vector<unsigned char> image(W * H * 3, 0);

    if (sphere) {
        Sphere S1(Vector(-20,0,0), 10., Vector(0., 0.5, 1.), false, true);
        Sphere S2(Vector(0,0,0), 10., Vector(1., 1., 1.), false, false);
        Sphere S3(Vector(20,0,0), 10., Vector(0.5, 0.5, 0.5), true, false);
        scene.addSphere(&S1);
        scene.addSphere(&S2);
        scene.addSphere(&S3);
    }
    else {
        TriangleMesh mesh = Vector(0.3, 0.2, 0.25);

        int init = 0;
        int max_init = (mesh.indices).size();
        mesh.readOBJ("cat.obj");
        mesh.scale_and_translate(0.6, Vector(0,-10,0));
        mesh.bvh(mesh.root, 0, (mesh.indices).size());
        scene.addMesh(&mesh);

    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        // int thread_id = omp_get_thread_num();
        for (int j = 0; j < W; j++) {
            Vector color(0,0,0);

            if (antialiasing) {
                for (int k=0; k < nb_paths; k++) {
                    // antialiasing
                    double u1 = uniform(engine);
                    double u2 = uniform(engine);
                    double randomX = sqrt(-2 * log(u1)) * cos(2. * M_PI * u2) * 0.75;
                    double randomY = sqrt(-2 * log(u1)) * sin(2. * M_PI * u2) * 0.75;

                    Vector ray_dir;
                    ray_dir[0] = j - W / 2. + 0.5 + randomX;
                    ray_dir[1] = -i + H / 2. + 0.5 + randomY;
                    ray_dir[2] = - W / (2. * tan(alpha / 2.));
                    ray_dir.normalize();
                    Ray r(camera_center, ray_dir);

                    color += scene.getColor(r, max_bounces);
                }
            }

            else {
                Vector ray_dir;
                ray_dir[0] = j - W / 2. + 0.5;
                ray_dir[1] = -i + H / 2. + 0.5;
                ray_dir[2] = - W / (2. * tan(alpha / 2.));
                ray_dir.normalize();
                Ray r(camera_center, ray_dir);

                for (int k=0; k < nb_paths; k++) {
                    color += scene.getColor(r, max_bounces);
                }
            }

            color = color / nb_paths;

            image[(i * W + j) * 3 + 0] = std::min(255.,std::pow(color[0], 0.45));
            image[(i * W + j) * 3 + 1] = std::min(255.,std::pow(color[1], 0.45));
            image[(i * W + j) * 3 + 2] = std::min(255.,std::pow(color[2], 0.45));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    
    // end timer and print the time it took to run the code
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}