#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <cmath>

double sqr(double x){
    return std::pow(x, 2);
    
}

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
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
// New
Vector operator-(const Vector& a){
    return Vector(-a[0], -a[1], -a[2]);
}

class Ray {
public:
	Ray(const Vector &origin, const Vector &direction)
    {
        U = direction;
        O = origin;
    }
    Vector O;
    Vector U;
};

class Sphere {
public:
    Sphere(const Vector& center, double radius, const Vector& rho, bool mirror = false, bool transparency = false, bool hollow = false)
    {
        C = center; R = radius; RHO = rho; is_mirror = mirror; is_transparent = transparency; is_hollow = hollow;
    }
	bool intersect(const Ray& r, Vector& P, Vector &N, double& t)
    // Ray-sphere intersection
    {
        double delta = sqr(dot(r.U, r.O - C)) - ((r.O - C).norm2() - R*R);
        
        if (delta >= 0)
        {
            double t1 = dot(r.U, C - r.O) - sqrt(delta);
            double t2 = dot(r.U, C - r.O) + sqrt(delta);
            if (t2 < 0)
                return false;
            
            if (t1 > 0)
            {
                t = t1;
            }
            else
            {
                t = t2;
            }
            P = r.O + (t *  r.U);
            N = P - C;
            N.normalize();
            return true;
        }
        return (delta >= 0);
    }
    
    Vector C;
    double R;
    Vector RHO;
    bool is_mirror;
    bool is_transparent;
    bool is_hollow;

};

class Scene{
public:
    void addSphere(const Sphere& s)
    {
        objects.push_back(s);
    }
	
	bool intersect(const Ray& r, Vector& P, Vector& N, double& t, int& sphere_id)
    {
        Vector localP, localN;
        double localt = 0.0;
        bool has_inter = false;
        
        t = std::numeric_limits<double>::max();
        
        for(int i = 0; i < objects.size(); i++)
        {
            if(objects[i].intersect(r, localP, localN, localt))
            {
                if(localt < t)
                {
                    t = localt;
                    P = localP;
                    N = localN;
                    sphere_id = i;
                    has_inter = true;
                }
            }
        }
        return has_inter;
    }

	Vector getColor(const Ray& r, int bounces)
    {
        if(bounces <=0)
            return Vector(10, 0, 0);
        
        Vector P, N;
        double t;
        int sphere_id;
        Vector color(255,255,255);
        
        if(intersect(r, P, N, t, sphere_id))
        {
            if(objects[sphere_id].is_hollow)
                //Hollowness
            {
                N = -N;
            }
            
            if(objects[sphere_id].is_mirror)
                //Mirors
            {
                Vector mirrorDirection = r.U - 2*dot(r.U, N)*N;
                Ray mirrorR(P+0.001*N, mirrorDirection);
                return getColor(mirrorR, bounces - 1);
            }

            if (objects[sphere_id].is_transparent)
                //Transparence
            {
                double n1 = 1;
                double n2 = 1.4;
                Vector NTransparency = N;
                
                if(dot(r.U, N) > 0)
                {
                    std::swap(n1, n2);
                    NTransparency = -NTransparency;
                }
                
                Vector tTangent, tNormal;
                tTangent = (n1/n2) * (r.U - (dot(r.U, NTransparency)*NTransparency));
                double rad = 1 - (sqr(n1/n2) * (1 - sqr(dot(r.U, NTransparency))));
                
                if(rad < 0)
                {
                    Vector mirrorDirection = r.U - 2*dot(r.U, N)*N;
                    Ray mirrorR(P - 0.001*N, mirrorDirection);
                    return getColor(mirrorR, bounces - 1);
                }
                else
                {
                    tNormal = -sqrt(rad) * NTransparency;
                    Ray refractedRay(P - NTransparency*0.001, tTangent + tNormal);
                    return getColor(refractedRay, bounces - 1);
                }
			}
			Vector shadowP, shadowN;
            double shadowt;
            int shadow_id;
            bool in_shadow = false;

            double d2 = (L - P).norm2();
            Vector lightDir = (L - P); lightDir.normalize();
            Ray shadowRay(P+0.001*N, lightDir);
            
            if (intersect(shadowRay, shadowP, shadowN, shadowt, shadow_id))
                //Handling shadows
            {
                if ((shadowt * shadowt) < d2)
                {
                    in_shadow = true;
                }
            }
            if (!in_shadow)
            {
                //Setting pixel color
                color = (I/(4.* M_PI * d2)) * (objects[sphere_id].RHO/M_PI) * std::max(0., dot(N, lightDir));
            }
            return color;
        }
        return color;
    }
	
	std::vector<Sphere> objects;
    double I;
    Vector L;

;


int main() {
    int W = 512;
    int H = 512;

    Scene s;
    s.I = 2E10;
    s.L = Vector(-10, 20, 40);

    Vector camera_centre(0, 0, 55);
    double alpha = (60 * M_PI) / 180; //F.O.V (in rad.)
    int max_bounces = 6;

// Settings for white disk and black background only
//    Sphere S(Vector(0,0,0), 12, Vector(255, 255, 255));
//    scene.addSphere(S);
    

// Scene setup
    Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));
    Sphere ceiling(Vector(0, 1000,0), 940, Vector(0.3, 0.5, 0.3));
    Sphere front_wall(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
    Sphere left_wall(Vector(-1000, 0, 0), 940 , Vector(0.5, 0.8, 0.1));
    Sphere right_wall(Vector(1000, 0, 0), 940,  Vector(0.9, 0.2 , 0.3));
    Sphere behind_wall(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
    
//Uncomment this line for a simple white sphere
//    Sphere S(Vector(0, 0, 0), 10, Vector(0.5, 0.5, 0.5));

    
//Uncomment this line for a mirror sphere
//    Sphere S(Vector(0, 0, 0), 10, Vector(0.5, 0.5, 0.5), true);
    
//Uncomment this line for a transparent sphere
    Sphere S(Vector(0, 0, 0), 12, Vector(0.5, 0.5, 0.5), false, true);
    
//Uncomment these lines for a hollow sphere
//        Sphere S1(Vector(0, 0, 0), 12, Vector(0.5, 0.5, 0.5), false, true);
//        Sphere S2(Vector(0, 0, 0), 11, Vector(0.5, 0.5, 0.5), false, true, true);
    

// Adding elements to the scene (order is important !)
    s.addSphere(S);
    
//Uncomment these lines for a hollow sphere and comment the above scene.addSphere(S);
//    scene.addSphere(S1);
//    scene.addSphere(S2);
    
    
    s.addSphere(floor);
    s.addSphere(ceiling);
    s.addSphere(left_wall);
    s.addSphere(right_wall);
    s.addSphere(front_wall);
    s.addSphere(behind_wall);

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector ray_dir;
            
            //Launcing rays
            ray_dir[0] = j - W / 2. +0.5;
            ray_dir[1] = -i + H / 2. +0.5;
            ray_dir[2] = - W / (2. * tan(alpha/2));
            ray_dir.normalize();
            
            Ray r(camera_centre, ray_dir);
            Vector color = s.getColor(r, max_bounces);
            
            //Gamma correction
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}
    