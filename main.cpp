#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <math.h>

#define M_PI  3.14159265358979323846

#include <iostream>
#include <chrono>
#include <random>

static std::default_random_engine engine(10) ; // random s e e d = 10
static std::uniform_real_distribution<double> uniform ( 0 , 1 );

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
static inline double sqr(double x) {
	return x*x;
}

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& b) {
	return Vector( -b[0], -b[1], -b[2]);
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
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
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
	Vector O,dirrection;
	Ray(const Vector& pos, const Vector& dir):O(pos), dirrection(dir) {}
};

class Sphere {
public:
	Vector center;
	double r;
	Vector albedo;
	bool isMirror,isTrans;

	Sphere(const Vector& center, double r, const Vector& color,bool ism = false, bool ist = false): 
	center(center), r(r), albedo(color), isMirror(ism), isTrans(ist) {};

	bool intersection(const Ray& ray, Vector &P, Vector &N, double &t){
		double dotval = dot(ray.dirrection,ray.O - center);
		double delta = dotval*dotval - ((ray.O - center).norm2() - r*r);
		if (delta<0){
			return false;
		}

		double sqrtdelta = sqrt(delta);
		double t1 = -dotval - sqrtdelta;
		double t2 = -dotval + sqrtdelta;

		if (t2<0){
			return false;
		}

		if (t1>0){
			t = t1;
			P = ray.O + t1 * ray.dirrection;
		}
		else{
			t = t2;
			P = ray.O + t2 * ray.dirrection;
		}

		N = P - center;
		N.normalize();
		return true;
	}

};


class Scene {
public:
	Scene() {};
	void addObject(const Sphere& s){
		objects.push_back(s);
	}
	
	std::vector<Sphere> objects;

	bool intersection(const Ray& ray, Vector &P, Vector &N, double &t, int &ID){
		t = 2E10;
		bool f_inter = false;
		for (int i =0; i < objects.size(); i++){
			Vector localP,localN;
			double localt;

			bool inter = objects[i].intersection(ray,localP,localN,localt);
			if (inter){
				f_inter = true;
				if (localt < t){
					ID = i;
					t = localt;
					P = localP;
					N = localN;
				}
			}
		}

		return f_inter;
	}
//------
// Get color function 
Vector getColor(const Ray& ray, int bounce) {
    if (bounce <= 0) return Vector(0, 0, 0);  // Prevent infinite recursion

    Vector P, N;
    double t;
    int ID;
    Vector color(0, 0, 0);

    double I = 1E10;
    Vector L(-10, 20, 40); // Light position

    if (intersection(ray, P, N, t, ID)) {
        Sphere& obj = objects[ID];

        if (obj.isMirror) {
            Vector reflectDir = ray.dirrection - 2 * dot(ray.dirrection, N) * N;
            Ray reflectRay(P + 0.001 * N, reflectDir);
            return getColor(reflectRay, bounce - 1);
        }
        if (obj.isTrans) {
            double n1 = 1.; // Air's refractive index
            double n2 = 1.4; // Sphere's refractive index
            Vector N_Trans = N;
            if (dot(ray.dirrection, N) > 0) {
                std::swap(n1, n2);
                N_Trans = -N_Trans;
            }

            double cosI = -dot(N_Trans, ray.dirrection);

			// Fresnel Implementation start
            // double r0 = (n1 - n2) / (n1 + n2);
            // r0 = r0 * r0;
            // double R = r0 + (1 - r0) * pow((1 - cosI), 5);

            // if (uniform(engine) < R) {
            //     Ray reflectedRay = Ray(P + 0.0001 * N, ray.dirrection - 2 * dot(ray.dirrection, N) * N);
            //     return getColor(reflectedRay, bounce - 1);
            // }
			// Fresnel implementaiton end.

            double radic = 1 - (n1 / n2) * (n1 / n2) * (1 - cosI * cosI);
            if (radic < 0) { // Total internal reflection
                Vector reflectDir = ray.dirrection - 2 * dot(ray.dirrection, N) * N;
                Ray reflectRay(P + 0.001 * N, reflectDir);
                return getColor(reflectRay, bounce - 1);
            }

            Vector Tt = (n1 / n2) * (ray.dirrection - dot(ray.dirrection, N_Trans) * N_Trans);
            Vector Tn = -sqrt(radic) * N_Trans;
            Vector T = Tt + Tn;
            Ray refractRay(P - 0.001 * N_Trans, T);
            return getColor(refractRay, bounce - 1);
        }

        Vector lightDir = L - P;
        double lightDist2 = lightDir.norm2();
        lightDir.normalize();
        Ray shadowRay(P + 0.001 * N, lightDir);
        Vector shadowP, shadowN;
        double shadowt;
        int shadowID;
        bool in_shadow = intersection(shadowRay, shadowP, shadowN, shadowt, shadowID) && shadowt * shadowt < lightDist2;
        if (!in_shadow) {
            color += I / (4 * M_PI * lightDist2) * obj.albedo / M_PI * std::max(0., dot(N, lightDir));
        }

        // Optionally add indirect lighting
        // Ray indirectRay(P + 0.0001 * N, random_cos(N));
        // color += obj.albedo * getColor(indirectRay, bounce - 1);

        return color;
    }

    return Vector(0, 0, 0); // No intersection found, return background color
}

};



// ------------------------------------------------------------------------------------
// MAIN

void boxMuller ( double stdev , double &x , double &y ) {
	double r1 = uniform ( engine ) ;
	double r2 = uniform ( engine ) ;
	x = sqrt(-2 * log ( r1 ) ) *cos ( 2 * M_PI*r2 ) *stdev ;
	y = sqrt(-2 * log ( r1 ) ) *sin ( 2 * M_PI*r2 ) *stdev ;
};

int main() {
	// Computing time 
    auto start_time = std::chrono::high_resolution_clock::now();

	int W = 512;
	int H = 512;

	double fov = 60 * M_PI / 180;
	Vector Q(0,0,55);

	Sphere S1       (Vector(20,0,0),    10, Vector(0.9,0.9,0.9));
	Sphere S2       (Vector(0,0,0),     10, Vector(0.7,0.3,0.1),true);
	Sphere S3       (Vector(-20,0,0),     10, Vector(0.7,0.3,0.1),false,true);

	Sphere Sfloor   (Vector(0,-1000,0), 990, Vector(0.2,0.9,0.1));
	Sphere Sceilling(Vector(0,1000,0),  940, Vector(0.4,0.4,0.2));
	Sphere Sleft    (Vector(-1000,0,0), 940, Vector(0.1,0.1,0.8));
	Sphere Sright   (Vector(1000,0,0),  940, Vector(0.7,0.2,0.5));
	Sphere Sback    (Vector(0,0,1000),  940, Vector(0.1,0.4,0.2));
	Sphere Sfront   (Vector(0,0,-1000), 940, Vector(0.9,0.3,0.5));

	Scene scene;
	scene.addObject(S1);
	scene.addObject(S2);
	scene.addObject(S3);

	scene.addObject(Sfloor);
	scene.addObject(Sceilling);
	scene.addObject(Sleft);
	scene.addObject(Sright);
	scene.addObject(Sback);
	scene.addObject(Sfront);


	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector u(j - W/2 + 0.5, H/2 - i - 0.5, -W / (2 * tan( fov / 2)));
			u.normalize();
			Ray ray(Q,u);
			Vector color = scene.getColor(ray,5);
			
			
			image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
			image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	// end timer and print the time it took to run the code

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

	return 0;
}

