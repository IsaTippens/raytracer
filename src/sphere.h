#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"

class Sphere : public Hittable {
    public:
        Sphere() {}
        Sphere(Point3 cen, double r) : center(cen), radius(r) {};

        virtual bool hit(
            const Ray& r, double t_min, double t_max, hit_record& rec 
        ) const override;

    public:
        Point3 center;
        double radius
};



#endif