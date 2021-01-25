#ifndef TEXTURE_H
#define TEXTURE_H

#include "utilities.h"

class Texture {
    public:
        virtual Color value(double u, double v, const Point3& p) const = 0;
};

class Solid_Color : public Texture {
    public:
        Solid_Color() {}
        Solid_Color(Color c) : color_value(c) {}

        Solid_Color(double red, double green, double blue)
            :   Solid_Color(Color(red, green, blue)) {}
        
        virtual Color value(double u, double v, const Point3& p) const override {
            return color_value;
        }

    private:
        Color color_value;
};

#endif