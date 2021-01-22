//[header]
// This program illustrates how the concept of vector and matrix can be implemented
// in C++. This is a light version of the implementation. It contains the most
// essential methods to manipulate vectors and matrices. It should be enough
// for most projects. Vectors and matrices are really the alphabet as we said
// in the lesson of any graphics application. It's really important you feel
// confortable with these techniques especially with the concepts of
// normalizing vectors, computing their length, computing the dot and cross products
// of two vectors, and the point- and vector-matrix multiplication (and knowing
// the difference between the two).
//[/header]
//[compile]
// c++ geometry.cpp  -o geometry -std=c++11
//[/compile]
//[ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//[/ignore]
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>

template <typename T>
class Vec3
{
public:
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(const T &xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    T x, y, z;

    T length()
    {
        return sqrt(x * x + y * y + z * z);
    }

    T dot(const Vec3<T> &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    Vec3<T> cross(const Vec3<T> &v) const
    {
        return Vec3<T>(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x);
    }

    Vec3<T> &normalize()
    {
        T len = dot(*this);
        if (len > 0)
        {
            T invLen = 1 / sqrt(len);
            x *= invLen, y *= invLen, z *= invLen;
        }
        return *this;
    }

    Vec3<T> operator+(const Vec3<T> &v) const
    {
        return Vec3<T>(x + v.x, y + v.y, z + v.z);
    }

    Vec3<T> operator-(const Vec3<T> &v)
    {
        return Vec3<T>(x - v.x, y - v.y, z - v.z);
    }

    Vec3<T> operator+(const T &r) const
    {
        return Vec3<T>(x * r, y * r, z * r);
    }
};

template <typename T>
class Matrix44
{
public:
    Matrix44() {}
    const T *operator[](uint8_t i) const { return m[i]; }
    T *operator[](uint8_t i) { return m[i]; }
    T m[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

    Matrix44 operator*(const Matrix44 &rhs) const
    {
        Matrix44 mult;
        for (uint8_t i = 0; i < 4; ++i)
        {
            for (uint8_t j = 0; j < 4; ++j)
            {
                mult[i][j] = m[i][0] * rhs[0][j] +
                             m[i][1] * rhs[1][j] +
                             m[i][2] * rhs[2][j] +
                             m[i][3] * rhs[3][j];
            }
        }
        return mult;
    }

    Vec3<T> multVecMatrix(const Vec3<T> &v)
    {
#ifdef ROWMAJOR
        return Vec3<T>(
            v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
            v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
            v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2]);
#else
        return Vec3<T>(
            v.x * m[0][0] + v.y * m[0][1] + v.z * m[0][2],
            v.x * m[1][0] + v.y * m[1][1] + v.z * m[1][2],
            v.x * m[2][0] + v.y * m[2][1] + v.z * m[2][2]);
#endif
    }

    void multVecMatrix(const Vec3<T> &src, Vec3<T> &dst) const
    {
        T a, b, c, w;
        a = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0] + m[3][0];
        b = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1] + m[3][1];
        c = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2] + m[3][2];
        w = src.x * m[0][3] + src.y * m[1][3] + src.z * m[2][3] + m[3][3];

        if (w != 1 && w != 0)
        {
            dst.x = a / w;
            dst.y = b / w;
            dst.z = c / w;
        }
    }

    void multDirMatrix(const Vec3<T> &src, Vec3<T> &dst) const
    {
        dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0];
        dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1];
        dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2];
    }

    Vec3 vecMatMult(
        Vec3 v,
        float AXx, float AXy, float AXz, float Tx,
        float AYx, float AYy, float AYz, float Ty,
        float AZz, float AZy, float AZz, float Tz)
    {
    return Vec3( 
        v.x * AXx + v.y * AYx + v.z * AZx + Tx, 
        v.x * AXy + v.y * AYy + v.z * AZy + Ty, 
        v.x * AXz + v.y * AZz + v.z * AZz + Tz
    }

    Matrix44 transpose() const
    {
        Matrix44 transpMat;
        for (uint8_t i = 0; i < 4; ++i)
        {
            for (uint8_t j = 0; j < 4; ++j)
            {
                transpMat[i][j] = m[j][i];
            }
        }

        return transpMat;
    }
};

typedef Matrix44<float> Matrix44f;
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#define MAX_ITER 10e8

int main(int argc, char **argv)
{
    clock_t start = clock();
    Vec3<float> v(1, 2, 3);
    Matrix44<float> M;
    float *tmp = &M.m[0][0];
    for (int i = 0; i < 16; i++)
        *(tmp + i) = drand48();
    for (int i = 0; i < MAX_ITER; ++i)
    {
        Vec3<float> vt = M.multVecMatrix(v);
    }
    fprintf(stderr, "Clock time %f\n", (clock() - start) / float(CLOCKS_PER_SEC));
    return 0;
}
