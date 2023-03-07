#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <string>

struct vector3d {
    double x,y,z;

    double& operator()(int i) {
      if (i == 0) return x;
      else if (i == 1) return y;
      else if (i == 2) return z;
      else std::cout << "Slozka neexistuje\n";
    }
};

vector3d operator*(const vector3d& v, double a) {
    vector3d r;
    r.x = v.x * a; r.y = v.y * a; r.z = v.z * a;
    return r;
}

vector3d operator+(const vector3d& v, const vector3d& u) {
    vector3d r;
    r.x = v.x + u.x; r.y = v.y + u.y; r.z = v.z + u.z;
    return r;
}

vector3d operator*(double a, const vector3d& v) {
    return v * a;
}

vector3d operator-(const vector3d& v, const vector3d& u) {
    vector3d r;
    r.x = v.x - u.x; r.y = v.y - u.y; r.z = v.z - u.z;
    return r;
}

std::ostream& operator<<(std::ostream& out, const vector3d& v) {
    return (out << v.x << ", " << v.y << ", " << v.z);
}

void vypis(double cislo) {
    std::cout << cislo << "\n";
}

void vypis(std::string s) {
    std::cout << s << "\n";
}

void vypis(vector3d v) {
    std::cout << v << "\n";
}
  
#endif
