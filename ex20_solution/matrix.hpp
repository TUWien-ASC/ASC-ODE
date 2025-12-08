#pragma once
#include <vector>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

struct Matrix {
    int r=0, c=0;
    std::vector<double> d;
    Matrix() = default;
    Matrix(int rr, int cc) { resize(rr,cc); }
    void resize(int rr, int cc) { r=rr; c=cc; d.assign(r*c, 0.0); }
    void zero() { std::fill(d.begin(), d.end(), 0.0); }
    double& operator()(int i,int j){ assert(i>=0 && i<r && j>=0 && j<c); return d[i*c + j]; }
    double operator()(int i,int j) const { assert(i>=0 && i<r && j>=0 && j<c); return d[i*c + j]; }
    int rows() const { return r; }
    int cols() const { return c; }
    void print(const char* name=nullptr) const {
        if(name) std::cout<<name<<":\n";
        for(int i=0;i<r;++i){
            for(int j=0;j<c;++j) std::cout<<std::setw(12)<<operator()(i,j)<<" ";
            std::cout<<"\n";
        }
    }
};
