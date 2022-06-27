#ifndef ANISOTROPICDIFFUSION_H
#define ANISOTROPICDIFFUSION_H

#include <Rcpp.h>

class AnisotropicDiffusion {
private:
    Rcpp::IntegerMatrix image;
    unsigned int iterations;
    double lambda, k;
public:
    AnisotropicDiffusion(const Rcpp::IntegerMatrix& _image, unsigned int _iterations, double _lambda, double _k);
    AnisotropicDiffusion() = default;
    ~AnisotropicDiffusion() = default;
private:
    void setTemperature(unsigned int x, unsigned int y, unsigned char temperature);
    unsigned char getTemperature(unsigned int x, unsigned int y);
public:
    void applyDiffusion();
public:
    inline Rcpp::IntegerMatrix& getImage() { return image; }
    inline unsigned int getIterations() const { return iterations; }
    inline double getLambda() const { return lambda; }
    inline double getK() const { return k; }
};

#endif
