#include "AnisotropicDiffusion.h"

#include <cmath>

AnisotropicDiffusion::AnisotropicDiffusion(const Rcpp::IntegerMatrix& _image, unsigned int _iterations, double _lambda, double _k)
    : image(_image), iterations(_iterations), lambda(_lambda), k(_k)
{
}

void AnisotropicDiffusion::setTemperature(unsigned int x, unsigned int y, unsigned char temperature)
{
    image(y,x) = temperature;
}

unsigned char AnisotropicDiffusion::getTemperature(unsigned int x, unsigned int y) {
    // Only take red channel (r = g = b)
    return image(y, x);
}

void AnisotropicDiffusion::applyDiffusion()
{
    if(iterations == 0)
      return;

    // Diffusion algorithm
    for(int i = 0; i < image.ncol(); i ++)
    {
        for(int j = 0; j < image.nrow(); j ++)
        {
            if (i == 0 || i == image.ncol() - 1 || j == 0 || j == image.nrow() - 1)
            {
                setTemperature(i, j, getTemperature(i, j));
                continue;
            }

            double NI = getTemperature(i, j - 1) - getTemperature(i, j);
            double SI = getTemperature(i, j + 1) - getTemperature(i, j);
            double EI = getTemperature(i + 1, j) - getTemperature(i, j);
            double WI = getTemperature(i - 1, j) - getTemperature(i, j);

            double cN = exp( -pow(std::abs(NI) / k, 2) );
            double cS = exp( -pow(std::abs(SI) / k, 2) );
            double cE = exp( -pow(std::abs(EI) / k, 2) );
            double cW = exp( -pow(std::abs(WI) / k, 2) );

            double newTemperature = getTemperature(i, j) + lambda * (cN * NI + cS * SI + cE * EI + cW * WI);
            setTemperature(i, j, newTemperature);
        }
    }

    // Repeat 'iterations' times
    iterations--;
    applyDiffusion();
}
