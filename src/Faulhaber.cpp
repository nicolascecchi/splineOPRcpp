#include "Faulhaber.h"
#include <cmath>
#include <stdexcept>

double Faulhaber(int n, int deg){
    switch (deg){
        case 1:
            return 0.5* (std::pow(n,2) + n);
        case 2:
            return 1.0/3.0 * (std::pow(n,3)+ 3./2. * std::pow(n,2) + 0.5*n);
        case 3:
            return 1./4. * (std::pow(n,4) + 2. * std::pow(n,3) + std::pow(n,2));
        case 4:
            return 1./5. * (std::pow(n,5) + 5./2. * std::pow(n,4) + 10./6. * std::pow(n,3) - 5./30. * n);
        default:
            throw std::out_of_range("Degree must be between 1 and 4");
    }
}
