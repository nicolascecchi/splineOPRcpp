#include "Faulhaber.h"
#include <cmath>
#include <stdexcept>

double Faulhaber(int n, int deg)
{
    switch (deg)
    {
        case 1:
            return 0.5 * n * (n + 1);
        case 2:
            return 1./6. * n * (n + 1) * (2*n + 1);
        case 3:
            return 1./4. * n*n*(n + 1)*(n + 1);
        case 4:
            return 1./5. * n * (n + 1) * (n + 0.5) * (n - (-3+sqrt(21.))/6.) * (n - (-3-sqrt(21.))/6.);     
        default:
            throw std::out_of_range("Degree must be between 1 and 4");
    }
}
