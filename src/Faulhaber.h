#ifndef FAULHABER_H
#define FAULHABER_H

#include <stdexcept> // For std::out_of_range (used in the function's exception)
#include <cmath>
/**
 * @brief Calculates the sum of powers (Faulhaber's formula) for degrees 1-4.
 * * @param n The upper limit of the sum.
 * @param deg The power (degree) to sum (must be 1, 2, 3, or 4).
 * @return The sum as a double.
 * @throws std::out_of_range If the degree is not between 1 and 4.
 */
double Faulhaber(int n, int deg);

#endif // FAULHABER_H
