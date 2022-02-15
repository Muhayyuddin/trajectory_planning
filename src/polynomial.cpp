/*************************************************************************\
   Copyright 2022: Muhayy Ud Din.
    All Rights Reserved.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 \*************************************************************************/
/* Author: Muhayy Ud Din */
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "polynomial.h"
#define _USE_MATH_DEFINES // for sin/log
#define DEBUG1
//#define DEBUG2

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief It compute the n uniform numbers between start and end number.
 *
 * @tparam T
 * @param start_in The starting number
 * @param end_in The ending number.
 * @param n uniform intervals between start and end.
 * @return std::vector<double>
 */
std::vector<double> Polynomial::TimeInterval(double start_in, double end_in, int n)
{

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(n);

    if (num == 0)
    {
        return linspaced;
    }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end);
    return linspaced;
}
/**
 * @brief This function will generat the system of linear equation for Quantic polynomial (of order 5)
 * i.e.  s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5.
 * Compute the coefficient using A^-1 * B = C.
 *
 * @param init_values Initial values of [s, s_dot, s_double_dot], where s is position, s_dot, is velocity and s_double_dot is acceleration.
 * @param final_values Final values of [s, s_dot, s_double_dot], where s is position, s_dot, is velocity and s_double_dot is acceleration.
 * @param ti initial time
 * @param tf final time.
 * @return vector<double> vector of length 6 contains computed coefficients of the polynomial.
 */

std::vector<double> Polynomial::GenerateQuanticPolynomial(std::vector<double> init_values, std::vector<double> final_values, double ti, double tf)
{
#ifdef DEBUG2
    std::cout << "Generating Matrix A: " << std::endl;
#endif
    // matrix representing equations.
    MatrixXd A = MatrixXd(6, 6);
    A << 1,   ti,  ti*ti,   ti*ti*ti, ti*ti*ti*ti,  ti*ti*ti*ti*ti,
	     0,   1,   2*ti,    3*ti*ti,  4*ti*ti*ti,   5*ti*ti*ti*ti,
	     0,   0,   2,       6*ti,     12*ti*ti,     20*ti*ti*ti,
         1,   tf,  tf*tf,   tf*tf*tf, tf*tf*tf*tf,  tf*tf*tf*tf*tf,
	     0,   1,   2*tf,    3*tf*tf,  4*tf*tf*tf,   5*tf*tf*tf*tf,
	     0,   0,   2,       6*tf,     12*tf*tf,     20*tf*tf*tf;

// A<< 1, 0, 0, 0, 0, 0,
//     0, 1, 0, 0, 0, 0,
//     0, 0, 1, 0, 0, 0,
//     0, 0, 0, 1, 0, 0,
//     0, 0, 0, 0, 1, 0,
//     0, 0, 0, 0, 0, 1;
#ifdef DEBUG1
    std::cout << "The matrix A is: " << std::endl;
    std::cout << A << std::endl;
#endif
#ifdef DEBUG2
    std::cout << "Generating Matrix B: " << std::endl;
#endif
    MatrixXd B = MatrixXd(6, 1);
    B << init_values[0], // initial position
        init_values[1],  // initial velocity
        init_values[2],  // initial acceleration
        final_values[0], // final position
        final_values[1], // final velocity
        final_values[2]; // final acceleration

#ifdef DEBUG1
    std::cout << "The matrix B is: " << std::endl;
    std::cout << B << std::endl;
#endif
#ifdef DEBUG2
    std::cout << "Computing Inverse of A: " << std::endl;
#endif

    MatrixXd Ai = A.inverse();

#ifdef DEBUG1
    std::cout << "The Inverse of matrix A is: " << std::endl;
    std::cout << Ai << std::endl;
#endif
#ifdef DEBUG2
    std::cout << "Computing Coefficients using C = A^-1 * B : " << std::endl;
#endif

    MatrixXd C = Ai * B;

#ifdef DEBUG1
    std::cout << "The computed Coefficients are:  " << std::endl;
    std::cout << C << std::endl;
#endif

    std::vector<double> result;
    for (int i = 0; i < C.size(); i++)
    {
        result.push_back(C.data()[i]);
    }

    return result;
}
/**
 * @brief This function will compute the points from the polynomials
 *
 * @param coefficients coefficients of the polynomials.
 * @param from_t starting point
 * @param to_t ending point
 * @param n number of points to generate between start and end.
 * @return std::vector<double>  vector of points from the polynomial
 */
std::vector<double> Polynomial::PolynomialPointsFromCoefficients(std::vector<double> coefficients, double from_t, double to_t, int n)
{

#ifdef DEBUG2
    std::cout << "Computing [ " << n << " ] Uniform points between [ " << from_t << " , " << to_t << " ]" << std::endl;
#endif
    std::vector<double> ts = TimeInterval(from_t, to_t, n);
    std::vector<double> polynomial;
    polynomial.resize(n);
    for (auto &p : polynomial)
    {
        p = 0;
    }
#ifdef DEBUG1
    std::cout << " Point Intervals:  " << polynomial.size() << std::endl;
    std::cout << " [ ";
    for (const auto &t : ts)
    {
        std::cout << t << " , ";
    }
    std::cout << " ] " << std::endl;
#endif
#ifdef DEBUG2
    std::cout << "Computing polynomials: " << std::endl;
#endif
    int power{0};
    for (size_t t{0}; t < ts.size(); ++t)
    {
        for (const auto &coefficient : coefficients)
        {
            polynomial[t] = polynomial[t] + (coefficient * std::pow(ts[t], power));
            ++power;
        }
        power = 0;
    }
#ifdef DEBUG2
    std::cout << "Computed polynomials." << std::endl;
#endif
#ifdef DEBUG1
    std::cout << " Polynomial vector size is :  " << polynomial.size() << std::endl;
    std::cout << " [ ";
    for (const auto &p : polynomial)
    {
        std::cout << p << " , ";
    }
    std::cout << " ] " << std::endl
              << std::endl;
#endif

    return polynomial;
}

/**
* @brief This function will compute the velocities at time t using vel = C_1 + 2C_2t + 3C_3t^2 + 4C_4t^3 + 5C_5t^4.
*
* @param coefficients The computed coefficients
* @param time the vector of uniform time intervals from start time to end time.
* @return std::vector<double> vector of computed velocities at given time.
    */
std::vector<double> Polynomial::ComputeVelocity(std::vector<double> coefficient, std::vector<double> time)
{
    std::vector<double> velocity;
    for (size_t t{0}; t < time.size(); ++t)
    {
        velocity.push_back(coefficient[1] + 2 * coefficient[2] * time[t] + 
                                            3 * coefficient[3] * std::pow(time[t], 2.0) + 
                                            4 * coefficient[4] * std::pow(time[t], 3.0)+ 
                                            5 * coefficient[5] * std::pow(time[t], 4.0));
    }
    return velocity;
}

/**
 * @brief This function will compute the acceleration at time t using acc = 2C_2 + 6C_3t^2 + 12 C_5t^3.
 *
 * @param coefficients The computed coefficients
 * @param time the vector of uniform time intervals from start time to end time.
 * @return std::vector<double> derivative results
 */
std::vector<double> Polynomial::ComputeAcceleration(std::vector<double> coefficient, std::vector<double> time)
{
    std::vector<double> acceleration;
    for (size_t t{0}; t < time.size(); ++t)
    {
        acceleration.push_back(2 * coefficient[2] + 6 * coefficient[3] * time[t] + 
                               12 * coefficient[4] * std::pow(time[t], 2.0) + 
                               20 * coefficient[5] * std::pow(time[t], 3.0));
    }
    return acceleration;
}
