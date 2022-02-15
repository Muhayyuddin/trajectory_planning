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
#define _USE_MATH_DEFINES // for sin/log
//#define DEBUG1
//#define DEBUG2

using Eigen::MatrixXd;
using Eigen::VectorXd;
class Polynomial
{

public:
   /**
    * @brief It compute the n uniform numbers between start and end number.
    *
    * @tparam T
    * @param start_in The starting number
    * @param end_in The ending number.
    * @param n uniform intervals between start and end.
    * @return std::vector<double>
    */
   std::vector<double> TimeInterval(double start_in, double end_in, int n);

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

   std::vector<double> GenerateQuanticPolynomial(std::vector<double> init_values, std::vector<double> final_values, double ti, double tf);

   /**
    * @brief This function will compute the points from the polynomials
    *
    * @param coefficients coefficients of the polynomials.
    * @param from_t starting point
    * @param to_t ending point
    * @param n number of points to generate between start and end.
    * @return std::vector<double>  vector of points from the polynomial
    */
   std::vector<double> PolynomialPointsFromCoefficients(std::vector<double> coefficients, double from_t, double to_t, int n);

   /**
    * @brief This function will compute the velocities at time t using vel = C_1 + 2C_2t + 3C_3t^2 + 4C_4t^3 + 5C_5t^4.
    *
    * @param coefficients The computed coefficients
    * @param time the vector of uniform time intervals from start time to end time.
    * @return std::vector<double> vector of computed velocities at given time.
    */
   std::vector<double> ComputeVelocity(std::vector<double> coefficient, std::vector<double> time);
   /**
    * @brief This function will compute the acceleration at time t using acc = 2C_2 + 6C_3t^2 + 12 C_5t^3.
    *
    * @param coefficients The computed coefficients
    * @param time the vector of uniform time intervals from start time to end time.
    * @return std::vector<double> vector of computed acceleration at given time.
    */
   std::vector<double> ComputeAcceleration(std::vector<double> coefficient, std::vector<double> time);
};