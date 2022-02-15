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
#include "polynomial.h"
#include "../visulization/matplotlib-cpp/matplotlibcpp.h"

#define DEBUG1
#define DEBUG2
namespace plt = matplotlibcpp;

struct JointTrajectoryConstraints
{

    std::vector<double> initial_joint_constraints; // initial position, velocity and acceleration
    std::vector<double> final_joint_constraints;   // final position, velocity and acceleration
    double ti;
    double tf;
};

int main()
{
    using Position = std::vector<std::vector<double>>;
    using Velocity = std::vector<std::vector<double>>;
    using Acceleration = std::vector<std::vector<double>>;
    std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};
    Position joint_position;
    Velocity joint_velocity;
    Acceleration joint_acceleration;

    std::vector<JointTrajectoryConstraints> joint_trajectory_constraints;
    Polynomial polynomialGenerator;
    // trajectory constraints for two link robot first trajectory point
    {
        JointTrajectoryConstraints joint1;
        joint1.initial_joint_constraints = {0.5, 0, 0};
        joint1.final_joint_constraints = {1.57, 0, 2};
        joint1.ti = 0;
        joint1.tf = 2;
        joint_trajectory_constraints.push_back(joint1);

        JointTrajectoryConstraints joint2;
        joint2.initial_joint_constraints = {-0.2, 0, 0};
        joint2.final_joint_constraints = {0.5, 1, 1};
        joint2.ti = 0;
        joint2.tf = 2;
        joint_trajectory_constraints.push_back(joint2);

        // JointTrajectoryConstraints joint3;
        // joint3.initial_joint_constraints = {-0.57, 0, 0};
        // joint3.final_joint_constraints = {0.5, 0, 0};
        // joint3.ti = 0;
        // joint3.tf = 1;
        // joint_trajectory_constraints.push_back(joint3);

        // JointTrajectoryConstraints joint4;
        // joint4.initial_joint_constraints = {-0.2, 0, 0};
        // joint4.final_joint_constraints = {1.1, 0, 0};
        // joint4.ti = 0;
        // joint4.tf = 1;
        // joint_trajectory_constraints.push_back(joint4);
    }

    std::vector<double> x = polynomialGenerator.TimeInterval(0, 2, 50);
    for (size_t i{0}; i < joint_trajectory_constraints.size(); ++i)
    {
        std::vector<double> coefficients = polynomialGenerator.GenerateQuanticPolynomial(joint_trajectory_constraints[i].initial_joint_constraints, joint_trajectory_constraints[i].final_joint_constraints, joint_trajectory_constraints[i].ti, joint_trajectory_constraints[i].tf);
        std::vector<double> pos = polynomialGenerator.PolynomialPointsFromCoefficients(coefficients, joint_trajectory_constraints[i].ti, joint_trajectory_constraints[i].tf, 50);
        joint_position.push_back(pos);
        std::vector<double> vel = polynomialGenerator.ComputeVelocity(coefficients, x);
        joint_velocity.push_back(vel);
        std::vector<double> ace = polynomialGenerator.ComputeAcceleration(coefficients, x);
        joint_acceleration.push_back(ace);
    }
    plt::figure();
    for (size_t joint{0}; joint < joint_position.size(); ++joint)
    {

        plt::title("Trajectory");
        std::map<std::string, std::string> legends;
        legends.insert(std::pair<std::string, std::string>("label", "Joint" + std::to_string(joint)));
        plt::plot(x, joint_position[joint], legends);
        plt::legend();
    }

    plt::figure();

    for (size_t joint{0}; joint < joint_position.size(); ++joint)
    {

        plt::title("Velocity");
        std::map<std::string, std::string> legends;
        legends.insert(std::pair<std::string, std::string>("label", "Joint" + std::to_string(joint)));
        plt::plot(x, joint_velocity[joint], legends);
        plt::legend();
    }
    plt::figure();

    for (size_t joint{0}; joint < joint_position.size(); ++joint)
    {

        plt::title("Acceleration");
        std::map<std::string, std::string> legends;
        legends.insert(std::pair<std::string, std::string>("label", "Joint" + std::to_string(joint)));
        plt::plot(x, joint_acceleration[joint], legends);
        plt::legend();
    }
    plt::show();
    // plt::savefig("standard.pdf"); // save the figure
    /*std::cout<<"Poly points are: [ ";
        for(const auto& polypoint:polypoints)
        {
            std::cout<<polypoint<<" , ";
        }
        std::cout<<" ]"<<std::endl;


        std::cout<<"Poly derivatives are: [ ";
        for(const auto& polypoint:polyderivatives)
        {
            std::cout<<polypoint<<" , ";
        }
        std::cout<<" ]"<<std::endl;

        std::cout<<"Poly second derivatives are: [ ";
        for(const auto& polypoint:polyderivativesAccel)
        {
            std::cout<<polypoint<<" , ";
        }
        std::cout<<" ]"<<std::endl;*/
    for (size_t joint{0}; joint < joint_position.size(); ++joint)
    {

        plt::figure();
        plt::title("position, velocity and acceleration");
        plt::plot(x, joint_position[joint], {{"label", "position "}});
        plt::legend();

        plt::plot(x, joint_velocity[joint], {{"label", "velocity "}});
        plt::legend();

        plt::plot(x, joint_acceleration[joint], {{"label", "acceleration "}});
        plt::legend();

        if (joint == joint_position.size() - 1)
            plt::show();
    }
    return 0;
}
