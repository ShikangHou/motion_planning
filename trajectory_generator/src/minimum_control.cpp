/*
 * @Author: HouShikang
 * @Date: 2024-05-06 13:50:50
 * @Description:
 */
#include "minimum_control.h"
#include <OsqpEigen/OsqpEigen.h>
#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;

int MinimumControl::factorial(int num)
{
    int result = 1;
    for (int i = num; i > 0; i--)
    {
        result *= i;
    }
    return result;
}

Eigen::VectorXd MinimumControl::timeAllocation(const Eigen::MatrixXd &path, bool isFixed, double max_vel,
                                               double max_acc)
{
    Eigen::VectorXd times = Eigen::VectorXd::Zero(path.rows() - 1);

    if (isFixed)
    {
        times.setOnes();
    }
    else
    {
        double t_acc = max_vel / max_acc;
        double distance_threshold = max_acc * t_acc * t_acc;
        for (size_t i = 1; i < path.rows(); i++)
        {
            double t = 0;
            double dist = (path.row(i) - path.row(i - 1)).norm();
            if (dist > distance_threshold)
            {
                t = 2 * t_acc + (dist - distance_threshold) / max_vel;
            }
            else
            {
                t = 2 * sqrt(dist / max_acc);
            }
            times(i - 1) = t;
        }
    }

    return times;
}

Eigen::MatrixXd MinimumControl::solveQPClosedForm(int order, const Eigen::MatrixXd &path, const Eigen::MatrixXd &vel,
                                                  Eigen::MatrixXd &acc, Eigen::VectorXd &time)
{
    Eigen::MatrixXd PolyCoeff;
    int poly_order = 2 * order - 1;
    int poly_param_num = poly_order + 1;
    int stage_num = time.size();

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(stage_num * poly_param_num, stage_num * poly_param_num);
    for (int i = 0; i < stage_num; i++)
    {
        double T = time(i);
        Eigen::MatrixXd sub_M = Eigen::MatrixXd::Zero(poly_param_num, poly_param_num);
        for (int j = 0; j < order; j++)
        {
            for (int k = j; k < poly_param_num; k++)
            {
                sub_M(j, k) = factorial(k) / factorial(k - j) * pow(0, k - j);
                sub_M(j + order, k) = factorial(k) / factorial(k - j) * pow(T, k - j);
            }
        }
        M.block(i * poly_param_num, i * poly_param_num, poly_param_num, poly_param_num) = sub_M;
    }

    // C_T
    int valid_num = order * (stage_num + 1);
    int fixed_num = stage_num + 2 * order - 1;
    Eigen::MatrixXd C_T = Eigen::MatrixXd::Zero(poly_param_num * stage_num, valid_num);
    C_T.topLeftCorner(order, order).setIdentity();
    C_T.block(0, 0, poly_param_num * stage_num, fixed_num).bottomRightCorner(order, order).setIdentity();
    for (int i = 0; i < stage_num - 1; i++)
    {
        C_T((2 * i + 1) * order, order + i) = 1;
        C_T((2 * i + 2) * order, order + i) = 1;
        for (int j = 1; j < order; j++)
        {
            C_T((2 * i + 1) * order + j, fixed_num + i * (order - 1) + j - 1) = 1;
            C_T((2 * i + 2) * order + j, fixed_num + i * (order - 1) + j - 1) = 1;
        }
    }

    // Q
    int derivatives_num = poly_param_num * stage_num;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(derivatives_num, derivatives_num);
    for (int i = 0; i < stage_num; i++)
    {
        double T = time(i);
        for (int j = order; j < poly_param_num; j++)
        {
            for (int k = order; k < poly_param_num; k++)
            {
                Q(i * poly_param_num + j, i * poly_param_num + k) = factorial(j) / factorial(j - order) * factorial(k) /
                                                                    factorial(k - order) / (j + k - poly_order) *
                                                                    pow(T, j + k - poly_order);
            }
        }
    }

    Eigen::MatrixXd d_fixed = Eigen::MatrixXd::Zero(fixed_num, 3);
    for (int i = 0; i < fixed_num; i++)
    {
        if (i == 0)
        {
            d_fixed.row(i) = path.row(0);
            continue;
        }
        if (i == 1 && order >= 2)
        {
            d_fixed.row(i) = vel.row(0);
            continue;
        }
        if (i == 2 && order >= 3)
        {
            d_fixed.row(i) = acc.row(0);
            continue;
        }
        if (i == fixed_num - order)
        {
            d_fixed.row(i) = path.bottomRows(1);
            continue;
        }
        if (i == fixed_num - order + 1)
        {
            d_fixed.row(i) = vel.bottomRows(1);
            continue;
        }
        if (i == fixed_num - order + 2)
        {
            d_fixed.row(i) = acc.bottomRows(1);
            continue;
        }
        if (i >= order && i < fixed_num - order)
        {
            d_fixed.row(i) = path.row(i - order + 1);
            continue;
        }
    }

    int varibles_num = valid_num - fixed_num;
    Eigen::MatrixXd R = C_T.transpose() * M.transpose().inverse() * Q * M.inverse() * C_T;
    Eigen::MatrixXd Rpp = R.bottomRightCorner(varibles_num, varibles_num);
    Eigen::MatrixXd Rpf = R.bottomLeftCorner(varibles_num, fixed_num);
    Eigen::MatrixXd d_p = -Rpp.inverse() * Rpf * d_fixed;
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(valid_num, 3);
    d.topRows(fixed_num) = d_fixed;
    d.bottomRows(varibles_num) = d_p;
    PolyCoeff = (M.inverse() * C_T * d).transpose();

    // std::cout << "M:" << endl << M << endl;
    // std::cout << "C_T:" << endl << C_T << endl;
    // std::cout << "Q:" << std::endl << Q << std::endl;
    // std::cout << "R:" << endl << R << endl;
    // std::cout << "polycoeff:" << std::endl << PolyCoeff << std::endl;

    return PolyCoeff.transpose();
}

Eigen::MatrixXd MinimumControl::solveQPNumerical(int order, const Eigen::MatrixXd &path, const Eigen::MatrixXd &vel,
                                                 Eigen::MatrixXd &acc, Eigen::VectorXd &time)
{

    const int poly_param_num = 2 * order;
    const int stage_num = time.size();
    const int variables_num = poly_param_num * stage_num;
    const int constraints_num = (stage_num + 1) * order + stage_num - 1;

    // Q
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(variables_num, variables_num);
    for (int i = 0; i < stage_num; i++)
    {
        double T = time(i);
        Eigen::MatrixXd subQ = Eigen::MatrixXd::Zero(poly_param_num, poly_param_num);
        for (int j = order; j < 2 * order; j++)
        {
            for (int k = order; k < 2 * order; k++)
            {
                subQ(j, k) = factorial(j) / factorial(j - order) * factorial(k) / factorial(k - order) /
                             (j + k - 2 * order + 1) * pow(T, j + k - 2 * order + 1);
            }
        }
        Q.block(i * poly_param_num, i * poly_param_num, poly_param_num, poly_param_num) = subQ;
    }

    // A
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(constraints_num, variables_num);
    const int derivative_num = 2 * order + stage_num - 1; // the number of derivative contraint
    const int continuity_num = (stage_num - 1) * order;   // the number of continuity contraint

    for (int i = 0; i < order; i++)
    {
        A(i, i) = factorial(i);
    }

    Eigen::MatrixXd Am = Eigen::MatrixXd::Zero(order, 2 * order);

    double T = time(time.size() - 1);
    for (int i = 0; i < order; i++)
    {
        for (int j = i; j < 2 * order; j++)
        {
            Am(i, j) = factorial(j) / factorial(j - i) * pow(T, j - i);
        }
    }
    A.block(derivative_num - order, variables_num - 2 * order, order, 2 * order) = Am;
    for (int i = 1; i < stage_num; i++)
    {
        A(order + i - 1, i * 2 * order) = 1;
    }
    for (int i = 0; i < stage_num - 1; i++)
    {
        double T = time(i);
        for (int j = 0; j < order; j++)
        {
            for (int k = j; k < 2 * order; k++)
            {
                A(derivative_num + i * order + j, i * 2 * order + k) = factorial(k) / factorial(k - j) * pow(T, k - j);
                A(derivative_num + i * order + j, (i + 1) * 2 * order + k) =
                    -factorial(k) / factorial(k - j) * pow(0, k - j);
            }
        }
    }

    Eigen::MatrixXd polyCoeff = Eigen::MatrixXd::Zero(variables_num, 3);

    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    solver.settings()->setVerbosity(false);
    solver.data()->setNumberOfVariables(variables_num);
    solver.data()->setNumberOfConstraints(constraints_num);
    Eigen::SparseMatrix<double> hession = Q.sparseView();
    solver.data()->setHessianMatrix(hession);
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(variables_num);
    solver.data()->setGradient(gradient);
    Eigen::SparseMatrix<double> linear_matrix = A.sparseView();
    solver.data()->setLinearConstraintsMatrix(linear_matrix);
    Eigen::VectorXd linear_d = Eigen::VectorXd::Zero(constraints_num);

    for (int dim = 0; dim < 3; dim++)
    {

        for (int i = 0; i < derivative_num; i++)
        {
            if (i == 0)
            {
                linear_d(0) = path(0, dim);
            }
            else if (i == 1 && order > 1)
            {
                linear_d(i) = vel(0, dim);
            }
            else if (i == 2 && order > 2)
            {
                linear_d(i) = acc(0, dim);
            }
            else if (i > 2 && order > i)
            {
                linear_d(i) = 0;
            }
            else if (i == derivative_num - order)
            {
                linear_d(i) = path(stage_num, dim);
            }
            else if (i == derivative_num - order + 1)
            {
                linear_d(i) = vel(1, dim);
            }
            else if (i == derivative_num - order + 2)
            {
                linear_d(i) = acc(1, dim);
            }
            else
                linear_d(i) = path(i - order + 1, dim);
        }

        if (dim == 0)
        {
            solver.data()->setLowerBound(linear_d);
            solver.data()->setUpperBound(linear_d);
            solver.initSolver();
        }
        else
        {
            solver.updateBounds(linear_d, linear_d);
        }

        solver.solve();

        Eigen::VectorXd QPsolution = solver.getSolution();
        polyCoeff.col(dim) = QPsolution;


        // std::cout << QPsolution.transpose() << std::endl;
        // std::cout << "Q:" << std::endl << Q << std::endl;

        // std::cout << "A" << A << std::endl;
        // std::cout << "linear" << linear_d.transpose() << std::endl;
    }

    // std::cout << "polyCoeff:" << std::endl << polyCoeff.transpose() << std::endl;

    return polyCoeff;
}
