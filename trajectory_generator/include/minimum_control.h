/*
 * @Author: HouShikang
 * @Date: 2024-05-06 13:51:06
 * @Description:
 */
#ifndef __MINIMUN_CONTROL_H
#define __MINIMUN_CONTROL_H
#include <Eigen/Eigen>

class MinimumControl
{
  private:
    int factorial(int num);

  public:
    MinimumControl(){};
    ~MinimumControl(){};
    Eigen::VectorXd timeAllocation(const Eigen::MatrixXd &path,bool isFixed, double max_vel = 0, double max_acc = 0);
    Eigen::MatrixXd solveQPClosedForm(int order, const Eigen::MatrixXd &path, const Eigen::MatrixXd &vel,
                                      Eigen::MatrixXd &acc, Eigen::VectorXd &time);
};

#endif //__MINIMUN_CONTROL_H