/*
 * @Author: HouShikang
 * @Date: 2024-04-27 06:30:08
 * @Description:
 */
#include "kinodynamic_astar.h"
#include <Eigen/Eigen>
#include <ros/ros.h>

using namespace Eigen;
using namespace std;

// config grid map
void KinoAstarPathFinder::initGridMap(double resolution, const Eigen::Vector3d &coord_l, const Eigen::Vector3d &coord_u,
                                      const Eigen::Vector3i &max_index)
{
    grid_map_->initGridMap(resolution, coord_l, coord_u, max_index);
}
void KinoAstarPathFinder::setObs(double px, double py, double pz)
{
    grid_map_->setObs(Eigen::Vector3d(px, py, pz));
}
void KinoAstarPathFinder::resetAllNodes()
{
    grid_map_->resetAllNodes();
}
Eigen::Vector3d KinoAstarPathFinder::coordRounding(const Eigen::Vector3d &coord)
{
    return grid_map_->gridIndex2coord(grid_map_->coord2gridIndex(coord));
}

// kinodynamic A star
void KinoAstarPathFinder::expand(const GridNodePtr &node_ptr)
{
    expanded_nodes_.clear();
    cost_set_.clear();
    state_set_.clear();

    for (int i = 0; i <= discretize_step_; i++)
    {
        for (int j = 0; j <= discretize_step_; j++)
        {
            for (int k = 0; k <= discretize_step_; k++)
            {
                Vector3d input;
                input(0) = -max_acc_x_ + i * 2 * max_acc_x_ / (double)discretize_step_;
                input(1) = -max_acc_y_ + j * 2 * max_acc_y_ / (double)discretize_step_;
                input(2) = -max_acc_z_ + k * 2 * max_acc_z_ / (double)discretize_step_;

                double delta_t = time_interval_ / time_step_;

                bool safety_check = true;
                Vector3d node_pos, node_vel;
                Eigen::Matrix<double, 6, 1> current_state = node_ptr->state_;
                for (int i = 0; i < time_step_; i++)
                {
                    current_state.head(3) =
                        current_state.head(3) + delta_t * current_state.tail(3) + 0.5 * delta_t * delta_t * input;
                    current_state.tail(3) = current_state.tail(3) + delta_t * input;

                    node_pos = current_state.head(3);
                    node_vel = current_state.tail(3);

                    // check safety
                    if (!grid_map_->isFree(node_pos))
                    {
                        safety_check = false;
                        break;
                    }
                }

                if (safety_check == false)
                    continue;

                // check velocity
                // vf = v0 + at, v0符合要求，a为常数，只需检查vf
                if (abs(node_vel(0)) > max_vel_x_ || abs(node_vel(1)) > max_vel_y_ || abs(node_vel(2)) > max_vel_z_)
                    continue;

                // check not in the same voxel
                Vector3i index = grid_map_->coord2gridIndex(node_pos);
                if ((index - node_ptr->index_).norm() == 0)
                    continue;

                // ROS_WARN("%d %d %d", index(0), index(1), index(2));
                GridNodePtr node = grid_map_->getNodePtr(index);
                state_set_.push_back(current_state);
                expanded_nodes_.push_back(node);
                double cost = node_ptr->gScore_ + (input.squaredNorm() + weight_time_) * time_interval_;
                cost_set_.push_back(cost);
                // ROS_WARN("cost %f", cost);
            }
        }
    }
}


bool KinoAstarPathFinder::analytical_expansion(const GridNodePtr &node_ptr)
{
    oneShot_set_.clear();
    double optimal_T = optimal_time_to_goal_;
    Matrix<double, 6, 1> ds = (goal_ptr_->state_ - node_ptr->state_);
    ds.head(3) = ds.head(3) - node_ptr->state_.tail(3) * optimal_T; // delta_p = p1 - p0 - V0 * t
    Matrix<double, 6, 1> alpha;

    alpha.head(3) =
        -12.0 / (optimal_T * optimal_T * optimal_T) * ds.head(3) + 6.0 / (optimal_T * optimal_T) * ds.tail(3);
    alpha.tail(3) = 6.0 / (optimal_T * optimal_T) * ds.head(3) - 2.0 / optimal_T * ds.tail(3);

    double delta_t = time_interval_ / time_step_;

    Vector3d node_pos, node_vel;
    Eigen::Matrix<double, 6, 1> current_state;
    Eigen::Matrix<double, 6, 1> node_state = node_ptr->state_;

    for (double t = delta_t; t <= optimal_T; t += delta_t)
    {
        Vector3d input = t * alpha.head(3) + alpha.tail(3);
        // check accelerate
        if (abs(input(0)) > max_acc_x_ || abs(input(1)) > max_acc_y_ || abs(input(2)) > max_acc_z_)
        {
            // return false;
        }

        current_state.head(3) =
            node_state.head(3) + t * node_state.tail(3) + t * t * t / 6.0 * alpha.head(3) + 0.5 * t * t * alpha.tail(3);
        current_state.tail(3) = node_state.tail(3) + 0.5 * t * t * alpha.head(3) + t * alpha.tail(3);

        node_pos = current_state.head(3);
        node_vel = current_state.tail(3);

        // check safety
        if (!grid_map_->isFree(node_pos))
        {
            return false;
        }
        // check velocity
        if (abs(node_vel(0)) > max_vel_x_ || abs(node_vel(1)) > max_vel_y_ || abs(node_vel(2)) > max_vel_z_)
        {
            // return false;
        }
        oneShot_set_.push_back(current_state);
    }
    return true;
}

// double
double KinoAstarPathFinder::getHeuristic(GridNodePtr node1, GridNodePtr goal_node, double &optimal_time)
{
    Vector3d delta_pos = (goal_node->state_ - node1->state_).head(3);
    Vector3d v0 = node1->state_.tail(3);
    Vector3d v1 = goal_node->state_.tail(3);
    double a0 = -36 * delta_pos.squaredNorm() / weight_time_;
    double a1 = 24 * delta_pos.dot(v0 + v1) / weight_time_;
    double a2 = -4 * (v0.squaredNorm() + v1.squaredNorm() + v0.dot(v1)) / weight_time_;
    double a3 = 0;

    Matrix4d A;
    A << 0, 0, 0, -a0, 1, 0, 0, -a1, 0, 1, 0, -a2, 0, 0, 1, -a3;
    Vector4cd eigen_vector = A.eigenvalues();

    double optimal_cost = 1000000;
    for (int i = 0; i < 4; i++)
    {
        double T = eigen_vector(i).real();
        if (T < 0)
            continue;
        // ROS_WARN("%f",t);

        double cost = weight_time_ * (T - a2 / T - a1 / 2 / T / T - a0 / 3 / T / T / T);
        if (cost < optimal_cost)
        {
            optimal_time = T;
            optimal_cost = cost;
        }
    }
    return 1.0001 * optimal_cost;
}

void KinoAstarPathFinder::search(const Eigen::Vector3d &start_pt, const Eigen::Vector3d &start_vel,
                                 const Eigen::Vector3d &end_pt, const Eigen::Vector3d &end_vel)
{
    ros::Time time0 = ros::Time::now();
    discretize_step_ = 3;
    max_acc_x_ = 0.8;
    max_acc_y_ = 0.8;
    max_acc_z_ = 0.8;
    max_vel_x_ = 1.0;
    max_vel_y_ = 1.0;
    max_vel_z_ = 1.0;

    time_interval_ = 1;
    time_step_ = 4;
    weight_time_ = 1;

    Vector3i start_index = grid_map_->coord2gridIndex(start_pt);
    Vector3i end_index = grid_map_->coord2gridIndex(end_pt);
    Eigen::Vector3d goal_pt = grid_map_->gridIndex2coord(end_index);
    GridNodePtr start_ptr = grid_map_->getNodePtr(start_index);
    start_ptr->state_ << start_pt, start_vel;
    goal_ptr_ = grid_map_->getNodePtr(end_index);
    goal_ptr_->state_ << goal_pt, end_vel;

    ROS_WARN("goal: %f %f %f", goal_ptr_->state_(0), goal_ptr_->state_(1), goal_ptr_->state_(2));

    priority_queue<GridNodePtr, vector<GridNodePtr>, NodeComparator> empty;
    open_set_.swap(empty);

    oneShot_set_.clear();

    open_set_.push(start_ptr);
    start_ptr->id_ = 1;

    terminal_ptr_ = nullptr;

    bool isNear = false;

    while (!open_set_.empty())
    {
        current_ptr_ = open_set_.top();
        open_set_.pop();
        current_ptr_->id_ = -1;

        if ((current_ptr_->coord_ - goal_ptr_->coord_).norm() < 5 * grid_map_->resolution_)
        {
            isNear = true;
        }
        else
        {
            isNear = false;
        }

        if (isNear)
        {
            getHeuristic(current_ptr_, goal_ptr_, optimal_time_to_goal_);
            if (analytical_expansion(current_ptr_))
            {
                ROS_WARN("one shot %d %d %d", current_ptr_->index_(0), current_ptr_->index_(1),
                         current_ptr_->index_(2));
                terminal_ptr_ = current_ptr_;
                ros::Time time1 = ros::Time::now();
                ROS_INFO("Find path successfully! The path cost is %f m, the time cost is %f ms!",
                         current_ptr_->gScore_ * grid_map_->resolution_, 1000 * (time1 - time0).toSec());
                return;
            }
        }

        expand(current_ptr_);

        for (size_t i = 0; i < expanded_nodes_.size(); i++)
        {

            GridNodePtr expand_ptr = expanded_nodes_[i];
            double time;
            if (expand_ptr->id_ == 0)
            {
                expand_ptr->gScore_ = cost_set_[i];
                expand_ptr->fScore_ = expand_ptr->gScore_ + getHeuristic(expand_ptr, goal_ptr_, time);
                expand_ptr->parent_ = current_ptr_;
                expand_ptr->state_ = state_set_[i];
                expand_ptr->id_ = 1;
                open_set_.push(expand_ptr);
            }
            else if (expand_ptr->id_ == 1)
            {
                if (cost_set_[i] < expand_ptr->gScore_)
                {
                    expand_ptr->gScore_ = cost_set_[i];
                    expand_ptr->fScore_ = expand_ptr->gScore_ + getHeuristic(expand_ptr, goal_ptr_, time);
                    expand_ptr->parent_ = current_ptr_;
                    expand_ptr->state_ = state_set_[i];
                }
            }
        }
    }
    ROS_WARN("No path!");
}

std::vector<Eigen::Vector3d> KinoAstarPathFinder::getPath()
{
    std::vector<Eigen::Vector3d> path;

    GridNodePtr node = terminal_ptr_;
    while (node != nullptr)
    {
        // ROS_WARN("%f %f %f %f %f %f", node->state_.head(3).x(), node->state_.head(3).y(),
        // node->state_.head(3).z(),
        //          node->state_.tail(3).x(), node->state_.tail(3).y(), node->state_.tail(3).z());
        // ROS_WARN("%f %f %f ",node->state_(0),node->state_(1),node->state_(2));
        path.push_back(node->state_.head(3));
        node = node->parent_;
    }
    return path;
}

std::vector<Eigen::Vector3d> KinoAstarPathFinder::getOneShotPath()
{
    std::vector<Eigen::Vector3d> path;

    for (const auto &oneShot : oneShot_set_)
    {
        path.push_back(oneShot.head(3));
    }
    return path;
}

std::vector<Eigen::Vector3d> KinoAstarPathFinder::getVisitedNodes()
{
    return grid_map_->getVisitedNodes();
}