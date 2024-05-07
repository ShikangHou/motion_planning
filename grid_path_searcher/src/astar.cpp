/*
 * @Author: HouShikang
 * @Date: 2024-04-26 12:41:16
 * @Description:
 */
#include "astar.h"
#include <iostream>
#include <ros/ros.h>

void AstarPathFinder::initGridMap(double resolution, const Eigen::Vector3d &coord_l, const Eigen::Vector3d &coord_u,
                                  const Eigen::Vector3i &max_index)
{
    grid_map_->initGridMap(resolution, coord_l, coord_u, max_index);
}

void AstarPathFinder::setObs(double px, double py, double pz)
{
    grid_map_->setObs(Eigen::Vector3d(px, py, pz));
}

void AstarPathFinder::resetAllNodes()
{
    grid_map_->resetAllNodes();
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d &coord)
{
    return grid_map_->gridIndex2coord(grid_map_->coord2gridIndex(coord));
}

void AstarPathFinder::expand(const GridNodePtr &node_ptr)
{
    expanded_nodes_.clear();
    cost_set_.clear();
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            for (int dz = -1; dz <= 1; dz++)
            {
                if (dx == 0 && dy == 0 && dz == 0)
                    continue;

                Eigen::Vector3i index = node_ptr->index_ + Eigen::Vector3i(dx, dy, dz);
                if (!grid_map_->isFree(index))
                    continue;

                GridNodePtr neiborPtr = grid_map_->getNodePtr(index);
                if (neiborPtr->id_ == -1) // already in closeSet
                    continue;

                double cost = node_ptr->gScore_ + sqrt(dx * dx + dy * dy + dz * dz);

                cost_set_.push_back(cost);
                expanded_nodes_.push_back(neiborPtr);
            }
        }
    }
}

double AstarPathFinder::getHeuristic(GridNodePtr node1, GridNodePtr node2)
{
    double heuristic;
    double dx = abs(node1->index_.x() - node2->index_.x());
    double dy = abs(node1->index_.y() - node2->index_.y());
    double dz = abs(node1->index_.z() - node2->index_.z());

    // Diagonal
    std::vector<double> xyz{dx, dy, dz};
    sort(xyz.begin(), xyz.end());
    heuristic = (sqrt(3) - sqrt(2)) * xyz[0] + (sqrt(2) - 1) * xyz[1] + xyz[2];

    // Euclidean
    // heuristic = sqrt(dx * dx + dy * dy + dz * dz);

    // manhattan 在某些情况下不一定满足最优性
    // heuristic = dx + dy + dz;

    // Dijkstra
    // heuristic = 0.;

    // tie breaker
    // heuristic *= 1.0001;

    return heuristic;
}

void AstarPathFinder::search(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt)
{
    ros::Time time = ros::Time::now();
    Eigen::Vector3i start_idx = grid_map_->coord2gridIndex(start_pt);
    Eigen::Vector3i end_idx = grid_map_->coord2gridIndex(end_pt);

    GridNodePtr start_node_ptr = grid_map_->getNodePtr(start_idx);
    goal_ptr_ = grid_map_->getNodePtr(end_idx);

    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, cmp> empty;
    open_set_.swap(empty);

    open_set_.push(start_node_ptr);
    start_node_ptr->id_ = 1;

    terminal_ptr_ = nullptr;

    if (!grid_map_->isFree(end_idx))
    {
        ROS_WARN("The goal is not free!");
        return;
    }

    while (!open_set_.empty())
    {
        current_ptr_ = open_set_.top();
        open_set_.pop();
        current_ptr_->id_ = -1;
        if (current_ptr_ == goal_ptr_)
        {
            terminal_ptr_ = goal_ptr_;
            ros::Time time1 = ros::Time::now();
            ROS_INFO("Find path successfully!");
            ROS_WARN("The path cost is %f m, the time cost is %f ms", goal_ptr_->gScore_ * grid_map_->resolution_,
                     (time1 - time).toSec() * 1000);
            return;
        }

        expand(current_ptr_);

        for (size_t i = 0; i < expanded_nodes_.size(); i++)
        {
            GridNodePtr neiborPtr = expanded_nodes_[i];
            double hScore = getHeuristic(neiborPtr, goal_ptr_);
            if (neiborPtr->id_ == 0) // have not been visited
            {
                neiborPtr->gScore_ = cost_set_[i];
                neiborPtr->fScore_ = neiborPtr->gScore_ + hScore;
                neiborPtr->id_ = 1;
                neiborPtr->parent_ = current_ptr_;
                open_set_.push(neiborPtr);
            }
            else if (neiborPtr->id_ == 1) // in open_set_
            {
                if (cost_set_[i] < neiborPtr->gScore_)
                {
                    neiborPtr->gScore_ = cost_set_[i];
                    neiborPtr->fScore_ = neiborPtr->gScore_ + hScore;
                    neiborPtr->parent_ = current_ptr_;
                }
            }
        }
    }

    ROS_INFO("No path!");
}

std::vector<Eigen::Vector3d> AstarPathFinder::getPath()
{
    std::vector<Eigen::Vector3d> path;

    GridNodePtr node = terminal_ptr_;
    while (node != nullptr)
    {
        path.push_back(node->coord_);
        node = node->parent_;
    }
    return path;
}

std::vector<Eigen::Vector3d> AstarPathFinder::getVisitedNodes()
{
    return grid_map_->getVisitedNodes();
}