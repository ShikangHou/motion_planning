/*
 * @Author: HouShikang
 * @Date: 2024-04-26 09:12:57
 * @Description:
 */
#include "grid_map.h"
#include <iostream>
#include <ros/ros.h>

Eigen::Vector3d GridMap::gridIndex2coord(const Eigen::Vector3i &index) const
{
    return Eigen::Vector3d(coord_lower_(0) + resolution_ * (index(0) + 0.5),
                           coord_lower_(1) + resolution_ * (index(1) + 0.5),
                           coord_lower_(2) + resolution_ * (index(2) + 0.5));
}

Eigen::Vector3i GridMap::coord2gridIndex(const Eigen::Vector3d &coord) const
{

    return Eigen::Vector3i(
        std::min(std::max((int)floor((coord(0) - coord_lower_(0)) / resolution_), 0), max_index_(0)),
        std::min(std::max((int)floor((coord(1) - coord_lower_(1)) / resolution_), 0), max_index_(1)),
        std::min(std::max((int)floor((coord(2) - coord_lower_(2)) / resolution_), 0), max_index_(2)));
}

bool GridMap::isOutlier(const Eigen::Vector3d &coord) const
{
    if (coord(0) < coord_lower_(0) || coord(1) < coord_lower_(1) || coord(2) < coord_lower_(2) ||
        coord(0) > coord_upper_(0) || coord(1) > coord_upper_(1) || coord(2) > coord_upper_(2))
    {
        return true;
    }
    return false;
}

bool GridMap::isOutlier(const Eigen::Vector3i &index) const
{
    if (index(0) < 0 || index(1) < 0 || index(2) < 0 || index(0) >= x_size_ || index(1) >= y_size_ ||
        index(2) >= z_size_)
    {
        return true;
    }
    return false;
}

bool GridMap::isFree(const Eigen::Vector3i &index) const
{
    if (!isOutlier(index))
    {
        if (data_[index(0) * yz_size_ + index(1) * z_size_ + index(2)] == 0)
            return true;
    }
    return false;
}

bool GridMap::isFree(const Eigen::Vector3d &coord) const
{
    if (!isOutlier(coord))
    {
        Eigen::Vector3i index = coord2gridIndex(coord);
        if (data_[index(0) * yz_size_ + index(1) * z_size_ + index(2)] == 0)
            return true;
    }
    return false;
}

std::vector<Eigen::Vector3d> GridMap::getVisitedNodes()
{
    std::vector<Eigen::Vector3d> nodes;

    for (int i = 0; i < x_size_; i++)
    {
        for (int j = 0; j < y_size_; j++)
        {

            for (int k = 0; k < z_size_; k++)
            {
                if (gridNodeMap[i][j][k]->id_ != 0)
                    nodes.push_back(gridNodeMap[i][j][k]->coord_);
            }
        }
    }
    return nodes;
}

void GridMap::initGridMap(double resolution, const Eigen::Vector3d &coord_l, const Eigen::Vector3d &coord_u,
                          const Eigen::Vector3i &max_index)
{
    resolution_ = resolution;

    inv_resolution_ = 1.0 / resolution;
    coord_lower_ = coord_l;
    coord_upper_ = coord_u;
    max_index_ = max_index;

    x_size_ = max_index_(0) + 1; // 0--max_index_
    y_size_ = max_index_(1) + 1;
    z_size_ = max_index_(2) + 1;

    yz_size_ = y_size_ * z_size_;
    xyz_size_ = x_size_ * yz_size_;

    data_.resize(xyz_size_);
    std::fill(data_.begin(), data_.end(), 0);

    gridNodeMap.resize(x_size_);
    for (int i = 0; i < x_size_; i++)
    {
        gridNodeMap[i].resize(y_size_);
        for (int j = 0; j < y_size_; j++)
        {
            gridNodeMap[i][j].resize(z_size_);
            for (int k = 0; k < z_size_; k++)
            {
                Eigen::Vector3i index_(i, j, k);
                Eigen::Vector3d coord = gridIndex2coord(index_);
                gridNodeMap[i][j][k] = std::make_shared<GridNode>(GridNode(index_, coord));
            }
        }
    }
}

void GridMap::resetAllNodes()
{
    for (int i = 0; i < x_size_; i++)
    {

        for (int j = 0; j < y_size_; j++)
        {
            for (int k = 0; k < z_size_; k++)
            {
                gridNodeMap[i][j][k]->resetNode();
            }
        }
    }
}

GridNodePtr GridMap::getNodePtr(const Eigen::Vector3i &index)
{
    assert(!isOutlier(index));
    return gridNodeMap[index(0)][index(1)][index(2)];
}

void GridMap::setObs(const Eigen::Vector3d &coord_obs)
{
    // if (!isOutlier(coord_obs)) // 通过coord2gridIndex中的范围约束坐标点
    // {
        Eigen::Vector3i index_ = coord2gridIndex(coord_obs);
        data_[yz_size_ * index_(0) + z_size_ * index_(1) + index_(2)] = 1;
    // }
}
