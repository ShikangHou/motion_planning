/*
 * @Author: HouShikang
 * @Date: 2024-04-26 09:13:06
 * @Description:
 */
#ifndef __GRID_PATH_SEARCHER_GRID_MAP_h
#define __GRID_PATH_SEARCHER_GRID_MAP_h
#include <Eigen/Eigen>
#include <memory>

class GridNode;
typedef std::shared_ptr<GridNode> GridNodePtr;

class GridNode
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    int id_; // -1:close set, 0:unknown, 1:open set
    Eigen::Vector3i index_;
    Eigen::Vector3d coord_;
    // Eigen::Vector3i dir; // expand direction
    Eigen::Matrix<double, 6, 1> state_;

    GridNodePtr parent_;

    double gScore_;
    double fScore_;

    GridNode(Eigen::Vector3i index, Eigen::Vector3d coord)
    {
        id_ = 0;
        index_ = index;
        coord_ = coord;
        // dir = Eigen::Vector3i::Zero();
        state_ = Eigen::Matrix<double, 6, 1>::Zero();
        parent_ = nullptr;

        gScore_ = 0;
        fScore_ = 0;
    }

    GridNode()
    {
        id_ = 0;
        // dir = Eigen::Vector3i::Zero();
        state_ = Eigen::Matrix<double, 6, 1>::Zero();
        parent_ = nullptr;

        gScore_ = 0;
        fScore_ = 0;
    }

    void resetNode()
    {
        parent_ = nullptr;
        id_ = 0;
        gScore_ = 0;
        fScore_ = 0;
    }
    ~GridNode(){};
};

class GridMap
{
  public:
    std::vector<std::vector<std::vector<std::shared_ptr<GridNode>>>> gridNodeMap;

    std::vector<uint8_t> data_;
    Eigen::Vector3i max_index_;
    Eigen::Vector3d coord_lower_;
    Eigen::Vector3d coord_upper_;

    int x_size_, y_size_, z_size_, yz_size_, xyz_size_;
    double resolution_, inv_resolution_;
    GridMap(){};
    ~GridMap(){};

    void initGridMap(double resolution, const Eigen::Vector3d &coord_l, const Eigen::Vector3d &coord_u,
                     const Eigen::Vector3i &max_index);
    void setObs(const Eigen::Vector3d &coord_obs);
    void resetAllNodes();
    Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i &index) const;
    Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d &coord) const;
    GridNodePtr getNodePtr(const Eigen::Vector3i &index);
    bool isFree(const Eigen::Vector3i &index) const;
    bool isFree(const Eigen::Vector3d &coord) const;
    bool isOutlier(const Eigen::Vector3d &coord) const;
    bool isOutlier(const Eigen::Vector3i &index) const;
    std::vector<Eigen::Vector3d> getVisitedNodes();
};

#endif //__GRID_PATH_SEARCHER_GRID_MAP_h