/*
 * @Author: HouShikang
 * @Date: 2024-04-27 06:30:22
 * @Description:
 */
#ifndef __KINODYNAMIC_ASTAR_H
#define __KINODYNAMIC_ASTAR_H
#include "grid_map.h"
#include <queue>

struct NodeComparator
{
    bool operator()(const GridNodePtr &ptr1, const GridNodePtr &ptr2)
    {
        return ptr1->fScore_ > ptr2->fScore_;
    }
};

class KinoAstarPathFinder
{
  private:
    std::unique_ptr<GridMap> grid_map_;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> open_set_;
    GridNodePtr current_ptr_;
    GridNodePtr goal_ptr_;
    GridNodePtr terminal_ptr_;
    std::vector<GridNodePtr> expanded_nodes_;
    std::vector<double> cost_set_;
    std::vector<Eigen::Matrix<double, 6, 1>> state_set_;
    std::vector<Eigen::Matrix<double, 6, 1>> oneShot_set_;

    int discretize_step_;
    double time_interval_, weight_time_;
    double optimal_time_to_goal_;
    int time_step_;
    double max_acc_x_, max_acc_y_, max_acc_z_;
    double max_vel_x_, max_vel_y_, max_vel_z_;

    void expand(const GridNodePtr &node_ptr);
    bool analytical_expansion(const GridNodePtr &node_ptr);
    double getHeuristic(GridNodePtr node1, GridNodePtr goal_node, double &optimal_time);

  public:
    KinoAstarPathFinder()
    {
        grid_map_ = std::unique_ptr<GridMap>(new GridMap);
    };
    ~KinoAstarPathFinder(){};
    void initGridMap(double resolution, const Eigen::Vector3d &coord_l, const Eigen::Vector3d &coord_u,
                     const Eigen::Vector3i &max_index);
    void setObs(double px, double py, double pz);
    void resetAllNodes();
    Eigen::Vector3d coordRounding(const Eigen::Vector3d &coord);

    void search(const Eigen::Vector3d &start_pt, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &end_pt,
                const Eigen::Vector3d &end_vel);
    std::vector<Eigen::Vector3d> getPath();
    std::vector<Eigen::Vector3d> getVisitedNodes();
    std::vector<Eigen::Vector3d> getOneShotPath();
};

#endif //__KINODYNAMIC_ASTAR_H