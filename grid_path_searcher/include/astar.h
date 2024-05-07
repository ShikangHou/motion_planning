#ifndef __GRID_PATH_SEARCHER_ASATR_H
#define __GRID_PATH_SEARCHER_ASATR_H
#include "grid_map.h"
#include <queue>
#include <vector>

struct cmp
{
    bool operator()(const GridNodePtr &ptr1, const GridNodePtr &ptr2)
    {
        return ptr1->fScore_ > ptr2->fScore_;
    }
};

class AstarPathFinder
{
protected:
    std::unique_ptr<GridMap> grid_map_;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, cmp> open_set_;
    GridNodePtr current_ptr_;
    GridNodePtr goal_ptr_;
    GridNodePtr terminal_ptr_;
    std::vector<GridNodePtr> expanded_nodes_;
    std::vector<double> cost_set_;

    void expand(const GridNodePtr &node_ptr);
    double getHeuristic(GridNodePtr node1, GridNodePtr node2);

public:
    AstarPathFinder() { grid_map_ = std::unique_ptr<GridMap>(new GridMap); };
    ~AstarPathFinder(){};
    
    void initGridMap(double resolution, const Eigen::Vector3d &coord_l, const Eigen::Vector3d &coord_u, const Eigen::Vector3i &max_index);
    void setObs(double px, double py, double pz);
    void resetAllNodes();
    Eigen::Vector3d coordRounding(const Eigen::Vector3d &coord);

    void search(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);
    std::vector<Eigen::Vector3d> getPath();
    std::vector<Eigen::Vector3d> getVisitedNodes();
};
#endif //__GRID_PATH_SEARCHER_ASATR_H