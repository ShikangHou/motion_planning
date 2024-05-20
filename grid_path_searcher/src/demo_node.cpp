#include <fstream>
#include <iostream>
#include <math.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>

#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include "astar.h"
#include "kinodynamic_astar.h"

using namespace std;
using namespace Eigen;

// simulation param from launch file
double _resolution, _inv_resolution, _cloud_margin;
double _x_size, _y_size, _z_size;

// useful global variables
bool _has_map = false;

bool kinodynamic;

Vector3d _start_pt, _start_vel;
Vector3d _map_lower, _map_upper;
Vector3i _map_max_id;

// ros related
ros::Subscriber _map_sub, _pts_sub;
ros::Publisher _grid_path_vis_pub, _visited_nodes_vis_pub, _grid_map_vis_pub;

AstarPathFinder *_astar_path_finder;
KinoAstarPathFinder *_kino_astar_path_finder;

void rcvWaypointsCallback(const nav_msgs::Path &wp);
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map);

void visGridPath(vector<Vector3d> nodes);
void visOneShotPath(vector<Vector3d> nodes);
void visVisitedNode(vector<Vector3d> nodes);
void pathFinding(const Vector3d start_pt, const Vector3d target_pt);

void rcvWaypointsCallback(const nav_msgs::Path &wp)
{
    if (wp.poses[0].pose.position.z < 0.0 || _has_map == false)
        return;

    Vector3d target_pt;
    target_pt << wp.poses[0].pose.position.x, wp.poses[0].pose.position.y, wp.poses[0].pose.position.z;

    ROS_INFO("[node] receive the planning target");
    pathFinding(_start_pt, target_pt);
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map)
{
    if (_has_map)
        return;

    pcl::PointCloud<pcl::PointXYZ> cloud;
    pcl::PointCloud<pcl::PointXYZ> cloud_vis;
    sensor_msgs::PointCloud2 map_vis;

    pcl::fromROSMsg(pointcloud_map, cloud);

    if ((int)cloud.points.size() == 0)
        return;

    pcl::PointXYZ pt;

    if (kinodynamic)
    {
        for (int idx = 0; idx < (int)cloud.points.size(); idx++)
        {
            pt = cloud.points[idx];

            // set obstalces into grid map for path planning
            _kino_astar_path_finder->setObs(pt.x, pt.y, pt.z);

            // for visualize only
            Vector3d cor_round = _kino_astar_path_finder->coordRounding(Vector3d(pt.x, pt.y, pt.z));
            pt.x = cor_round(0);
            pt.y = cor_round(1);
            pt.z = cor_round(2);
            cloud_vis.points.push_back(pt);
        }
    }
    else
    {
        for (int idx = 0; idx < (int)cloud.points.size(); idx++)
        {
            pt = cloud.points[idx];

            // set obstalces into grid map for path planning
            _astar_path_finder->setObs(pt.x, pt.y, pt.z);

            // for visualize only
            Vector3d cor_round = _astar_path_finder->coordRounding(Vector3d(pt.x, pt.y, pt.z));
            pt.x = cor_round(0);
            pt.y = cor_round(1);
            pt.z = cor_round(2);
            cloud_vis.points.push_back(pt);
        }
    }

    cloud_vis.width = cloud_vis.points.size();
    cloud_vis.height = 1;
    cloud_vis.is_dense = true;

    pcl::toROSMsg(cloud_vis, map_vis);

    map_vis.header.frame_id = "/world";
    _grid_map_vis_pub.publish(map_vis);

    _has_map = true;
}

void pathFinding(const Vector3d start_pt, const Vector3d target_pt)
{
    std::vector<Eigen::Vector3d> grid_path, visited_nodes;
    // Call A* to search for a path
    if (kinodynamic)
    {
        _kino_astar_path_finder->search(start_pt, _start_vel, target_pt, Eigen::Vector3d(0, 0, 0));
        grid_path = _kino_astar_path_finder->getPath();
        visited_nodes = _kino_astar_path_finder->getVisitedNodes();
        auto oneShot = _kino_astar_path_finder->getOneShotPath();
        _kino_astar_path_finder->resetAllNodes();
        visOneShotPath(oneShot);
    }
    else
    {
        _astar_path_finder->search(start_pt, target_pt);
        grid_path = _astar_path_finder->getPath();
        visited_nodes = _astar_path_finder->getVisitedNodes();
        _astar_path_finder->resetAllNodes();
    }

    visGridPath(grid_path);
    visVisitedNode(visited_nodes);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "demo_node");
    ros::NodeHandle nh("~");

    nh.param("map/cloud_margin", _cloud_margin, 0.0);
    nh.param("map/resolution", _resolution, 0.2);

    nh.param("map/x_size", _x_size, 50.0);
    nh.param("map/y_size", _y_size, 50.0);
    nh.param("map/z_size", _z_size, 5.0);

    nh.param("planning/start_x", _start_pt(0), 0.0);
    nh.param("planning/start_y", _start_pt(1), 0.0);
    nh.param("planning/start_z", _start_pt(2), 0.0);
    nh.param("planning/start_vx", _start_vel(0), 0.0);
    nh.param("planning/start_vy", _start_vel(1), 0.0);
    nh.param("planning/start_vz", _start_vel(2), 0.0);
    nh.param("kinodynamic", kinodynamic, true);

    _map_lower << -_x_size / 2.0, -_y_size / 2.0, 0.0;
    _map_upper << +_x_size / 2.0, +_y_size / 2.0, _z_size;

    _inv_resolution = 1.0 / _resolution;

    _map_max_id << (int)(_x_size * _inv_resolution) - 1, (int)(_y_size * _inv_resolution) - 1,
        (int)(_z_size * _inv_resolution) - 1;

    _astar_path_finder = new AstarPathFinder();
    _astar_path_finder->initGridMap(_resolution, _map_lower, _map_upper, _map_max_id);

    _kino_astar_path_finder = new KinoAstarPathFinder();
    _kino_astar_path_finder->initGridMap(_resolution, _map_lower, _map_upper, _map_max_id);

    _map_sub = nh.subscribe("map", 1, rcvPointCloudCallBack);
    _pts_sub = nh.subscribe("waypoints", 1, rcvWaypointsCallback);

    _grid_map_vis_pub = nh.advertise<sensor_msgs::PointCloud2>("grid_map_vis", 1);
    _grid_path_vis_pub = nh.advertise<visualization_msgs::Marker>("grid_path_vis", 1);
    _visited_nodes_vis_pub = nh.advertise<visualization_msgs::Marker>("visited_nodes_vis", 1);

    ros::Rate rate(100);
    while (ros::ok())
    {
        ros::spinOnce();
        rate.sleep();
    }

    delete _astar_path_finder;
    delete _kino_astar_path_finder;
    return 0;
}

void visGridPath(vector<Vector3d> nodes)
{
    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "demo_node/astar_path";

    node_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;

    node_vis.color.a = 0.8;
    node_vis.color.r = 1.0;
    node_vis.color.g = 0.0;
    node_vis.color.b = 0.0;

    node_vis.scale.x = _resolution / 2;
    node_vis.scale.y = _resolution / 2;
    node_vis.scale.z = _resolution / 2;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        node_vis.points.push_back(pt);
    }

    _grid_path_vis_pub.publish(node_vis);
}

void visOneShotPath(vector<Vector3d> nodes)
{
    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "demo_node/oneShot";

    node_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;

    node_vis.color.a = 0.8;
    node_vis.color.r = 0.0;
    node_vis.color.g = 1.0;
    node_vis.color.b = 0.0;

    node_vis.scale.x = _resolution / 2;
    node_vis.scale.y = _resolution / 2;
    node_vis.scale.z = _resolution / 2;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        node_vis.points.push_back(pt);
    }

    _grid_path_vis_pub.publish(node_vis);
}

void visVisitedNode(vector<Vector3d> nodes)
{
    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "demo_node/expanded_nodes";
    node_vis.type = visualization_msgs::Marker::CUBE_LIST;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;
    node_vis.color.a = 0.5;
    node_vis.color.r = 0.0;
    node_vis.color.g = 0.0;
    node_vis.color.b = 1.0;

    node_vis.scale.x = _resolution;
    node_vis.scale.y = _resolution;
    node_vis.scale.z = _resolution;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        node_vis.points.push_back(pt);
    }

    _visited_nodes_vis_pub.publish(node_vis);
}