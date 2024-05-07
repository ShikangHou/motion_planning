#include <Eigen/Eigen>
#include <cmath>
#include <fstream>
#include <geometry_msgs/Point.h>
#include <nav_msgs/Path.h>
#include <random>
#include <ros/console.h>
#include <ros/ros.h>
#include <vector>
#include <visualization_msgs/Marker.h>

#include "minimum_control.h"

using namespace std;
using namespace Eigen;

double visualization_traj_width;
double max_vel, max_acc;
int dev_order; // cost函数对应的导数阶数, =3: 最小化jerk =4: 最小化snap

ros::Subscriber way_pts_sub;
ros::Publisher waypoint_traj_vis_pub, waypoint_path_vis_pub;

int poly_coeff_num;
MatrixXd poly_coeff;
Vector3d start_position;
Vector3d start_velocity;
vector<Vector3d> path_set;

MinimumControl *minco;

void visWayPointTraj(MatrixXd polyCoeff, VectorXd time);

void visWayPointPath(MatrixXd path);

Vector3d getPosPoly(MatrixXd polyCoeff, int k, double t);

/*!
 * 订阅rviz发布的waypoints
 * @param wp
 */
void rcvWaypointsCallBack(const nav_msgs::Path &wp)
{
    if (wp.poses[0].pose.position.z > 8.0)
    {
        path_set.clear();
    }
    else
    {
        Vector3d pos(wp.poses[0].pose.position.x, wp.poses[0].pose.position.y, wp.poses[0].pose.position.z);
        path_set.push_back(pos);
        MatrixXd path = MatrixXd::Zero(path_set.size() + 1, 3);
        path.row(0) = start_position.transpose();
        for (int i = 0; i < path_set.size(); i++)
        {
            path.row(i + 1) = path_set[i].transpose();
        }
        VectorXd time = minco->timeAllocation(path,false,3,1);
        MatrixXd vel = MatrixXd::Zero(2, 3);
        vel.row(0) = start_velocity.transpose();
        MatrixXd acc = MatrixXd::Zero(2, 3);
        poly_coeff = minco->solveQPClosedForm(dev_order, path, vel, acc, time);
        visWayPointPath(path);
        visWayPointTraj(poly_coeff,time);
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "traj_node");
    ros::NodeHandle nh("~");

    nh.param("planning/max_vel", max_vel, 1.0); // 当前机器人能运行的最大速度
    nh.param("planning/max_acc", max_acc, 1.0); // 当前机器人能运行的最大加速度
    nh.param("planning/dev_order", dev_order, 3);
    nh.param("vis/vis_traj_width", visualization_traj_width, 0.15);

    //_poly_numID is the maximum order of polynomial
    poly_coeff_num = 2 * dev_order;

    path_set.clear();
    minco = new MinimumControl();

    // state of start point
    start_position(0) = 0;
    start_position(1) = 0;
    start_position(2) = 0;

    start_velocity(0) = 0.1;
    start_velocity(1) = 0.2;
    start_velocity(2) = 0.3;

    way_pts_sub = nh.subscribe("waypoints", 1, rcvWaypointsCallBack);

    waypoint_traj_vis_pub = nh.advertise<visualization_msgs::Marker>("vis_trajectory", 1);
    waypoint_path_vis_pub = nh.advertise<visualization_msgs::Marker>("vis_waypoint_path", 1);

    ros::Rate rate(100);
    while (ros::ok())
    {
        ros::spinOnce();
        rate.sleep();
    }
    return 0;
}

void visWayPointTraj(MatrixXd polyCoeff, VectorXd time)
{
    visualization_msgs::Marker _traj_vis;

    _traj_vis.header.stamp = ros::Time::now();
    _traj_vis.header.frame_id = "/map";

    _traj_vis.ns = "traj_node/trajectory_waypoints";
    _traj_vis.id = 0;
    _traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    _traj_vis.action = visualization_msgs::Marker::ADD;
    _traj_vis.scale.x = visualization_traj_width;
    _traj_vis.scale.y = visualization_traj_width;
    _traj_vis.scale.z = visualization_traj_width;
    _traj_vis.pose.orientation.x = 0.0;
    _traj_vis.pose.orientation.y = 0.0;
    _traj_vis.pose.orientation.z = 0.0;
    _traj_vis.pose.orientation.w = 1.0;

    _traj_vis.color.a = 1.0;
    _traj_vis.color.r = 1.0;
    _traj_vis.color.g = 0.0;
    _traj_vis.color.b = 0.0;

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();

    _traj_vis.points.clear();
    Vector3d pos;
    geometry_msgs::Point pt;

    for (int i = 0; i < time.size(); i++)
    {
        for (double t = 0.0; t < time(i); t += 0.01, count += 1)
        {
            pos = getPosPoly(polyCoeff, i, t);
            cur(0) = pt.x = pos(0);
            cur(1) = pt.y = pos(1);
            cur(2) = pt.z = pos(2);
            _traj_vis.points.push_back(pt);

            if (count)
                traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    waypoint_traj_vis_pub.publish(_traj_vis);
}

void visWayPointPath(MatrixXd path)
{
    visualization_msgs::Marker points, line_list;
    int id = 0;
    points.header.frame_id = line_list.header.frame_id = "/map";
    points.header.stamp = line_list.header.stamp = ros::Time::now();
    points.ns = line_list.ns = "wp_path";
    points.action = line_list.action = visualization_msgs::Marker::ADD;
    points.pose.orientation.w = line_list.pose.orientation.w = 1.0;
    points.pose.orientation.x = line_list.pose.orientation.x = 0.0;
    points.pose.orientation.y = line_list.pose.orientation.y = 0.0;
    points.pose.orientation.z = line_list.pose.orientation.z = 0.0;

    points.id = id;
    line_list.id = id;

    points.type = visualization_msgs::Marker::SPHERE_LIST;
    line_list.type = visualization_msgs::Marker::LINE_LIST;

    points.scale.x = 0.3;
    points.scale.y = 0.3;
    points.scale.z = 0.3;
    points.color.a = 1.0;
    points.color.r = 0.0;
    points.color.g = 0.0;
    points.color.b = 0.0;

    line_list.scale.x = 0.15;
    line_list.scale.y = 0.15;
    line_list.scale.z = 0.15;
    line_list.color.a = 1.0;

    line_list.color.r = 0.0;
    line_list.color.g = 1.0;
    line_list.color.b = 0.0;

    line_list.points.clear();

    for (int i = 0; i < path.rows(); i++)
    {
        geometry_msgs::Point p;
        p.x = path(i, 0);
        p.y = path(i, 1);
        p.z = path(i, 2);

        points.points.push_back(p);

        if (i < (path.rows() - 1))
        {
            geometry_msgs::Point p_line;
            p_line = p;
            line_list.points.push_back(p_line);
            p_line.x = path(i + 1, 0);
            p_line.y = path(i + 1, 1);
            p_line.z = path(i + 1, 2);
            line_list.points.push_back(p_line);
        }
    }

    waypoint_path_vis_pub.publish(points);
    waypoint_path_vis_pub.publish(line_list);
}

/*!
 * 求解第k个轨迹段t时刻对应的位置
 * @param polyCoeff 多项式系数矩阵
 * @param k 轨迹段序号
 * @param t 时刻
 * @return
 */
Vector3d getPosPoly(MatrixXd polyCoeff, int k, double t)
{
    Vector3d ret;

    for (int dim = 0; dim < 3; dim++)
    {
        VectorXd coeff = (polyCoeff.col(dim)).segment(k * poly_coeff_num, poly_coeff_num);
        VectorXd time = VectorXd::Zero(poly_coeff_num);

        for (int j = 0; j < poly_coeff_num; j++)
            if (j == 0)
                time(j) = 1.0;
            else
                time(j) = pow(t, j);

        double temp_pose = 0.0;
        for (int i = 0; i < time.rows(); ++i)
        {
            temp_pose = temp_pose + coeff(i) * time(i);
        }
        ret(dim) = temp_pose;
    }

    return ret;
}