//
// Created by caoqi on 2018/8/31.
//

// 3D:  1.36939, -1.17123, 7.04869
// obs: 0.180123 -0.156584

#include "sfm/bundle_adjustment.h"
/*
 * This function computes the Jacobian entries for the given camera and
 * 3D point pair that leads to one observation.
 *
 * The camera block 'cam_x_ptr' and 'cam_y_ptr' is:
 * - ID 0: Derivative of focal length f
 * - ID 1-2: Derivative of distortion parameters k0, k1
 * - ID 3-5: Derivative of translation t0, t1, t2
 * - ID 6-8: Derivative of rotation w0, w1, w2
 *
 * The 3D point block 'point_x_ptr' and 'point_y_ptr' is:
 * - ID 0-2: Derivative in x, y, and z direction.
 *
 * The function that leads to the observation is given as follows:
 *
 *   u = f * D(x,y) * x  (image observation x coordinate)
 *   v = f * D(x,y) * y  (image observation y coordinate)
 *
 * with the following definitions:
 *
 *   xc = R0 * X + t0  (homogeneous projection)
 *   yc = R1 * X + t1  (homogeneous projection)
 *   zc = R2 * X + t2  (homogeneous projection)
 *   x = xc / zc  (central projection)
 *   y = yc / zc  (central projection)
 *   D(x, y) = 1 + k0 (x^2 + y^2) + k1 (x^2 + y^2)^2  (distortion)
 */

/**
 * /description 给定一个相机参数和一个三维点坐标，求解雅各比矩阵，即公式中的df(theta)/dtheta
 * @param cam       相机参数
 * @param point     三维点坐标
 * @param cam_x_ptr 重投影坐标x 相对于相机参数的偏导数，相机有9个参数： [0] 焦距f; [1-2] 径向畸变系数k1, k2; [3-5] 平移向量 t1, t2, t3
 *                                                               [6-8] 旋转矩阵（角轴向量）
 * @param cam_y_ptr    重投影坐标y 相对于相机参数的偏导数，相机有9个参数
 * @param point_x_ptr  重投影坐标x 相对于三维点坐标的偏导数
 * @param point_y_ptr  重投影坐标y 相对于三维点坐标的偏导数
 */
void jacobian(sfm::ba::Camera const &cam,
              sfm::ba::Point3D const &point,
              double *cam_x_ptr, double *cam_y_ptr,
              double *point_x_ptr, double *point_y_ptr)
{
    const double f = cam.focal_length;
    const double *R = cam.rotation;
    const double *t = cam.translation;
    const double *X = point.pos;
    const double k0 = cam.distortion[0];
    const double k1 = cam.distortion[1];

    // TODO
    // 相机坐标系下的三维点坐标 Xc = R * X + t
    const double xc = R[0] * X[0] + R[1] * X[1] + R[2] * X[2] + t[0];
    const double yc = R[3] * X[0] + R[4] * X[1] + R[5] * X[2] + t[1];
    const double zc = R[6] * X[0] + R[7] * X[1] + R[8] * X[2] + t[2];

    const double x = xc / zc;
    const double y = yc / zc;
    const double r2 = x * x + y * y;
    const double r4 = r2 * r2;
    const double distort = 1.0 + (k0 * r2 + k1 * r4) * r2;
    double K[3][3] = {
        {f * distort, 0, 0},
        {0, f * distort, 0},
        {0, 0, 1}};

    // 重投影坐标 [u, v, 1]' = K * (R * X + t)
    const double u = f * distort * x;
    const double v = f * distort * y;

    // 投影像素坐标
    // u = f * ( 1 + k0 (x^2 + y^2) + k1 (x^2 + y^2)^2 ) * x + cx
    // v = f * ( 1 + k0 (x^2 + y^2) + k1 (x^2 + y^2)^2 ) * y + cy

    // 相机焦距的偏导数
    cam_x_ptr[0] = distort * x;
    cam_y_ptr[0] = distort * y;

    // 相机径向畸变的偏导数
    cam_x_ptr[1] = f * r2 * x;
    cam_x_ptr[2] = f * r4 * x;
    cam_y_ptr[1] = f * r2 * y;
    cam_y_ptr[2] = f * r4 * y;

    // 相机平移向量的偏导数
    // u = K * (R * X + t) [0]
    // v = K * (R * X + t) [1]

    const double x_deriv_xc = 1 / zc;
    const double x_deriv_yc = 0;
    const double x_deriv_zc = -x / zc;
    const double y_deriv_xc = 0;
    const double y_deriv_yc = 1 / zc;
    const double y_deriv_zc = -y / zc;

    // 计算u, v关于x, y的偏导数
    const double u_deriv_x = f * distort;
    const double v_deriv_y = f * distort;

    // 计算中间变量distort关于r2的偏导数
    const double distort_deriv_r2 = k0 + 2 * k1 * r2;

    // 计算中间变量r2关于xc, yc, zc的偏导数
    const double r2_deriv_xc = 2 * x / zc;
    const double r2_deriv_yc = 2 * y / zc;
    const double r2_deriv_zc = -2 * r2 / zc;

    // 计算中间变量distort关于xc, yc, zc的偏导数
    const double distort_deriv_xc = distort_deriv_r2 * r2_deriv_xc;
    const double distort_deriv_yc = distort_deriv_r2 * r2_deriv_yc;
    const double distort_deriv_zc = distort_deriv_r2 * r2_deriv_zc;

    // 计算(u,v)关于xc, yc, zc的偏导数
    const double u_deriv_xc = f * x * distort_deriv_xc + u_deriv_x * x_deriv_xc;
    const double u_deriv_yc = f * x * distort_deriv_yc + u_deriv_x * x_deriv_yc;
    const double u_deriv_zc = f * x * distort_deriv_zc + u_deriv_x * x_deriv_zc;

    const double v_deriv_xc = f * y * distort_deriv_xc + v_deriv_y * y_deriv_xc;
    const double v_deriv_yc = f * y * distort_deriv_yc + v_deriv_y * y_deriv_yc;
    const double v_deriv_zc = f * y * distort_deriv_zc + v_deriv_y * y_deriv_zc;

    /* 计算关于平移向量的t0, t1, t2的偏导数*/
    const double xc_deriv_t0 = 1;
    const double yc_deriv_t1 = 1;
    const double zc_deriv_t2 = 1;

    cam_x_ptr[3] = u_deriv_xc * xc_deriv_t0;
    cam_x_ptr[4] = u_deriv_yc * yc_deriv_t1;
    cam_x_ptr[5] = u_deriv_zc * zc_deriv_t2;
    cam_y_ptr[3] = v_deriv_xc * xc_deriv_t0;
    cam_y_ptr[4] = v_deriv_yc * yc_deriv_t1;
    cam_y_ptr[5] = v_deriv_zc * zc_deriv_t2;

    // 相机旋转矩阵的偏导数
    const double rx = R[0] * X[0] + R[1] * X[1] + R[2] * X[2];
    const double ry = R[3] * X[0] + R[4] * X[1] + R[5] * X[2];
    const double rz = R[6] * X[0] + R[7] * X[1] + R[8] * X[2];
    const double xc_deriv_w0 = 0;
    const double xc_deriv_w1 = rz;
    const double xc_deriv_w2 = -ry;
    const double yc_deriv_w0 = -rz;
    const double yc_deriv_w1 = 0;
    const double yc_deriv_w2 = rx;
    const double zc_deriv_w0 = ry;
    const double zc_deriv_w1 = -rx;
    const double zc_deriv_w2 = 0;

    cam_x_ptr[6] = u_deriv_yc * yc_deriv_w0 + u_deriv_zc * zc_deriv_w0;
    cam_x_ptr[7] = u_deriv_xc * xc_deriv_w1 + u_deriv_zc * zc_deriv_w1;
    cam_x_ptr[8] = u_deriv_xc * xc_deriv_w2 + u_deriv_yc * yc_deriv_w2;
    cam_y_ptr[6] = v_deriv_yc * yc_deriv_w0 + v_deriv_zc * zc_deriv_w0;
    cam_y_ptr[7] = v_deriv_xc * xc_deriv_w1 + v_deriv_zc * zc_deriv_w1;
    cam_y_ptr[8] = v_deriv_xc * xc_deriv_w2 + v_deriv_yc * yc_deriv_w2;

    // 三维点的偏导数
    const double xc_deriv_X = R[0];
    const double xc_deriv_Y = R[1];
    const double xc_deriv_Z = R[2];
    const double yc_deriv_X = R[3];
    const double yc_deriv_Y = R[4];
    const double yc_deriv_Z = R[5];
    const double zc_deriv_X = R[6];
    const double zc_deriv_Y = R[7];
    const double zc_deriv_Z = R[8];

    point_x_ptr[0] = u_deriv_xc * xc_deriv_X + u_deriv_yc * yc_deriv_X + u_deriv_zc * zc_deriv_X;
    point_x_ptr[1] = u_deriv_xc * xc_deriv_Y + u_deriv_yc * yc_deriv_Y + u_deriv_zc * zc_deriv_Y;
    point_x_ptr[2] = u_deriv_xc * xc_deriv_Z + u_deriv_yc * yc_deriv_Z + u_deriv_zc * zc_deriv_Z;
    point_y_ptr[0] = v_deriv_xc * xc_deriv_X + v_deriv_yc * yc_deriv_X + v_deriv_zc * zc_deriv_X;
    point_y_ptr[1] = v_deriv_xc * xc_deriv_Y + v_deriv_yc * yc_deriv_Y + v_deriv_zc * zc_deriv_Y;
    point_y_ptr[2] = v_deriv_xc * xc_deriv_Z + v_deriv_yc * yc_deriv_Z + v_deriv_zc * zc_deriv_Z;
}

int main(int argc, char *argv[])
{

    sfm::ba::Camera cam;
    cam.focal_length = 0.919654;
    cam.distortion[0] = -0.108298;
    cam.distortion[1] = 0.103775;

    cam.rotation[0] = 0.999999;
    cam.rotation[1] = -0.000676196;
    cam.rotation[2] = -0.0013484;
    cam.rotation[3] = 0.000663243;
    cam.rotation[4] = 0.999949;
    cam.rotation[5] = -0.0104095;
    cam.rotation[6] = 0.00135482;
    cam.rotation[7] = 0.0104087;
    cam.rotation[8] = 0.999949;

    cam.translation[0] = 0.00278292;
    cam.translation[1] = 0.0587996;
    cam.translation[2] = -0.127624;

    sfm::ba::Point3D pt3D;
    pt3D.pos[0] = 1.36939;
    pt3D.pos[1] = -1.17123;
    pt3D.pos[2] = 7.04869;

    double cam_x_ptr[9] = {0};
    double cam_y_ptr[9] = {0};
    double point_x_ptr[3] = {0};
    double point_y_ptr[3] = {0};

    jacobian(cam, pt3D, cam_x_ptr, cam_y_ptr, point_x_ptr, point_y_ptr);

    std::cout << "Result is :" << std::endl;
    std::cout << "cam_x_ptr: ";
    for (int i = 0; i < 9; i++)
    {
        std::cout << cam_x_ptr[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "cam_y_ptr: ";
    for (int i = 0; i < 9; i++)
    {

        std::cout << cam_y_ptr[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "point_x_ptr: ";
    std::cout << point_x_ptr[0] << " " << point_x_ptr[1] << " " << point_x_ptr[2] << std::endl;

    std::cout << "point_y_ptr: ";
    std::cout << point_y_ptr[0] << " " << point_y_ptr[1] << " " << point_y_ptr[2] << std::endl;

    std::cout << "\nResult should be :\n"
              << "cam_x_ptr: 0.195942 0.0123983 0.000847141 0.131188 0.000847456 -0.0257388 0.0260453 0.95832 0.164303\n"
              << "cam_y_ptr: -0.170272 -0.010774 -0.000736159 0.000847456 0.131426 0.0223669 -0.952795 -0.0244697 0.179883\n"
              << "point_x_ptr: 0.131153 0.000490796 -0.0259232\n"
              << "point_y_ptr: 0.000964926 0.131652 0.0209965\n";

    return 0;
}
