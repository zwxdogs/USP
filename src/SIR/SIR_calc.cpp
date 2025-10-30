// Calculate Spatial Impulse Response (SIR).
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <vector>
#include <nanobind/stl/vector.h>
#include <algorithm>
#include <math.h>
#define _USE_MATH_DEFINES
#include <iostream>

namespace nb = nanobind;
using arr_np_c_in = nb::ndarray<double, nb::numpy, nb::c_contig, nb::shape<-1, 3>>;
using arr_np_c_out = nb::ndarray<double, nb::numpy, nb::c_contig>;
using vec_d = std::vector<double>;
using vec_i = std::vector<int>;
using vec_vec_d = std::vector<vec_d>;
using vec_vec_i = std::vector<vec_i>;

double
convert_angle(double angle, double y)
{
    // convert angle to [0, 2pi]
    double convert_angle = angle;
    if ( y < 0 ) {
        convert_angle = 2 * M_PI + angle;
    }
    return convert_angle;
}

vec_i
inside_demarcate(const arr_np_c_in& boundary_line, const arr_np_c_in& boundary_corners)
{
    // injudge which direction about the lines is inside the polygon.

    // calculate centroid of the polygon.
    double centroid_x = 0.0;
    double centroid_y = 0.0;
    for ( int i = 0; i < boundary_corners.shape(0); ++i ) {
        centroid_x += boundary_corners(i, 0);
        centroid_y += boundary_corners(i, 1);
    }
    centroid_x /= boundary_corners.shape(0);
    centroid_y /= boundary_corners.shape(0);

    vec_i inside_direction;
    // injudge and store.
    for ( int i = 0; i < boundary_line.shape(0); ++i ) {
        double line_a = boundary_line(i, 0);
        double line_b = boundary_line(i, 1);
        double line_c = boundary_line(i, 2);

        double injudge_value = line_a * centroid_x + line_b * centroid_y + line_c;
        if ( injudge_value > 0 ) {
            inside_direction.emplace_back(1);
        } else {
            inside_direction.emplace_back(-1);
        }
    }
    return inside_direction;
}

bool
point_inside_polygon(const vec_i& inside_direction,
                     const arr_np_c_in& boundary_line,
                     double x_point,
                     double y_point)
{
    // injudge whether the point is inside the polygon.
    for ( int i = 0; i < boundary_line.shape(0); ++i ) {
        double line_a = boundary_line(i, 0);
        double line_b = boundary_line(i, 1);
        double line_c = boundary_line(i, 2);

        double injudge_value = line_a * x_point + line_b * y_point + line_c;
        if ( injudge_value * inside_direction[i] <= 0 ) {
            return false;
        }
    }
    return true;
}

vec_i
update_active_edges(const arr_np_c_in& calc_line,
                    double x_p,
                    double y_p,
                    double z_p,
                    double dt,
                    double step,
                    double c0)
{
    // certain active edges at the discontinuites.
    vec_i active_edges;

    double r = c0 * step * dt;
    double r_proj_pow = r * r - z_p * z_p;
    for ( int i = 0; i < calc_line.shape(0); ++i ) {
        if ( calc_line(i, 1) == 0 ) {
            // k = inf
            double x_l = -calc_line(i, 2) / calc_line(i, 0);
            if ( std::pow(x_l - x_p, 2) <= r_proj_pow ) {
                active_edges.emplace_back(i);
            }
        } else {
            double line_k = -calc_line(i, 0) / calc_line(i, 1);
            double line_b = -calc_line(i, 2) / calc_line(i, 1);

            double eq_a = line_k * line_k + 1;
            double eq_b = 2 * line_k * line_b - 2 * x_p - 2 * line_k * y_p;
            double eq_c = y_p * y_p + x_p * x_p + line_b * line_b - 2 * line_b * y_p - r_proj_pow;
            double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;
            // Delta >= 0 means intersection exists.
            if ( eq_delta >= 0 ) {
                active_edges.emplace_back(i);
            }
        }
    }

    return active_edges;
}

vec_i
find_discontinuities(const arr_np_c_in& calc_line,
                     const arr_np_c_in& calc_corners,
                     double x_p,
                     double y_p,
                     double z_p,
                     double dt,
                     double c0)
{
    int N_dc = calc_corners.shape(0) + calc_line.shape(0);
    vec_i dc_step(N_dc, 0.0);
    // calculate all of the step (t / dt) for discontinuties points
    // for corners
    for ( int i = 0; i < calc_corners.shape(0); ++i ) {
        double length = std::sqrt(
          std::pow(calc_corners(i, 0) - x_p, 2)
          + std::pow(calc_corners(i, 1) - y_p, 2)
          + std::pow(z_p, 2));

        double tmp_t = length / c0;
        dc_step[i] = static_cast<int>(std::ceil(tmp_t / dt));
    }
    // for lines
    for ( int i = 0; i < calc_line.shape(0); ++i ) {
        // the minimum distance from scan point to the line.
        double tmp_t = 0;
        if ( calc_line(i, 1) == 0 ) {
            // k = inf
            double x_dc = -calc_line(i, 2) / calc_line(i, 0);
            double y_dc = y_p;
            double length = std::sqrt(
              std::pow(x_dc - x_p, 2)
              + std::pow(y_dc - y_p, 2)
              + std::pow(z_p, 2));

            tmp_t = length / c0;
            dc_step[calc_corners.shape(0) + i] = static_cast<int>(std::ceil(tmp_t / dt));
        } else {
            double line_k = -calc_line(i, 0) / calc_line(i, 1);
            double line_b = -calc_line(i, 2) / calc_line(i, 1);
            double x_dc = (line_k * y_p + x_p - line_k * line_b) / (line_k * line_k + 1);
            double y_dc = line_k * x_dc + line_b;
            double length = std::sqrt(
              std::pow(x_dc - x_p, 2)
              + std::pow(y_dc - y_p, 2)
              + std::pow(z_p, 2));

            tmp_t = length / c0;
            dc_step[calc_corners.shape(0) + i] = static_cast<int>(std::ceil(tmp_t / dt));
        }
    }

    return dc_step;
}

vec_d
circle_line_intersection(double x_p, double y_p, const arr_np_c_in& calc_line, int idx, double pow_r)
{
    // calculate intersection (convert to angle) points between circle and line.
    vec_d result;
    if ( calc_line(idx, 1) == 0 ) {
        // k = inf
        double x_is = -calc_line(idx, 2) / calc_line(idx, 0);
        // calculate y_is by quadratic equation.
        double eq_a = 1;
        double eq_b = -2 * y_p;
        double eq_c = y_p * y_p + std::pow(x_is - x_p, 2) - pow_r;
        double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;

        double y_is_1 = (-eq_b + std::sqrt(eq_delta)) / (2 * eq_a);
        double y_is_2 = (-eq_b - std::sqrt(eq_delta)) / (2 * eq_a);
        double angle_1 = convert_angle(std::atan2(y_is_1 - y_p, x_is - x_p), y_is_1 - y_p);
        double angle_2 = convert_angle(std::atan2(y_is_2 - y_p, x_is - x_p), y_is_2 - y_p);

        result.emplace_back(angle_1);
        result.emplace_back(angle_2);
    } else {
        double line_k = -calc_line(idx, 0) / calc_line(idx, 1);
        double line_b = -calc_line(idx, 2) / calc_line(idx, 1);
        // calculate x_is by quadratic equation.
        double eq_a = line_k * line_k + 1;
        double eq_b = 2 * line_k * line_b - 2 * x_p - 2 * line_k * y_p;
        double eq_c = y_p * y_p + x_p * x_p + line_b * line_b - 2 * line_b * y_p - pow_r;
        double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;

        double x_is_1 = (-eq_b + std::sqrt(eq_delta)) / (2 * eq_a);
        double x_is_2 = (-eq_b - std::sqrt(eq_delta)) / (2 * eq_a);
        double y_is_1 = line_k * x_is_1 + line_b;
        double y_is_2 = line_k * x_is_2 + line_b;
        double angle_1 = convert_angle(std::atan2(y_is_1 - y_p, x_is_1 - x_p), y_is_1 - y_p);
        double angle_2 = convert_angle(std::atan2(y_is_2 - y_p, x_is_2 - x_p), y_is_2 - y_p);

        result.emplace_back(angle_1);
        result.emplace_back(angle_2);
    }
    return result;
}

vec_vec_d
calc_polygon(arr_np_c_in boundary_line,
             arr_np_c_in calc_line,
             arr_np_c_in boundary_corners,
             arr_np_c_in calc_corners,
             arr_np_c_in scan_position,
             double dt,
             double c0)
{
    // injudge which direction about the lines is inside the polygon.
    vec_i inside_direction = inside_demarcate(boundary_line, boundary_corners);

    vec_vec_d SIR_result;
    // loop all of the scan positions
    for ( int i = 0; i < scan_position.shape(0); ++i ) {
        double x_p = scan_position(i, 0);
        double y_p = scan_position(i, 1);
        double z_p = scan_position(i, 2);

        // find all of the discontinuity points' time.
        vec_i dc_step = find_discontinuities(calc_line, calc_corners, x_p, y_p, z_p, dt, c0);
        // sort and remove repetition.
        std::sort(dc_step.begin(), dc_step.end());
        dc_step.erase(std::unique(dc_step.begin(), dc_step.end()), dc_step.end());

        // ----------------------------------------------------------------------------------------------------- //
        // calculate SIR.
        int max_step = dc_step.back() + 1;
        vec_d SIR_data(max_step, 0.0);

        int step_index = 0;  // index for calculate intersection steps.
        vec_i active_edges;  // store active edges index.

        for ( int j = dc_step.front(); j <= max_step; ++j ) {
            double pow_r = std::pow(c0 * (j * dt), 2) - std::pow(z_p, 2);
            // loop lines for calculate intersection points and angles.
            vec_d angles;  // Angles (band line) relative to x axis (origin is scan point).
            if ( j == dc_step[step_index] ) {
                // In discontinuities.
                active_edges = update_active_edges(calc_line, x_p, y_p, z_p, dt, j, c0);
                // loop active edges only.
                for ( auto idx : active_edges ) {
                    vec_d intersections_angle = circle_line_intersection(x_p, y_p, calc_line, idx, pow_r);
                    angles.emplace_back(intersections_angle[0]);
                    if ( intersections_angle[1] == intersections_angle[0] ) {
                        continue;
                    }
                    angles.emplace_back(intersections_angle[1]);
                }
                // sort angles.
                std::sort(angles.begin(), angles.end());
                angles.erase(std::unique(angles.begin(), angles.end()), angles.end());

                // find which point in the aperture and calculate sum of the angle difference.
                double angle_difference = 0;
                double N_outside = 0;
                vec_vec_i inside_angle_index;
                if ( angles.size() == 2 ) {
                    // special case: only two intersection points.
                    double angle_now = angles[0];
                    double angle_next = angles[1];
                    double angle_mid_1 = angle_now + M_PI / 2;
                    double point_x_mid_1 = x_p + std::cos(angle_mid_1) * std::sqrt(pow_r);
                    double point_y_mid_1 = y_p + std::sin(angle_mid_1) * std::sqrt(pow_r);

                    double angle_mid_2 = angle_next + M_PI / 2;
                    double point_x_mid_2 = x_p + std::cos(angle_mid_2) * std::sqrt(pow_r);
                    double point_y_mid_2 = y_p + std::sin(angle_mid_2) * std::sqrt(pow_r);
                    if ( point_inside_polygon(inside_direction, boundary_line, point_x_mid_1, point_y_mid_1) ) {
                        angle_difference += (angle_next - angle_now);
                    } else {
                        N_outside++;
                    }
                    if ( point_inside_polygon(inside_direction, boundary_line, point_x_mid_2, point_y_mid_2) ) {
                        angle_difference += (angle_next - angle_now);
                    } else {
                        N_outside++;
                    }

                } else {
                    for ( int k = 1; k < angles.size(); ++k ) {
                        // loop angles.
                        double angle_now = angles[k - 1];
                        double angle_next = angles[k];
                        double angle_mid = (angle_now + angle_next) / 2;
                        double point_x_mid = x_p + std::cos(angle_mid) * std::sqrt(pow_r);
                        double point_y_mid = y_p + std::sin(angle_mid) * std::sqrt(pow_r);

                        if ( point_inside_polygon(inside_direction, boundary_line, point_x_mid, point_y_mid) ) {
                            angle_difference += (angle_next - angle_now);
                            vec_i temp_inside_angle;
                            temp_inside_angle.emplace_back(k - 1);
                            temp_inside_angle.emplace_back(k);
                            inside_angle_index.emplace_back(temp_inside_angle);
                        } else {
                            N_outside++;
                        }
                    }
                }
                // all inside case.
                if ( N_outside == 0 ) {
                    angle_difference = 2 * M_PI;
                }
                // write data in this step.
                SIR_data[j - 1] = angle_difference / (2 * M_PI) * c0;
                step_index++;
            } else {
            }
        }
        // write data in this scan position.
        SIR_result.emplace_back(SIR_data);
    }

    return SIR_result;
}

NB_MODULE(SIR_calc, m)
{
    m.doc() = "Calculate Spatial Impulse Response (SIR) for polygon.";

    m.def("calc_polygon", &calc_polygon, "Calculate Spatial Impulse Response (SIR) for polygon.");
}