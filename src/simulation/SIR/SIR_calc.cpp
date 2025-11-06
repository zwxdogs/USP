// Calculate Spatial Impulse Response (SIR).
// line store: y - kx + b and x = c. use array[N_lines, 3]. First column: 0 or 1 for define whether the slope is infinte;
// second column: k or 0, depending on slope; third column: b or c, depending on slope.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <vector>
#include <nanobind/stl/vector.h>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <set>
#include <cfloat>
#include <iostream>
#include "../apodization.h"

// #define EPS_ANGLE 1e-10

namespace nb = nanobind;
using arr_np_c_in = nb::ndarray<double, nb::numpy, nb::c_contig, nb::shape<-1, 3>>;
using arr_np_c_out = nb::ndarray<double, nb::numpy, nb::c_contig>;
using vec_d = std::vector<double>;
using vec_i = std::vector<int>;
using set_i = std::set<int>;
using set_d = std::set<double>;
using vec_set_i = std::vector<set_i>;

double
convert_tan2angle(double angle)
{
    // convert angle to [0, 2pi]
    double convert_angle = angle;
    if ( angle < 0 ) {
        convert_angle = 2 * M_PI + angle;
    }
    return convert_angle;
}

void
inside_demarcate(const arr_np_c_in& boundary_line, const arr_np_c_in& boundary_corners, int* inside_direction)
{
    // Injudge which direction about the lines is inside the polygon.

    // Calculate centroid of the polygon.
    // Because ultrasound aperture is square in my subject, so this code just uses average to calculate centroid.
    auto* line_data = boundary_line.data();
    auto* corner_data = boundary_corners.data();

    double centroid_x = 0.0;
    double centroid_y = 0.0;
    for ( int i = 0; i < boundary_corners.shape(0); ++i ) {
        centroid_x += corner_data[i * 3 + 0];
        centroid_y += corner_data[i * 3 + 1];
    }
    centroid_x /= boundary_corners.shape(0);
    centroid_y /= boundary_corners.shape(0);

    // injudge and store.
    for ( int i = 0; i < boundary_line.shape(0); ++i ) {
        double injudge_value = 0.0;
        if ( static_cast<int>(line_data[i * 3 + 0] == 1) ) {
            injudge_value = centroid_x - line_data[i * 3 + 2];
        } else {
            injudge_value = line_data[i * 3 + 1] * centroid_x + line_data[i * 3 + 2] - centroid_y;
        }

        inside_direction[i] = (injudge_value > 0) ? 1 : -1;
    }
}

bool
point_inside_polygon(const int* inside_direction,
                     const arr_np_c_in& boundary_line,
                     double x_point,
                     double y_point)
{
    auto* line_data = boundary_line.data();
    // injudge whether the point is inside the polygon.
    for ( int i = 0; i < boundary_line.shape(0); ++i ) {
        double injudge_value = 0.0;
        if ( static_cast<int>(line_data[i * 3 + 0] == 1) ) {
            injudge_value = x_point - line_data[i * 3 + 2];
        } else {
            injudge_value = line_data[i * 3 + 1] * x_point + line_data[i * 3 + 2] - y_point;
        }

        if ( injudge_value * inside_direction[i] < 0 ) {
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
    auto* line_data = calc_line.data();

    // certain active edges at the discontinuites.
    vec_i active_edges;

    double r = c0 * step * dt;
    double r_proj_pow = r * r - z_p * z_p;
    for ( int i = 0; i < calc_line.shape(0); ++i ) {
        // line parameters: y = kx + b or x = c
        if ( static_cast<int>(line_data[i * 3 + 0] == 1) ) {
            // k = inf
            double x_l = line_data[i * 3 + 2];
            if ( ((x_l - x_p) * (x_l - x_p) < r_proj_pow) || fabs((x_l - x_p) * (x_l - x_p) - r_proj_pow) < DBL_EPSILON ) {
                active_edges.emplace_back(i);
            }
        } else {
            double line_k = line_data[i * 3 + 1];
            double line_b = line_data[i * 3 + 2];

            double eq_a = line_k * line_k + 1;
            double eq_b = 2 * line_k * line_b - 2 * x_p - 2 * line_k * y_p;
            double eq_c = y_p * y_p + x_p * x_p + line_b * line_b - 2 * line_b * y_p - r_proj_pow;
            double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;
            // Delta >= 0 means intersection exists.
            if ( (eq_delta > 0) || (fabs(eq_delta) < DBL_EPSILON) ) {
                active_edges.emplace_back(i);
            }
        }
    }

    return active_edges;
}

set_i
find_discontinuities(const arr_np_c_in& calc_line,
                     const arr_np_c_in& calc_corners,
                     double x_p,
                     double y_p,
                     double z_p,
                     double dt,
                     double c0)
{
    auto* line_data = calc_line.data();
    auto* corner_data = calc_corners.data();

    int N_dc = calc_corners.shape(0) + calc_line.shape(0);
    set_i dc_step;
    // calculate all of the step (t / dt) for discontinuties points
    // for corners
    for ( int i = 0; i < calc_corners.shape(0); ++i ) {
        double length = std::sqrt(
          (corner_data[i * 3 + 0] - x_p) * (corner_data[i * 3 + 0] - x_p)
          + (corner_data[i * 3 + 1] - y_p) * (corner_data[i * 3 + 1] - y_p)
          + z_p * z_p);

        double tmp_t = length / c0;
        dc_step.insert(static_cast<int>(std::ceil(tmp_t / dt)));
    }
    // for lines
    for ( int i = 0; i < calc_line.shape(0); ++i ) {
        // the minimum distance from scan point to the line.
        double tmp_t = 0;
        double x_dc = 0.0;
        double y_dc = 0.0;

        if ( static_cast<int>(line_data[i * 3 + 0] == 1) ) {
            // k = inf
            x_dc = line_data[i * 3 + 2];
            y_dc = y_p;
        } else {
            double line_k = line_data[i * 3 + 1];
            double line_b = line_data[i * 3 + 2];

            x_dc = (line_k * y_p + x_p - line_k * line_b) / (line_k * line_k + 1);
            y_dc = line_k * x_dc + line_b;
        }

        double length = std::sqrt(
          (x_dc - x_p) * (x_dc - x_p)
          + (y_dc - y_p) * (y_dc - y_p)
          + z_p * z_p);
        tmp_t = length / c0;
        dc_step.insert(static_cast<int>(std::ceil(tmp_t / dt)));
    }

    return dc_step;
}

void
circle_line_intersection(double x_p, double y_p, const arr_np_c_in& calc_line, int idx, double pow_r, double* result)
{
    auto* line_data = calc_line.data();
    // calculate intersection (convert to angle) points between circle and line.
    if ( static_cast<int>(line_data[idx * 3 + 0] == 1) ) {
        // k = inf
        double x_is = line_data[idx * 3 + 2];
        // calculate y_is by quadratic equation.
        double eq_a = 1;
        double eq_b = -2 * y_p;
        double eq_c = y_p * y_p + (x_is - x_p) * (x_is - x_p) - pow_r;
        double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;

        double y_is_1 = (-eq_b + std::sqrt(eq_delta)) / (2 * eq_a);
        double y_is_2 = (-eq_b - std::sqrt(eq_delta)) / (2 * eq_a);
        double angle_1 = convert_tan2angle(std::atan2(y_is_1 - y_p, x_is - x_p));
        double angle_2 = convert_tan2angle(std::atan2(y_is_2 - y_p, x_is - x_p));

        result[0] = angle_1;
        result[1] = angle_2;
    } else {
        double line_k = line_data[idx * 3 + 1];
        double line_b = line_data[idx * 3 + 2];
        // calculate x_is by quadratic equation.
        double eq_a = line_k * line_k + 1;
        double eq_b = 2 * line_k * line_b - 2 * x_p - 2 * line_k * y_p;
        double eq_c = y_p * y_p + x_p * x_p + line_b * line_b - 2 * line_b * y_p - pow_r;
        double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;

        double x_is_1 = (-eq_b + std::sqrt(eq_delta)) / (2 * eq_a);
        double x_is_2 = (-eq_b - std::sqrt(eq_delta)) / (2 * eq_a);
        double y_is_1 = line_k * x_is_1 + line_b;
        double y_is_2 = line_k * x_is_2 + line_b;
        double angle_1 = convert_tan2angle(std::atan2(y_is_1 - y_p, x_is_1 - x_p));
        double angle_2 = convert_tan2angle(std::atan2(y_is_2 - y_p, x_is_2 - x_p));

        result[0] = angle_1;
        result[1] = angle_2;
    }
}

double
calc_angle_difference(double angle_now,
                      double angle_next,
                      double x_p,
                      double y_p,
                      double pow_r,
                      const int* inside_direction,
                      const arr_np_c_in& boundary_line)
{
    // calculate angle difference between two angles.
    double angle_mid = (angle_now + angle_next) / 2;
    double point_x_mid = x_p + std::cos(angle_mid) * std::sqrt(pow_r);
    double point_y_mid = y_p + std::sin(angle_mid) * std::sqrt(pow_r);
    double tmp_angle_difference = 0;

    if ( point_inside_polygon(inside_direction, boundary_line, point_x_mid, point_y_mid) ) {
        tmp_angle_difference = (angle_next - angle_now);
    }

    return tmp_angle_difference;
}

arr_np_c_out
calc_polygon(arr_np_c_in boundary_line,  // line for Ax + By + C = 0, store by A、B、C.
             arr_np_c_in calc_line,
             arr_np_c_in boundary_corners,
             arr_np_c_in calc_corners,
             arr_np_c_in scan_position,
             double dt,
             double c0,
             int apo_switch)
{
    vec_set_i dc_step;
    int max_step = 0;
    // loop to find max_step for pre-allocation.
    for ( int i = 0; i < scan_position.shape(0); ++i ) {
        double x_p = scan_position(i, 0);
        double y_p = scan_position(i, 1);
        double z_p = scan_position(i, 2);

        // find all of the discontinuity points' time.
        set_i step_it = find_discontinuities(calc_line, calc_corners, x_p, y_p, z_p, dt, c0);
        if ( *step_it.begin() == 0 ) {
            step_it.erase(step_it.begin());
        }
        // calculate max step.
        if ( *step_it.rbegin() > max_step ) {
            max_step = *step_it.rbegin();
        }

        dc_step.emplace_back(step_it);
    }
    max_step += 5;  // extra 5 steps.
    // pre-allocate data. Cols is scan positions, rows is time steps.
    double* SIR_data = new double[(max_step + 1) * scan_position.shape(0)]();

    // injudge which direction about the lines is inside the polygon.
    int* inside_direction = new int[boundary_line.shape(0)]();
    inside_demarcate(boundary_line, boundary_corners, inside_direction);

    // loop all of the scan positions
    for ( int i = 0; i < scan_position.shape(0); ++i ) {
        double x_p = scan_position(i, 0);
        double y_p = scan_position(i, 1);
        double z_p = scan_position(i, 2);

        set_i dc_step_this = dc_step[i];

        // calculate SIR.
        set_i::iterator step_it = dc_step_this.begin();  // index for calculate intersection steps.
        vec_i active_edges;                              // store active edges index.
        bool full_arc_flag = false;                      // flag for full arc case.
        // apodization flag.
        int apo_start = 0;
        int apo_end = 0;
        for ( int j = *dc_step_this.begin(); j <= max_step; ++j ) {
            double pow_r = (c0 * (j * dt)) * (c0 * (j * dt)) - z_p * z_p;
            // loop lines for calculate intersection points and angles.
            set_d angles;  // Angles (band line) relative to x axis (origin is scan point).

            // test plus if and else.
            if ( j == *step_it ) {
                active_edges = update_active_edges(calc_line, x_p, y_p, z_p, dt, j, c0);
                step_it++;
                full_arc_flag = false;
            }
            // loop active edges only.
            double* intersections_angle = new double[2];
            for ( auto idx : active_edges ) {
                circle_line_intersection(x_p, y_p, calc_line, idx, pow_r, intersections_angle);
                angles.insert(intersections_angle[0]);
                angles.insert(intersections_angle[1]);
            }
            delete[] intersections_angle;

            // calculate angle difference.
            double angle_diff_all = 0.0;
            if ( full_arc_flag == true ) {
                // full arc case.
                angle_diff_all = 2 * M_PI;
            } else {
                set_d::iterator angle_it = angles.begin();
                for ( angle_it++; angle_it != angles.end(); ++angle_it ) {
                    double angle_diff_it = calc_angle_difference(*std::prev(angle_it, 1),
                                                                 *angle_it,
                                                                 x_p,
                                                                 y_p,
                                                                 pow_r,
                                                                 inside_direction,
                                                                 boundary_line);
                    angle_diff_all += angle_diff_it;
                }
                // calcualte last angle.
                double angle_diff_last = calc_angle_difference(*angles.rbegin(),
                                                               *angles.begin() + 2 * M_PI,
                                                               x_p,
                                                               y_p,
                                                               pow_r,
                                                               inside_direction,
                                                               boundary_line);
                angle_diff_all += angle_diff_last;
            }
            // change full_arc_flag.
            if ( fabs(angle_diff_all - 2 * M_PI) < DBL_EPSILON ) {
                full_arc_flag = true;
            }
            // write data in this step.
            double SIR_data_value = angle_diff_all / (2 * M_PI) * c0;
            // apodization process.
            if ( SIR_data_value && (apo_start == 0) ) {
                apo_start = j;
            }
            if ( (SIR_data_value < DBL_EPSILON) && (apo_start != 0) && (apo_end == 0) ) {
                apo_end = j;
            }
            // write SIR data value.
            SIR_data[i + j * scan_position.shape(0)] = SIR_data_value;
        }
        // apodization process.
        int apo_length = apo_end - apo_start;  // include start and not include end.
        double* apo_vec = new double[apo_length];
        apo_calc(apo_length, apo_switch, apo_vec);

        for ( int j = apo_start; j < apo_end; ++j ) {
            SIR_data[i + j * scan_position.shape(0)] *= apo_vec[j - apo_start];
        }

        delete[] apo_vec;
    }
    delete[] inside_direction;

    nb::capsule owner(SIR_data,
                      [](void* p) noexcept {
                          delete[] (double*)p;
                      });
    std::initializer_list<size_t> shape = {(size_t)(max_step + 1),
                                           (size_t)scan_position.shape(0)};

    return arr_np_c_out(SIR_data, shape, owner);
}

NB_MODULE(SIR_calc, m)
{
    m.doc() = "Calculate Spatial Impulse Response (SIR) for polygon.";

    m.def("calc_polygon", &calc_polygon, "Calculate Spatial Impulse Response (SIR) for polygon.");
}