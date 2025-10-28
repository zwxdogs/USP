// Calculate Spatial Impulse Response (SIR).
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <vector>
#include <nanobind/stl/vector.h>
#include <algorithm>
#include <math.h>
#define _USE_MATH_DEFINES

namespace nb = nanobind;
using arr_np_c_in = nb::ndarray<double, nb::numpy, nb::c_contig, nb::shape<-1, 3>>;
using arr_np_c_out = nb::ndarray<double, nb::numpy, nb::c_contig>;
using vec_vec_d = std::vector<std::vector<double>>;
using vec_d = std::vector<double>;
using vec_i = std::vector<int>;

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
inside_demarcate(arr_np_c_in probe_line, arr_np_c_in probe_corners)
{
    // injudge which direction about the lines is inside the polygon.

    // calculate centroid of the polygon.
    double centroid_x = 0.0;
    double centroid_y = 0.0;
    for ( int i = 0; i < probe_corners.shape(0); ++i ) {
        centroid_x += probe_corners(i, 0);
        centroid_y += probe_corners(i, 1);
    }
    centroid_x /= probe_corners.shape(0);
    centroid_y /= probe_corners.shape(0);

    vec_i inside_direction;
    // injudge and store.
    for ( int i = 0; i < probe_line.shape(0); ++i ) {
        double line_a = probe_line(i, 0);
        double line_b = probe_line(i, 1);
        double line_c = probe_line(i, 2);

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
point_inside_polygon(vec_i inside_direction,
                     arr_np_c_in probe_line,
                     double x_point,
                     double y_point)
{
    // injudge whether the point is inside the polygon.
    for ( int i = 0; i < probe_line.shape(0); ++i ) {
        double line_a = probe_line(i, 0);
        double line_b = probe_line(i, 1);
        double line_c = probe_line(i, 2);

        double injudge_value = line_a * x_point + line_b * y_point + line_c;
        if ( injudge_value * inside_direction[i] <= 0 ) {
            return false;
        }
    }
    return true;
}

vec_vec_d
calc_polygon(arr_np_c_in probe_line,
             arr_np_c_in probe_corners,
             arr_np_c_in scan_position,
             double dt,
             double c0)
{
    // injudge which direction about the lines is inside the polygon.
    vec_i inside_direction = inside_demarcate(probe_line, probe_corners);

    vec_vec_d SIR_result;
    // loop all of the scan positions
    for ( int i = 0; i < scan_position.shape(0); ++i ) {
        double x_p = scan_position(i, 0);
        double y_p = scan_position(i, 1);
        double z_p = scan_position(i, 2);

        // calculate all of the t (step) for discontinuties points
        int N_dc = probe_corners.shape(0) + probe_line.shape(0);
        vec_d dc_t(N_dc, 0.0);
        for ( int j = 0; j < probe_corners.shape(0); ++j ) {
            double length = std::sqrt(
              std::pow(probe_corners(j, 0) - x_p, 2)
              + std::pow(probe_corners(j, 1) - y_p, 2)
              + std::pow(z_p, 2));

            double tmp_t = length / c0;
            dc_t[j] = tmp_t;
        }
        for ( int j = 0; j < probe_line.shape(0); ++j ) {
            // the minimum distance from scan point to the line.
            if ( probe_line(j, 1) == 0 ) {
                // k = inf
                double x_dc = -probe_line(j, 2) / probe_line(j, 0);
                double y_dc = y_p;
                double length = std::sqrt(
                  std::pow(x_dc - x_p, 2)
                  + std::pow(y_dc - y_p, 2)
                  + std::pow(z_p, 2));

                double tmp_t = length / c0;
                dc_t[probe_corners.shape(0) + j] = tmp_t;
            } else {
                double line_k = -probe_line(j, 0) / probe_line(j, 1);
                double line_b = -probe_line(j, 2) / probe_line(j, 1);
                double x_dc = (line_k * y_p + x_p - line_k * line_b) / (line_k * line_k + 1);
                double y_dc = line_k * x_dc + line_b;
                double length = std::sqrt(
                  std::pow(x_dc - x_p, 2)
                  + std::pow(y_dc - y_p, 2)
                  + std::pow(z_p, 2));

                double tmp_t = length / c0;
                dc_t[probe_corners.shape(0) + j] = tmp_t;
            }
        }
        vec_i dc_step(N_dc, 0);
        for ( int j = 0; j < N_dc; ++j ) {
            dc_step[j] = static_cast<int>(std::ceil(dc_t[j] / dt));
        }
        // sort dc_t and remove repetition.
        std::sort(dc_step.begin(), dc_step.end());
        dc_step.erase(std::unique(dc_step.begin(), dc_step.end()), dc_step.end());
        int max_step = dc_step.back();

        // ----------------------------------------------------------------------------------------------------- //
        // calculate SIR.
        vec_d SIR_data(max_step, 0.0);

        for ( int j = dc_step.front(); j <= max_step; ++j ) {
            double pow_r = std::pow(c0 * (j * dt), 2) - std::pow(z_p, 2);
            // loop lines for calculate intersection points and angles.
            vec_d angles;  // Angles relative to x axis (origin is scan point).
            for ( int k = 0; k < probe_line.shape(0); ++k ) {
                if ( probe_line(k, 1) == 0 ) {
                    // k = inf
                    double x_is = -probe_line(k, 2) / probe_line(k, 0);
                    // calculate y_is by quadratic equation.
                    double eq_a = 1;
                    double eq_b = -2 * y_p;
                    double eq_c = y_p * y_p + std::pow(x_is - x_p, 2) - pow_r;
                    double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;

                    if ( eq_delta >= 0 ) {
                        double y_is_1 = (-eq_b + std::sqrt(eq_delta)) / (2 * eq_a);
                        double y_is_2 = (-eq_b - std::sqrt(eq_delta)) / (2 * eq_a);
                        double angle1 = convert_angle(std::atan2(y_is_1 - y_p, x_is - x_p), y_is_1 - y_p);
                        double angle2 = convert_angle(std::atan2(y_is_2 - y_p, x_is - x_p), y_is_2 - y_p);
                        angles.emplace_back(angle1);
                        angles.emplace_back(angle2);
                    }
                } else {
                    double line_k = -probe_line(k, 0) / probe_line(k, 1);
                    double line_b = -probe_line(k, 2) / probe_line(k, 1);
                    // calculate x_is by quadratic equation.
                    double eq_a = line_k * line_k + 1;
                    double eq_b = 2 * line_k * line_b - 2 * x_p - 2 * line_k * y_p;
                    double eq_c = y_p * y_p + x_p * x_p + line_b * line_b - 2 * line_b * y_p - pow_r;
                    double eq_delta = eq_b * eq_b - 4 * eq_a * eq_c;

                    if ( eq_delta >= 0 ) {
                        double x_is_1 = (-eq_b + std::sqrt(eq_delta)) / (2 * eq_a);
                        double x_is_2 = (-eq_b - std::sqrt(eq_delta)) / (2 * eq_a);
                        double y_is_1 = line_k * x_is_1 + line_b;
                        double y_is_2 = line_k * x_is_2 + line_b;
                        double angle1 = convert_angle(std::atan2(y_is_1 - y_p, x_is_1 - x_p), y_is_1 - y_p);
                        double angle2 = convert_angle(std::atan2(y_is_2 - y_p, x_is_2 - x_p), y_is_2 - y_p);
                        angles.emplace_back(angle1);
                        angles.emplace_back(angle2);
                    }
                }
            }
            std::sort(angles.begin(), angles.end());
            angles.erase(std::unique(angles.begin(), angles.end()), angles.end());

            // find which point in the aperture and calculate sum of the angle difference.
            double angle_difference = 0;
            for ( int k = 1; k < angles.size(); ++k ) {
                // loop angles.
                double angle_now = angles[k - 1];
                double angle_next = angles[k];
                double angle_mid = (angle_now + angle_next) / 2;

                double point_x_mid = x_p + std::cos(angle_mid) * std::sqrt(pow_r);
                double point_y_mid = y_p + std::sin(angle_mid) * std::sqrt(pow_r);

                if ( point_inside_polygon(inside_direction, probe_line, point_x_mid, point_y_mid) ) {
                    angle_difference += (angle_next - angle_now);
                }
            }
            // write data in this step.
            SIR_data[j - 1] = angle_difference * c0 / (2 * M_PI);
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