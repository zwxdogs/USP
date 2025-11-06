// For calculate apodization weights.
#include <algorithm>
#pragma once

void
apo_calc(int length, int apo_swtich, double* apo_vec)
{
    switch ( apo_swtich ) {
        case 0:
            std::fill(apo_vec, apo_vec + length, 1.0);
            break;

        case 1:
            // hanning
            for ( int i = 0; i < length; ++i ) {
                apo_vec[i] = 0.5 - 0.5 * std::cos(2 * M_PI * i / (length - 1));
            }
            break;

        case 2:
            // hamming
            for ( int i = 0; i < length; ++i ) {
                apo_vec[i] = 0.53836 - 0.46164 * std::cos(2 * M_PI * i / (length - 1));
            }
            break;

        case 3:
            // Gaussian
            for ( int i = 0; i < length; ++i ) {
                double sigma = 0.4;
                double x = (static_cast<double>(i) - (length - 1) / 2.0) / (sigma * (length - 1) / 2.0);
                apo_vec[i] = std::exp(-0.5 * x * x);
            }
            break;

        default:
            std::fill(apo_vec, apo_vec + length, 1.0);
            break;
    }
};