//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//

#include <cmath>

class PredictorControl {
    public:
    PredictorControl( double base_step,
                      int target_count,
                      double min_size,
                      double max_size,
                      double exponent)
    // : base_step(base_step)
    : target_count(target_count)
    , min_size(min_size)
    , max_size(max_size)
    , exponent(exponent)
    , last_step(base_step)
    , next_step(base_step)
    , current_count(0)
    {

    }
    
    void reset(double step_size) {
        current_count = 0;
        last_step = step_size;
        next_step = step_size;
    }

    int update() {
        current_count++;
        return current_count;
    }

    double size() const {
        double factor = std::pow(target_count / current_count, exponent);
        double step = last_step * factor;

        if (step < min_size)
            step = min_size;
        else if (step > max_size)
            step = max_size;
        return step;
    }

    
    double last_step;
    double next_step;
    double current_count;

    private:
    double target_count;
    double min_size;
    double max_size;
    double exponent;
    // double base_step;
    //
};