//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <vector>
#include <algorithm>

// a simple structure for the sparse mass matrix

struct SparseMatrix
{
    // a triplet for the sparse mass matrix
    struct triplet_t {
        int i,j;
        double val;
        triplet_t() = default;
        triplet_t(int r, int c, double v)
          : i(r), j(c), val(v) {}

        inline bool operator < (const triplet_t& b) const {
            if (i < b.i) return true;
            if (i > b.i) return false;
            if (j < b.j) return true;
            if (j > b.j) return false;
            return (val < b.val);
        }
        inline bool accum(triplet_t& b) {
            if (b.i == i && b.j == j) {
                val += b.val;
                b.val = 0.0;
                b.i = b.j = -1;
                return true;
            }
            return false;
        }
    };
    // all triplets
    std::vector<triplet_t> triplets;

    struct row_t {
        size_t pos = 0;
        size_t count = 0;
    };
    // row access data
    std::vector<row_t> rows;

    // append a new triplet
    inline void append(int i, int j, double value) {
        triplets.push_back({ i, j, value });
    }

    // get a value, returns 0.0 if not found
    inline double get(int i, int j) const {
        if (i < 0 || i >= static_cast<int>(rows.size()))
            return 0.0;

        const auto& row_data = rows[static_cast<size_t>(i)];
        if (row_data.count > 0) {
            for (size_t c = 0; c < row_data.count; ++c) {
                const auto& it = triplets[c + row_data.pos];
                if (it.j == j)
                    return it.val; // found
                if (it.j > j)
                    break; // no need to search further
            }
        }
        return 0.0;
    }

    // call this after all append calls
    inline void finish() {
        if (triplets.size() == 0)
            return;
        // sort triplets first by row, then by column
        std::sort(triplets.begin(), triplets.end());
        size_t i_current = 0;
        size_t count = 1;
        for (size_t i = 1; i < triplets.size(); ++i) {
            auto& it = triplets[i];
            if (!triplets[i_current].accum(it)) {
                i_current = i;
                ++count;
            }
        }
        // exclude duplicates
        std::vector<triplet_t> aux;
        aux.swap(triplets);
        triplets.resize(count);
        count = 0;
        for (const auto& it : aux)
            if (it.i > -1)
                triplets[count++] = it;
        // make row structure for fast access
        size_t max_row = static_cast<size_t>(triplets.back().i);
        rows.resize(max_row + 1);
        for (size_t i = 0; i < triplets.size(); ++i) {
            const auto& it = triplets[i];
            auto& row_data = rows[static_cast<size_t>(it.i)];
            if (row_data.count == 0)
                row_data.pos = i; // on first access, save the position
            ++(row_data.count); // increase counter
        }
    }
};

