//
// Created by knirschl on 22.05.23.
//

#ifndef BA_NJSIMPLE_H
#define BA_NJSIMPLE_H

#include "tree.h"
#include "../misc/common_types.h"

template<typename T>
vector_t<T> meanDistance(matrix_t<T> const &distMatrix) {
    const size_t n = distMatrix.size();
    vector_t<T> r;
    r.reserve(n);
    for (size_t i{0}; i < n; ++i) {
        T sumK{0};
        for (size_t k{0}; k < n; ++k) {
            sumK += distMatrix[i][k];
        }
        r[i] = (1.0 / (n - 2)) * sumK; // 1.0 for double division
    }

    return r;
}

template<typename T>
matrix_t<T> interMatrix(matrix_t<T> const &distMatrix, vector_t<T> const &r) {
    const size_t n = distMatrix.size();
    matrix_t<T> iM;
    iM.reserve(n);
    for (size_t i{0}; i < n; ++i) {
        vector_t<T> row(n);
        for (size_t j{0}; j < n; ++j) {
            row[j] = distMatrix[i][j] - (r[i] + r[j]);
        }
        iM.push_back(row);
    }

    return iM;
}

template<typename T>
std::pair<size_t, size_t> findMin(matrix_t<T> const &matrix) {
    T min = std::numeric_limits<T>::infinity();
    const size_t n = matrix.size();
    std::pair<size_t, size_t> minDistPos;
    for (size_t i{0}; i < n; ++i) {
        for (size_t j{0}; j < n; ++j) {
            if (i == j) {
                continue;
            }
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                minDistPos.first = i;
                minDistPos.second = j;
            }
        }
    }

    return minDistPos;
}

template<typename T>
matrix_t<T>
calculateNewMatrix(matrix_t<T> const &distMatrix, std::pair<size_t, size_t> const &minDistPos) {
    const size_t n = distMatrix.size();
    matrix_t<T> newMatrix{};
    newMatrix.reserve(n - 1);
    // pre-initialize to prevent creating one too many
    for (size_t i{0}; i < n - 1; ++i) {
        newMatrix.emplace_back(n - 1);
    }
    // copy old values of rows/columns that won't be changed
    size_t xOff{0};
    for (size_t i{0}; i < n; ++i) {
        // skip row if one of the merged ones
        if (i == minDistPos.first || i == minDistPos.second) {
            xOff++;
        } else {
            size_t yOff{0};
            for (size_t j{0}; j < n; ++j) {
                // skip column if one of the merged ones
                if (j == minDistPos.first || j == minDistPos.second) {
                    yOff++;
                } else {
                    // index at newMatrix depends on how many rows/columns were already skipped
                    newMatrix[i - xOff][j - yOff] = distMatrix[i][j];
                }
            }
        }
    }

    // calculate new distances from merged clusters to everything else, put in last row/column
    size_t off{0};
    for (size_t k{0}; k < n - 2; ++k) {
        // calculate offset (how much needs to be skipped) via black magic
        while (k + off == minDistPos.first || k + off == minDistPos.second) {
            off++;
        }
        /*
         * d_u,k = (d_i,k + d_j,k - d_i,j) / 2
         * u = (i,j) the merged clusters, k is the current cluster
         */
        newMatrix[k][n - 2] =
                (distMatrix[minDistPos.first][k + off] + distMatrix[minDistPos.second][k + off] -
                 distMatrix[minDistPos.first][minDistPos.second]) / 2.0;
        newMatrix[n - 2][k] =
                (distMatrix[minDistPos.first][k + off] + distMatrix[minDistPos.second][k + off] -
                 distMatrix[minDistPos.first][minDistPos.second]) / 2.0;
    }
    newMatrix[n - 2][n - 2] = 0; // set new diagonal entry to 0

    return newMatrix;
}

// const& to be sdure if bigger datatypes than simple numbers are used
template<typename T>
std::pair<T, T> calculateBranchLengths(T const &dist, T const &r1, T const &r2) {
    return std::pair<T, T>{(dist + r1 - r2) / 2.0, dist - ((dist + r1 - r2) / 2.0)};
}

vector_t<std::shared_ptr<Tree>> copyUnjoinedTrees(vector_t<std::shared_ptr<Tree>> const &trees,
                                                  std::pair<size_t, size_t> const &minDistPos) {
    vector_t<std::shared_ptr<Tree>> copy;
    const size_t n = trees.size();
    for (size_t i{0}; i < n; ++i) {
        // skip and add offset if this tree is one of the two that will be merged
        if (i == minDistPos.first || i == minDistPos.second) {
            continue;
        }
        copy.emplace_back(trees[i]);
    }

    return copy;
}

template<typename T>
int neighborJoining(matrix_t<T> const &distMatrix, std::shared_ptr<Tree> &tree,
                    std::vector<int> active_nodes) {
    const size_t n = distMatrix.size();
    /*static_assert([n, distMatrix]() {
        bool quadratic = true;
        for (auto row : distMatrix) {
            quadratic &= row.size() == n;
        }
        return quadratic;
    }, "distMatrix needs to be quadratic");
    static_assert(n == trees.size(), "Same amount of leafs as rows in matrix needed");*/
    if (n > 2) {
        // calculations
        vector_t<T> r = meanDistance(distMatrix);
        std::pair<size_t, size_t> minDistPos = findMin(interMatrix(distMatrix, r));
        matrix_t<T> newDistMatrix = calculateNewMatrix(distMatrix, minDistPos);

        //tree building
        std::pair<T, T> branchLengths = calculateBranchLengths(
                distMatrix[minDistPos.first][minDistPos.second], r[minDistPos.first],
                r[minDistPos.second]);
        // add joined trees (neighbors)
        // delete minDistPos and add new node idx to activeNodes
        active_nodes.emplace_back(
                tree->make_node(active_nodes[minDistPos.first], branchLengths.first,
                                active_nodes[minDistPos.second],
                                branchLengths.second));
        active_nodes.erase(active_nodes.begin() + minDistPos.first);
        active_nodes.erase(active_nodes.begin() + minDistPos.second - 1);

        return neighborJoining(newDistMatrix, tree, active_nodes);
    } else if (n == 2) {
        // base case
        // connecting branch defined by two branches of same length to parent -> "rooted" tree...
        tree->root = tree->make_node(active_nodes[0], distMatrix[0][1], active_nodes[1],
                                     distMatrix[0][1]);
        return tree->root;
    } else {
        // edge case, only one or none clusters
        return tree->tree[active_nodes[0]].idx;
    }
}

#endif //BA_NJSIMPLE_H
