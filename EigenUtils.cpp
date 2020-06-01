//
// Created by 徐溶延 on 2020/6/1.
//

#include "EigenUtils.h"

namespace xry_mesh {

    template<class Scalar>
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    vt2v(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &vt) {
        const int dim = vt.cols();
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> v(dim * vt.rows());
        for (size_t i = 0; i < vt.rows(); i++) {
            for (size_t j = 0; j < dim; j++) {
                v[3 * i + j] = vt(i, j);
            }
        }
        return v;
    }

    template<class Scalar>
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    v2vt(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &v, int dim) {
        assert(v.rows() % dim == 0);
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> vt(v.rows() / dim, dim);
        for (size_t i = 0; i < v.rows(); i += dim) {
            for (size_t j = 0; j < dim; j++) {
                vt(i / dim, j) = v[i + j];
            }
        }
        return vt;
    }

} // namespace