//
// Created by 徐溶延 on 2020/6/1.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_EIGENUTILS_H
#define LOCALLYINJECTIVEMAPPINGS_EIGENUTILS_H

#include <Eigen/Dense>

namespace xry_mesh {

    /**
     * 将vt矩阵转换为v矩阵
     * @tparam Scalar
     * @param vt = {{x1, y1, z1};{x2, y2, z2}; ...  ;{xn, yn, zn}}
     * @return v = {x1, y1, z1, x2, y2, z2, ..., xn, yn, zn}
     */
    template<class Scalar>
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    vt2v(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &vt);

    /**
     * 将v矩阵转换为vt矩阵
     * @tparam Scalar
     * @param v = {x1, y1, z1, x2, y2, z2, ..., xn, yn, zn}
     * @return vt = {{x1, y1, z1};{x2, y2, z2}; ...  ;{xn, yn, zn}}
     */
    template<class Scalar>
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    v2vt(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &v, int dim);


} // namespace xry_mesh

#endif //LOCALLYINJECTIVEMAPPINGS_EIGENUTILS_H
