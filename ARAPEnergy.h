//
// Created by 徐溶延 on 2020/5/30.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_ARAPENERGY_H
#define LOCALLYINJECTIVEMAPPINGS_ARAPENERGY_H

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SurfaceMesh/SurfaceMesh.h>

#include "SurfaceMeshUtils.h"

/**
 * 该类用于表示ARAP能量，不进行ARAP的全局求解，该负责类提供梯度算子、局部旋转矩阵R_{i}的求解
 * @tparam Scalar float or double
 */
template<typename Scalar>
class ARAPEnergy {
    /**
     * Whether the vertex matrix is transposed
     */
    enum MODE {
        NOT_TRANSPOSE,
        TRANSPOSE,
    };

public:
    ARAPEnergy() {

    }

    ARAPEnergy(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &V, const Eigen::Matrix3Xi &F, MODE mode = TRANSPOSE):V_(V),
                                                                                                                    F_(F),
                                                                                                                    mode(mode) {
    }

    ARAPEnergy(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &v, const Eigen::Matrix3Xi &f,
               const std::vector<Eigen::Matrix<Scalar, 2, 3>> &idealElems, MODE mode = TRANSPOSE) :V_(v), F_(f),
                                                                                                   ideal_elems_(idealElems),
                                                                                                   mode(mode) {
    }

    void rebuild(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &V, const Eigen::Matrix3Xi &F, MODE mode = TRANSPOSE) {
        reset();
        this->V_ = V;
        this->F_ = F;
        this->mode = mode;
    }

    void rebuild(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &v, const Eigen::Matrix3Xi &f,
                 const std::vector<Eigen::Matrix<Scalar, 2, 3>> &idealElems, MODE mode = TRANSPOSE) {
        reset();
        this->V_ = v;
        this->F_ = f;
        this->mode = mode;
        this->ideal_elems_ = idealElems;
    }

    int precompute() {
        computeAreas();
        computeIdeals();
        computeInvDeltas();
        computeGradientOperator();
        return 0;
    }

    void computeLocalR(std::vector<Eigen::Matrix<Scalar, 2, 2>> &R) {
        // TODO 似乎没有必要在Energy中保存R，R交给外部算法维护是否可行
        R_.clear();
        R_.resize(F_.cols());
        Eigen::Matrix<Scalar, 2, 2> m;
        for (size_t i = 0; i < F_.cols(); i++) {
            for (size_t j = 0; j < 2; j++) {
                m.col(j) = V_.col(F_(j + 1, i)) - V_.col(F_(0, i));
            }
            optimizeRotation(m * inv_deltas_[i], R_[i]);
        }
        R = R_;
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> getRCalc() {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> R;
        if (mode == TRANSPOSE) {
            R.resize(2 * F_.cols(), 2);
            for (size_t i = 0; i < R_.size(); i++) {
                for (size_t j = 0; j < 2; j++) {
                    for (size_t k = 0; k < 3; k++) {
                        R(2 * i + j, k) = R_[i](j, k);
                    }
                }
            }
        } else {
            R.resize(4 * F_.cols(), 1);
            for (size_t i = 0; i < R_.size(); i++) {
                for (size_t j = 0; j < 2; j++) {
                    for (size_t k = 0; k < 3; k++) {
                        R(4 * i + 2 * j + k, 0) = R_[i](j, k);
                    }
                }
            }
        }
        return R;
    }

    const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &getV() const {
        return V_;
    }

    void setV(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &v) {
        V_ = v;
    }

    const Eigen::Matrix3Xi &getF() const {
        return F_;
    }

    void setF(const Eigen::Matrix3Xi &f) {
        F_ = f;
    }

    const Eigen::SparseMatrix<float> &getG() const {
        return G_;
    }

    const std::vector<Eigen::Matrix<Scalar, 3, 3>> &getR() const {
        return R_;
    }

private:
    MODE mode;
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> V_;
    Eigen::Matrix3Xi F_;
    Eigen::SparseMatrix<float> G_;
    std::vector<Scalar> areas;
    std::vector<Eigen::Matrix<Scalar, 2, 2>> R_;
    std::vector<Eigen::Matrix<Scalar, 2, 3>> ideal_elems_;
    std::vector<Eigen::Matrix<Scalar, 2, 2>> inv_deltas_;
    Eigen::SparseMatrix<Scalar> L_;

    /**
     * locally compute rotation matrix of a triangle
     * @param J PQ^{-1}
     * @param R rotation matrix
     * @return 1 if success, 0 if failed
     */
    int optimizeRotation(const Eigen::Matrix<Scalar, 2, 2> &J,
                         Eigen::Matrix<Scalar, 2, 2> &R) {
        Eigen::JacobiSVD<Eigen::Matrix<Scalar, 2, 2>> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
        R = svd.matrixU() * svd.matrixV().transpose();
        if (fabs(svd.singularValues()(2)) < 1e-8) {
            R.setIdentity();
            return 1;
        }
        if (R.determinant() < 0) {
            Eigen::MatrixXd svdV = svd.matrixV();
            svdV.col(2) = -svdV.col(2);
            R = svd.matrixU() * svdV.transpose();
        }
        return 0;
    }

    /**
     * flatten triangle to plane to get ideal elements
     * @param p1
     * @param p2
     * @param p3
     */
    Eigen::Matrix<Scalar, 2, 3> flattenTriangle(const Eigen::Matrix<Scalar, 3, 1> &p1,
                         const Eigen::Matrix<Scalar, 3, 1> &p2,
                         const Eigen::Matrix<Scalar, 3, 1> &p3) {
        const Scalar aa = (p2 - p1).squaredNorm();
        const Scalar bb = (p3 - p2).squaredNorm();
        const Scalar cc = (p3 - p1).squaredNorm();

        const Scalar a = std::sqrt(aa);
        const Scalar c = std::sqrt(cc);
        assert(a != 0 && c != 0);
        assert((aa + cc - bb) / (2 * a * c) >= -1 && (aa + cc - bb) / (2 * a * c) <= 1);
        const Scalar angle = std::acos((aa + cc - bb) / (2 * a * c));

        Eigen::Matrix<Scalar, 2, 3> p2d;
        p2d << 0, a, c * std::cos(angle),
                0, 0, c * std::sin(angle);
        return p2d;
    }

    void computeIdeals() {
        ideal_elems_.clear();
        for (size_t i = 0; i < F_.cols(); i++) {
            Eigen::Matrix<Scalar, 2, 3> ideal = flattenTriangle(V_.col(F_(0, i)),
                                                                V_.col(F_(1, i)),
                                                                V_.col(F_(2, i)));
            ideal_elems_.emplace_back(ideal);
        }
    }

    void computeInvDeltas() {
        inv_deltas_.clear();
        assert(ideal_elems_.size() == F_.cols());
        for (size_t i = 0; i < F_.cols(); i++) {
            Eigen::Matrix<Scalar, 2, 2> m;
            Eigen::Matrix<Scalar, 2, 3> ideal = ideal_elems_[i];
            for (size_t j = 0; j < 2; j++) {
                m.col(j) = ideal.col(j + 1) - ideal.col(0);
            }
            inv_deltas_.emplace_back(m.inverse());
        }
    }

    void computeAreas() {
        areas.clear();
        for (size_t i = 0; i < F_.cols(); i++) {
            areas.emplace_back(xry_mesh::computeArea<Scalar>(V_.col(F_(0, i)),
                                                             V_.col(F_(1, i)),
                                                             V_.col(F_(2, i))));
        }
    }

    void computeGradientOperator() {
        std::vector<Eigen::Triplet<Scalar>> triplets;
        for (size_t i = 0, col = 0; i < F_.cols(); i++) {
            for (size_t j = 0; j < 2; j++, col++) {
                double sum = 0;
                for (size_t k = 0; k < 2; k++) {
                    const double val = inv_deltas_[i](k, j) * sqrt(areas[i]);
                    sum -= val;
                    if (mode == NOT_TRANSPOSE) { // if vertex matrix is v
                        for (size_t t = 0; t < 2; t++) {
                            triplets.emplace_back(2 * F_(k + 1, i), 2 * col, val);
                        }
                    } else { // if vertex matrix is v^{T}
                        triplets.emplace_back(F_(k + 1, i), col, val);
                    }
                }
                if (mode == NOT_TRANSPOSE) {
                    for (size_t t = 0; t < 2; t++) {
                        triplets.emplace_back(2 * F_(0, i), 2 * col, sum);
                    }
                } else {
                    triplets.emplace_back(F_(0, i), col, sum);
                }
            }
        }
        Eigen::SparseMatrix<Scalar> Gt;
        if (mode == TRANSPOSE) {
            Gt.resize(V_.cols(), 2 * F_.cols());
        } else {
            Gt.resize(2 * V_.cols(), 4 * F_.cols());
        }
        Gt.setFromTriplets(triplets.begin(), triplets.end());
        G_ = Gt.transpose();
    }

    void reset() {
        G_.setZero();
        areas.clear();
        R_.clear();
        ideal_elems_.clear();
        inv_deltas_.clear();
    }

};


#endif //LOCALLYINJECTIVEMAPPINGS_ARAPENERGY_H
