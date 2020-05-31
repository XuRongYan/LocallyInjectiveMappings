//
// Created by 徐溶延 on 2020/5/30.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_ARAPENERGY_H
#define LOCALLYINJECTIVEMAPPINGS_ARAPENERGY_H

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SurfaceMesh/SurfaceMesh.h>

/**
 * 该类用于表示ARAP能量，不进行ARAP的全局求解，该负责类提供梯度算子、局部旋转矩阵R_{i}的求解
 * @tparam Scalar float or double
 */
template<class Scalar>
class ARAPEnergy {
    /**
     * Whether the vertex matrix is transposed
     */
    enum MODE {
        NOT_TRANSPOSE,
        TRANSPOSE,
    };

public:
    ARAPEnergy();

    ARAPEnergy(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &V, const Eigen::Matrix3Xi &F, MODE mode = TRANSPOSE);

    ARAPEnergy(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &v, const Eigen::Matrix3Xi &f,
               const std::vector<Eigen::Matrix<Scalar, 2, 3>> &idealElems, MODE mode = TRANSPOSE);

    int precompute();

    void computeLocalR(std::vector<Eigen::Matrix<Scalar, 2, 2>> &R);

    const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &getV() const;

    void setV(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &v);

    const Eigen::Matrix3Xi &getF() const;

    void setF(const Eigen::Matrix3Xi &f);

    const Eigen::SparseMatrix<float> &getG() const;

    const std::vector<Eigen::Matrix<Scalar, 3, 3>> &getR() const;

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
                         Eigen::Matrix<Scalar, 2, 2> &R);

    /**
     * flatten triangle to plane to get ideal elements
     * @param p1
     * @param p2
     * @param p3
     */
    Eigen::Matrix<Scalar, 2, 3> flattenTriangle(const Eigen::Matrix<Scalar, 3, 1> &p1,
                         const Eigen::Matrix<Scalar, 3, 1> &p2,
                         const Eigen::Matrix<Scalar, 3, 1> &p3);

    void computeIdeals();

    void computeInvDeltas();

    void computeAreas();

    void computeGradientOperator();

};


#endif //LOCALLYINJECTIVEMAPPINGS_ARAPENERGY_H
