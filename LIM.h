//
// Created by 徐溶延 on 2020/5/30.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_LIM_H
#define LOCALLYINJECTIVEMAPPINGS_LIM_H

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SurfaceMesh/SurfaceMesh.h>

#include "ARAPEnergy.h"

// TODO 将该类写为Builder模式，更方便参数设置
/**
 * This is the class of Locally Injective Mapping(LIM)
 */
class LIM {
    typedef Surface_Mesh::Scalar Scalar;
public:
    LIM();

    LIM(const Surface_Mesh::SurfaceMesh &mesh);

    LIM(const Surface_Mesh::SurfaceMesh &mesh,
        const std::vector<std::pair<int, Eigen::Matrix3f>> &posConstrains);

    void rebuild(const Surface_Mesh::SurfaceMesh &mesh);

    void rebuild(const Surface_Mesh::SurfaceMesh &mesh,
                 const std::vector<std::pair<int, Eigen::Matrix3f>> &posConstrains);

    void solve();

    const std::vector<std::pair<int, Eigen::Matrix3f>> &getPosConstrains() const;

    void setPosConstrains(const std::vector<std::pair<int, Eigen::Matrix3f>> &posConstrains);

    const ARAPEnergy<Surface_Mesh::Scalar> &getArapEnergy() const;

    void setArapEnergy(const ARAPEnergy<Surface_Mesh::Scalar> &arapEnergy);

    const Surface_Mesh::SurfaceMesh &getMesh() const;

    void setMesh(const Surface_Mesh::SurfaceMesh &mesh);

private:
    Scalar s_j_, epsilon_, mu_;
    Surface_Mesh::SurfaceMesh mesh_;
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> V_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vt_calc_, d_i_, t_i_; //用于计算的v向量
    Eigen::Matrix3Xi F_;
    std::vector<std::pair<int, Eigen::Matrix3f>> pos_constrains_;
    ARAPEnergy<Surface_Mesh::Scalar> arap_energy_;
    Eigen::SparseMatrix<Scalar> H_rigid_, H_pos_, H_barries, J_rigid_, J_pos_, J_barries; //Hessian and Jacobian
    Eigen::SparseMatrix<Scalar> M_; //面积公式的二次型矩阵
    std::vector<Scalar> weight1, weight2; //hessian 和 jacobian的权值矩阵

    void precompute();

    void computePosConstrainDi();

    void computeTi(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &C);

    void computePosJacobian();

    void computePosHessian();

    void computeRigidJacobian();

    void computeRigidHessian();

    void computeBarriesJacobian();

    void computeBarriesHessian();

    void computeEpsilon();

    void computeM(const std::vector<Scalar> &weights);

    void computeW1();

    void computeW2();

    /**
     * compute barries function
     * @param V
     * @return
     */
    Scalar barries();

    Scalar phi(Scalar x);

    Scalar g(Scalar x);

    Scalar gPrime(Scalar x);

    Scalar gPrime2(Scalar x);

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> lambdaPrime() const;

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> lambdaPrime2() const;

    Scalar c(const Eigen::Matrix<Scalar, 3, 1> &p1,
             const Eigen::Matrix<Scalar, 3, 1> &p2,
             const Eigen::Matrix<Scalar, 3, 1> &p3);

    void reset();

};


#endif //LOCALLYINJECTIVEMAPPINGS_LIM_H
