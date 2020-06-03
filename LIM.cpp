//
// Created by 徐溶延 on 2020/5/30.
//

#include "LIM.h"

#include <map>

#include <dbg.h>

#include "SurfaceMeshUtils.h"
#include "EigenUtils.h"


LIM::LIM() {}

LIM::LIM(const Surface_Mesh::SurfaceMesh &mesh) : mesh_(mesh) {

}

LIM::LIM(const Surface_Mesh::SurfaceMesh &mesh,
         const std::vector<std::pair<int, Eigen::Matrix<Scalar, 3, 1>>> &posConstrains)
        : mesh_(mesh), pos_constrains_(posConstrains) {

}

void LIM::rebuild(const Surface_Mesh::SurfaceMesh &mesh) {
    reset();
    this->mesh_ = mesh;
    precompute();
}

void LIM::rebuild(const Surface_Mesh::SurfaceMesh &mesh,
                  const std::vector<std::pair<int, Eigen::Matrix<Scalar, 3, 1>>> &posConstrains) {
    reset();
    this->mesh_ = mesh;
    this->pos_constrains_ = posConstrains;

}

void LIM::precompute() {
    assert(!mesh_.empty());
    s_j_ = 0.5;
    V_ = xry_mesh::getGeometryMatrix<Scalar>(mesh_);
    //取前两列
    Vt_calc_ = xry_mesh::vt2v<Scalar>(V_.transpose().block(0, 0, V_.cols(), 2));
    F_ = xry_mesh::getTopoMatrix(mesh_);
    arap_energy_.rebuild(V_, F_);
    arap_energy_.precompute();
    xry_mesh::computeAreas<Scalar>(mesh_);
    computePosConstrainDi();
    //computeM();
    computeEpsilon();
    computeW1();
    computeW2();
}

void LIM::solve() {

}

void LIM::computePosHessian() {
    Eigen::SparseMatrix<Scalar> C(Vt_calc_.size(), Vt_calc_.size());
    C.setIdentity();
    H_pos_ = 2 * C.transpose() * C;
}

void LIM::computePosJacobian() {
    Eigen::SparseMatrix<Scalar> C(Vt_calc_.size(), Vt_calc_.size());
    C.setIdentity();
    computePosConstrainDi();
    computeTi(C);
    Eigen::SparseMatrix<Scalar> Ct = C.transpose();
    auto J_pos_dense = 2 * Ct * C * Vt_calc_ - 2 * Ct * t_i_;
    J_pos_ = J_pos_dense.sparseView();
}

void LIM::computeRigidJacobian() {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> R;
    Eigen::SparseMatrix<Scalar> G;
    G = arap_energy_.getG();
    R = arap_energy_.getRCalc();
    J_rigid_ = 2 * G.transpose() * G - 2 * G.transpose() * R;

}

void LIM::computeRigidHessian() {
    Eigen::SparseMatrix<Scalar> G;
    G = arap_energy_.getG();
    H_rigid_ = G.transpose() * G;
}

void LIM::computeBarriesJacobian() {

}

void LIM::computeBarriesHessian() {

}

void LIM::computePosConstrainDi() {
    d_i_.resize(Vt_calc_.size());
    d_i_.setZero();
    for (const auto &pair : pos_constrains_) {
        for (size_t j = 0; j < 3; j++) {
            d_i_[3 * pair.first + j] = pair.second(j, 0);
        }
    }
}

void LIM::computeTi(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &C) {
    // TODO 初始化mu_
    const Scalar mumu = mu_ * mu_;
    t_i_ = C * Vt_calc_ + 1.0 / (1 + mumu) * (d_i_ - C * Vt_calc_);
}

void LIM::computeM(const std::vector<Scalar> &weights) {

}

void LIM::computeEpsilon() {
    Scalar MIN = INT_MAX;
    auto areas = mesh_.get_face_property<Scalar>("f:areas");
    for (const auto &f : mesh_.faces()) {
        auto area = areas[f];
        MIN = std::min(MIN, area);
    }
    epsilon_ = 1e-5 * MIN;
}

void LIM::computeW1() {
    auto areas = mesh_.get_face_property<Scalar>("f:areas");
    weight1.clear();
    weight1.resize(mesh_.n_faces());
    for (const auto &f : mesh_.faces()) {
        weight1[f.idx()] = -gPrime(areas[f]) / std::pow(g(areas[f]), 2);
    }
}

void LIM::computeW2() {
    auto areas = mesh_.get_face_property<Scalar>("f:areas");
    weight2.clear();
    weight2.resize(mesh_.n_faces());
    for (const auto &f : mesh_.faces()) {
        weight2[f.idx()] = 2 * std::pow(gPrime(areas[f]), 2) - gPrime2(areas[f]) * gPrime(areas[f])
                                                               / std::pow(g(areas[f]), 3);
    }
}

LIM::Scalar LIM::barries() {
    Scalar res = 0;
    for (size_t i = 0; i < F_.cols(); i++) {
        std::vector<Eigen::Matrix<Scalar, 3, 1>> points;
        for (size_t j = 0; j < 3; j++) {
            points.emplace_back(V_.col(F_(j, i)));
        }
        res += phi(c(points[0], points[1], points[2]));
    }
    return res;
}

LIM::Scalar LIM::phi(LIM::Scalar x) {
    if (x <= 0) return INT_MAX;
    else if (x > 0 && x < s_j_) return 1.0 / g(x) - 1;
    else return 0;
}

LIM::Scalar LIM::g(LIM::Scalar x) {
    return 1.0 / std::pow(s_j_, 3) * std::pow(x, 3)
           + 3.0 / std::pow(s_j_, 2) * std::pow(x, 2)
           + 3.0 / s_j_ * x;
}

LIM::Scalar LIM::gPrime(LIM::Scalar x) {
    return 3.0 / std::pow(s_j_, 3) * std::pow(x, 2)
           + 6.0 / std::pow(s_j_, 2) * x
           + 3.0 / s_j_;
}

LIM::Scalar LIM::gPrime2(LIM::Scalar x) {
    return 6.0 / std::pow(s_j_, 3) * x
           + 6.0 / std::pow(s_j_, 2);
}

Eigen::Matrix<LIM::Scalar, Eigen::Dynamic, Eigen::Dynamic> LIM::lambdaPrime() const {
    Eigen::SparseMatrix<Scalar> M_t = M_.transpose();
    return 0.5 * (M_ + M_t) * Vt_calc_;
}

Eigen::Matrix<LIM::Scalar, Eigen::Dynamic, Eigen::Dynamic> LIM::lambdaPrime2() const {
    Eigen::SparseMatrix<Scalar> M_t = M_.transpose();
    return 0.5 * (M_ + M_t);
}

LIM::Scalar LIM::c(const Eigen::Matrix<Scalar, 3, 1> &p1, const Eigen::Matrix<Scalar, 3, 1> &p2,
                   const Eigen::Matrix<Scalar, 3, 1> &p3) {
    return xry_mesh::computeArea<Scalar>(p1, p2, p3);
}

//-------------------------------------------getter and setter---------------------------------------//

const Surface_Mesh::SurfaceMesh &LIM::getMesh() const {
    return mesh_;
}

void LIM::setMesh(const Surface_Mesh::SurfaceMesh &mesh) {
    mesh_ = mesh;
}

const std::vector<std::pair<int, Eigen::Matrix<LIM::Scalar, 3, 1>>> &LIM::getPosConstrains() const {
    return pos_constrains_;
}

void LIM::setPosConstrains(const std::vector<std::pair<int, Eigen::Matrix<Scalar, 3, 1>>> &posConstrains) {
    pos_constrains_ = posConstrains;
}

const ARAPEnergy<Surface_Mesh::Scalar> &LIM::getArapEnergy() const {
    return arap_energy_;
}

void LIM::setArapEnergy(const ARAPEnergy<Surface_Mesh::Scalar> &arapEnergy) {
    arap_energy_ = arapEnergy;

}

void LIM::reset() {

}

