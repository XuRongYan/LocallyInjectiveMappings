//
// Created by 徐溶延 on 2020/5/30.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_LIMOPTIMIZER_H
#define LOCALLYINJECTIVEMAPPINGS_LIMOPTIMIZER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "dbg.h"
#include "EigenUtils.h"

template<typename Scalar>
class LIMOptimizer {
    typedef void (*funcJH)(Eigen::SparseMatrix<Scalar> &J, Eigen::SparseMatrix<Scalar> &H);

    typedef Scalar (*funcEnergy)(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>);

    typedef Scalar (*funcError)(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Vt_i_);
public:
    const Scalar mu_max_, sigma_max_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vt_i_;
    Scalar mu_, sigma_;
    int max_iter;
    funcEnergy funcE_;
    funcError funcError_;
    Eigen::SparseMatrix<Scalar> J_, H_;

    LIMOptimizer() {}

    LIMOptimizer(const Scalar muMax, const Scalar sigmaMax, const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &vtI,
                 funcEnergy funcE) : mu_max_(muMax),
                                                                                    sigma_max_(sigmaMax), Vt_i_(vtI),
                                                                                    funcE_(funcE) {
        sigma_ = sigma_max_;
        mu_ = mu_max_;
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> solve(funcJH getJandH) {
        Scalar err0 = -1, err1 = 0;
        size_t iter = 0;
        while (fabs(err1 - err0) > 1e-6 && iter < max_iter) {
            err0 = err1;
            getJandH(J_, H_);
            lineSearchMu(H_, mu_);
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1> p_i;
            solvePi(p_i);
            lineSearchSigma(p_i, sigma_);
            Vt_i_ = Vt_i_ - sigma_ * p_i;
            err1 = funcError_(Vt_i_);
            dbg(iter, err1);
            iter++;
        }
        return Vt_i_;
    }


private:
    Scalar lineSearchSigma(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &p_i, Scalar sigma_i) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> V_i = Vt_i_ - sigma_i * p_i;
        if (funcE_(V_i) >= funcE_(Vt_i_)) {
            while (funcE_(V_i) >= funcE_(Vt_i_)) {
                sigma_i /= 2;
            }
        } else {
            while (funcE_(V_i) < funcE_(Vt_i_) && sigma_i < sigma_max_) {
                sigma_i *= 2;
            }
            if (sigma_i > sigma_max_) sigma_i = sigma_max_;
        }
        return sigma_i;
    }

    Scalar lineSearchMu(const Eigen::SparseMatrix<Scalar> &H, Scalar mu_i) {
        Eigen::SparseMatrix<Scalar> I(H.size());
        I.setIdentity();
        Eigen::SparseMatrix<Scalar> matrix = H + mu_i * I;
        if (xry_mesh::isSparseMatrixInvertible(matrix)) {
            while (xry_mesh::isSparseMatrixInvertible(matrix)) {
                mu_i /= 2;
            }
            mu_i *= 2;
        } else {
            while (!xry_mesh::isSparseMatrixInvertible(matrix)) {
                mu_i *= 2;

            }
            mu_i /= 2;
        }
        return mu_i;
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> solvePi(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &p_i) {
        Eigen::SparseMatrix<Scalar> I(H_.size());
        I.setIdentity();
        Eigen::SparseMatrix<Scalar> M = H_ + mu_ * I;
        Eigen::PartialPivLU<Eigen::SparseMatrix<Scalar>> partialPivLu(M);
        p_i = partialPivLu.solve(J_);
    }
};


#endif //LOCALLYINJECTIVEMAPPINGS_LIMOPTIMIZER_H
