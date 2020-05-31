//
// Created by 徐溶延 on 2020/5/30.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_SURFACEMESHUTILS_H
#define LOCALLYINJECTIVEMAPPINGS_SURFACEMESHUTILS_H
#include <vector>

#include <Eigen/Dense>
#include <SurfaceMesh/SurfaceMesh.h>

namespace xry_mesh {
    /**
     * compute a triangle's area
     * @tparam Scalar double or float
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    template<class Scalar>
    Scalar computeArea(const Eigen::Matrix<Scalar, 3, 1> &p1,
                      const Eigen::Matrix<Scalar, 3, 1> &p2,
                      const Eigen::Matrix<Scalar, 3, 1> &p3);

    /**
     * compute a triangle's area
     * @tparam Scalar double or float
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    template<class Scalar>
    Scalar computeArea(const Surface_Mesh::SurfaceMesh &mesh,
                       const Surface_Mesh::SurfaceMesh::Face &f);

    /**
     * compute areas for each triangle of mesh
     * @tparam Scalar double or float
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    template<class Scalar>
    void computeAreas(const Surface_Mesh::SurfaceMesh &mesh);

    /**
     * get vertices geometry position matrix
     * @tparam Scalar
     * @param mesh
     * @return
     */
    template<class Scalar>
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> getGeometryMatrix(const Surface_Mesh::SurfaceMesh &mesh);

    /**
     * get face id matrix
     * @param mesh
     * @return
     */
    Eigen::Matrix3Xi getTopoMatrix(const Surface_Mesh::SurfaceMesh &mesh);

} // namespace xry_mesh


#endif //LOCALLYINJECTIVEMAPPINGS_SURFACEMESHUTILS_H
