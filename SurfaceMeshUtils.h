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
    template<typename Scalar>
    Scalar computeArea(const Eigen::Matrix<Scalar, 3, 1> &p1,
                      const Eigen::Matrix<Scalar, 3, 1> &p2,
                      const Eigen::Matrix<Scalar, 3, 1> &p3) {
        return 0.5 * (p2 - p1).dot(p3 - p1);
    }

    /**
     * compute a triangle's area
     * @tparam Scalar double or float
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    template<typename Scalar>
    Scalar computeArea(const Surface_Mesh::SurfaceMesh &mesh,
                       const Surface_Mesh::SurfaceMesh::Face &f) {
        std::vector<Eigen::Matrix<Scalar, 3, 1>> points;
        for (auto v : mesh.vertices(f)) {
            points.emplace_back(mesh.position(v));
        }
        assert(points.size() == 3);
        return computeArea<Scalar>(points[0], points[1], points[2]);
    }

    /**
     * compute areas for each triangle of mesh
     * @tparam Scalar double or float
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    template<class Scalar>
    void computeAreas(Surface_Mesh::SurfaceMesh &mesh) {
        auto areas = mesh.add_face_property<Scalar>("f:areas");
        for (auto f : mesh.faces()) {
            areas[f] = computeArea<Scalar>(mesh, f);
        }
    }

    /**
     * get vertices geometry position matrix
     * @tparam Scalar
     * @param mesh
     * @return
     */
    template<typename Scalar>
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> getGeometryMatrix(const Surface_Mesh::SurfaceMesh &mesh) {
        std::vector<Scalar> points;
        for (const auto &p : mesh.points()) {
            points.push_back(p.y());
            points.push_back(p.y());
            points.push_back(p.y());
        }
        return Eigen::Map<Eigen::Matrix<Scalar, 3, Eigen::Dynamic>>(points.data(), 3, points.size() / 3);
    }

    /**
     * get face id matrix
     * @param mesh
     * @return
     */
    Eigen::Matrix3Xi getTopoMatrix(const Surface_Mesh::SurfaceMesh &mesh);

} // namespace xry_mesh


#endif //LOCALLYINJECTIVEMAPPINGS_SURFACEMESHUTILS_H
