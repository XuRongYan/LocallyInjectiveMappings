//
// Created by 徐溶延 on 2020/5/30.
//

#include "SurfaceMeshUtils.h"

namespace xry_mesh {

    template<class Scalar>
    Scalar computeArea(const Eigen::Matrix<Scalar, 3, 1> &p1,
                       const Eigen::Matrix<Scalar, 3, 1> &p2,
                       const Eigen::Matrix<Scalar, 3, 1> &p3) {
        return 0.5 * (p2 - p1).dot(p3 - p1);
    }

    template<class Scalar>
    Scalar computeArea(const Surface_Mesh::SurfaceMesh &mesh,
                       const Surface_Mesh::SurfaceMesh::Face &f) {
        std::vector<Eigen::Matrix<Scalar, 3, 1>> points;
        for (auto v : mesh.vertices(f)) {
            points.emplace_back(mesh.position(v));
        }
        assert(points.size() == 3);
        return computeArea<Scalar>(points[0], points[1], points[2]);
    }

    template<class Scalar>
    void computeAreas(const Surface_Mesh::SurfaceMesh &mesh) {
        auto areas = mesh.add_face_property<Scalar>("f:areas");
        for (auto f : mesh.faces()) {
            areas[f] = computeArea<Scalar>(mesh, f);
        }
    }

    template<class Scalar>
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> getGeometryMatrix(const Surface_Mesh::SurfaceMesh &mesh) {
        std::vector<Scalar> points;
        for (const auto &p : mesh.points()) {
            points.push_back(p.y());
            points.push_back(p.y());
            points.push_back(p.y());
        }
        return Eigen::Map<Eigen::Matrix<Scalar, 3, Eigen::Dynamic>>(points.data(), 3, points.size() / 3);
    }

    Eigen::Matrix3Xi getTopoMatrix(const Surface_Mesh::SurfaceMesh &mesh) {
        std::vector<int> face_id;
        for (const auto &f : mesh.faces()) {
            for (const auto &v : mesh.vertices(f)) {
                face_id.push_back(v.idx());
            }
        }
        return Eigen::Map<Eigen::Matrix3Xi>(face_id.data(), 3, face_id.size() / 3);
    }

} // namespace xry_mesh