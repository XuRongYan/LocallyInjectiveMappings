//
// Created by 徐溶延 on 2020/5/30.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_LIM_H
#define LOCALLYINJECTIVEMAPPINGS_LIM_H
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SurfaceMesh/SurfaceMesh.h>

/**
 * This is the class of Locally Injective Mapping(LIM)
 */
class LIM {
public:
    LIM();

    LIM(const Surface_Mesh::SurfaceMesh &mesh);

    void solve();

    const Surface_Mesh::SurfaceMesh &getMesh() const;

    void setMesh(const Surface_Mesh::SurfaceMesh &mesh);

private:
    Surface_Mesh::SurfaceMesh mesh_;
    std::vector<std::pair<int, Eigen::Matrix3f>> pos_constrains_;
};


#endif //LOCALLYINJECTIVEMAPPINGS_LIM_H
