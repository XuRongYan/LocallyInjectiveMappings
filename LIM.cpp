//
// Created by 徐溶延 on 2020/5/30.
//

#include "LIM.h"

LIM::LIM() {}

LIM::LIM(const Surface_Mesh::SurfaceMesh &mesh) : mesh_(mesh) {}

void LIM::solve() {

}

const Surface_Mesh::SurfaceMesh &LIM::getMesh() const {
    return mesh_;
}

void LIM::setMesh(const Surface_Mesh::SurfaceMesh &mesh) {
    mesh_ = mesh;
}
