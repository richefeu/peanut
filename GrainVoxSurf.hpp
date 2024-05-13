#pragma once

#include <iostream>
#include <set>

#include "mat9.hpp"

struct GrainVoxSurf {
  std::vector<vec3i> voxPos;
  std::vector<int> color;

  vec3r baryCenter;

  std::vector<vec3r> contactPos;
  std::vector<vec3r> contactNormal;
  std::vector<double> contactArea;

  void computeBaryCenter() {
    if (voxPos.empty()) {
      std::cout << "@GrainVoxSurf::computeBaryCenter, EMPTY!\n";
      return;
    }

    baryCenter.reset();
    for (size_t i = 0; i < voxPos.size(); i++) {
      baryCenter.x += (float)voxPos[i].x;
      baryCenter.y += (float)voxPos[i].y;
      baryCenter.z += (float)voxPos[i].z;
    }
    baryCenter /= (float)voxPos.size();
  }

  void computeContacts() {
    std::set<int> cols;
    for (size_t i = 0; i < color.size(); i++) {
      if (color[i] > 0) {
        cols.insert(color[i]);
      }
    }
    //int nbCtc = (int)cols.size();
    // std::cout << "nbCtc = " << nbCtc << '\n';

    for (auto it = cols.begin(); it != cols.end(); ++it) {
      int current_color = *it;
      std::vector<vec3r> pset;
      vec3r c;
      for (size_t i = 0; i < voxPos.size(); i++) {
        if (color[i] == current_color) {
          vec3r vv(voxPos[i].x, voxPos[i].y, voxPos[i].z);
          pset.push_back(vv);
          c += vv;
        }
      }
      c /= (double)pset.size();
      mat9r C = CovarianceMatrix(pset);
      mat9r V;
      vec3r D;
      C.sorted_sym_eigen(V, D);
      vec3r n(V.xx,V.yx,V.zx);
      
      contactPos.push_back(c);
      contactNormal.push_back(n);
      contactArea.push_back((double)pset.size());
    }
  }
};
