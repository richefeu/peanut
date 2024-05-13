//   Jean-Yves Delenne & Vicent Richefeu
//   ...Peanut Code -> segment spheres in contact
//    ,+.
//   ((|))
//    )|(
//   ((|))
//    `-'
//
// Compilation :
// g++-13 -o detectContacts detectContacts.cpp -I ~/toofus `pkg-config --cflags --libs libtiff-4` -I/usr/X11R6/include
// -L/usr/X11R6/lib -lX11 -lpthread

#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

#define cimg_use_tiff
#include "CImg/CImg.h"
using namespace cimg_library;

struct ProfilePoint {
  float offset;
  float greyLevel;
};

struct GrainSeed {
  int x;
  int y;
  int z;
  float dst{0.0};
  GrainSeed(int t_x, int t_y, int t_z, int t_dst) : x(t_x), y(t_y), z(t_z), dst(t_dst) {}
};

struct Neighbor {
  int i;
  int j;

  float nx;
  float ny;
  float nz;

  int state{-1};  // -1 = unknown, 0 = not touching, 1 = contact
  float surfaceEstimation{0.0};
};

// Coordinates in voxels
struct Coord {
  int x, y, z;
  Coord(int t_x, int t_y, int t_z) : x(t_x), y(t_y), z(t_z) {}
};



struct GrainVoxSurf {
  std::vector<Coord> voxPos;
  std::vector<int> color; // TODO std::vector<std::vector<int>> color;
};

//
// Algorithme de Bresenham 3D pour déterminer les valeurs sur une droite
//
std::vector<ProfilePoint> getVoxelValuesAlongLine(const CImg<float>& image, int x1, int y1, int z1, int x2, int y2,
                                                  int z2) {

  std::vector<ProfilePoint> result;
  int x0 = x1;
  int y0 = y1;
  int z0 = z1;

  // Calcule le vecteur branche normé
  double ux = x2 - x0;
  double uy = y2 - y0;
  double uz = z2 - z0;
  double norm = sqrt(ux * ux + uy * uy + uz * uz);
  ux /= norm;
  uy /= norm;
  uz /= norm;

  // position sur la droite (projection sur le vecteur u)
  float s = (x1 - x0) * ux + (y1 - y0) * uy + (z1 - z0) * uz;
  result.push_back({s, image(x1, y1, z1, 0)});  // position 0 au début
  int dx = abs(x2 - x1);
  int dy = abs(y2 - y1);
  int dz = abs(z2 - z1);
  int xs;
  int ys;
  int zs;
  if (x2 > x1) {
    xs = 1;
  } else {
    xs = -1;
  }
  if (y2 > y1) {
    ys = 1;
  } else {
    ys = -1;
  }
  if (z2 > z1) {
    zs = 1;
  } else {
    zs = -1;
  }

  if (dx >= dy && dx >= dz) {  // Axe x
    int p1 = 2 * dy - dx;
    int p2 = 2 * dz - dx;
    while (x1 != x2) {
      x1 += xs;
      if (p1 >= 0) {
        y1 += ys;
        p1 -= 2 * dx;
      }
      if (p2 >= 0) {
        z1 += zs;
        p2 -= 2 * dx;
      }
      p1 += 2 * dy;
      p2 += 2 * dz;

      s = (x1 - x0) * ux + (y1 - y0) * uy + (z1 - z0) * uz;
      result.push_back({s, image(x1, y1, z1, 0)});
    }

  } else if (dy >= dx && dy >= dz) {  // Axe y
    int p1 = 2 * dx - dy;
    int p2 = 2 * dz - dy;
    while (y1 != y2) {
      y1 += ys;
      if (p1 >= 0) {
        x1 += xs;
        p1 -= 2 * dy;
      }
      if (p2 >= 0) {
        z1 += zs;
        p2 -= 2 * dy;
      }
      p1 += 2 * dx;
      p2 += 2 * dz;
      s = (x1 - x0) * ux + (y1 - y0) * uy + (z1 - z0) * uz;
      result.push_back({s, image(x1, y1, z1, 0)});
    }

  } else {  // Axe z
    int p1 = 2 * dy - dz;
    int p2 = 2 * dx - dz;
    while (z1 != z2) {
      z1 += zs;
      if (p1 >= 0) {
        y1 += ys;
        p1 -= 2 * dz;
      }
      if (p2 >= 0) {
        x1 += xs;
        p2 -= 2 * dz;
      }
      p1 += 2 * dy;
      p2 += 2 * dx;
      s = (x1 - x0) * ux + (y1 - y0) * uy + (z1 - z0) * uz;
      result.push_back({s, image(x1, y1, z1, 0)});
    }
  }

  return result;
}

//
// Efface la partie externe du cylindre contenant les grains
//
template <typename T>
void eraseCylOutside(CImg<T>& image, int x0, int y0, int z0, int x1, int y1, int z1, int rad, T eraseColor = 0) {

  double vx = x1 - x0;
  double vy = y1 - y0;
  double vz = z1 - z0;
  double norm = sqrt(vx * vx + vy * vy + vz * vz);
  vx /= norm;
  vy /= norm;
  vz /= norm;
  cimg_forXYZ(image, x, y, z) {
    double dx = x - x0;
    double dy = y - y0;
    double dz = z - z0;

    double projection = dx * vx + dy * vy + dz * vz;
    double proj_x = x0 + projection * vx;
    double proj_y = y0 + projection * vy;
    double proj_z = z0 + projection * vz;

    // Calculate the distance to the axis
    double distance_squared = (x - proj_x) * (x - proj_x) + (y - proj_y) * (y - proj_y) + (z - proj_z) * (z - proj_z);

    if (distance_squared > rad * rad || projection < 0 || projection > norm) {
      image(x, y, z, 0) = eraseColor;
    }
  }
}

//
// Fonction pour dessiner une sphère
// utile pour voir si l'on a repéré les centres des grains
//
template <typename T>
void draw_sphere(CImg<T>& image, int center_x, int center_y, int center_z, int radius, T color) {

  // Définir les bornes de la sphère à dessiner
  const int min_x = std::max(0, center_x - radius);
  const int max_x = std::min(image.width() - 1, center_x + radius);
  const int min_y = std::max(0, center_y - radius);
  const int max_y = std::min(image.height() - 1, center_y + radius);
  const int min_z = std::max(0, center_z - radius);
  const int max_z = std::min(image.depth() - 1, center_z + radius);

  // Dessiner les points de la sphère
  for (int x = min_x; x <= max_x; ++x) {
    for (int y = min_y; y <= max_y; ++y) {
      for (int z = min_z; z <= max_z; ++z) {
        int distance_squared =
            (x - center_x) * (x - center_x) + (y - center_y) * (y - center_y) + (z - center_z) * (z - center_z);
        if (distance_squared <= radius * radius) {
          image(x, y, z, 0) = color;
        }
      }
    }
  }
}

//
// Fonction pour ajuster une parabole aux données
// et déterminer la taille de la zone de contact
//
void polyregSecondOrder(std::vector<ProfilePoint>& prof, float& a, float& b, float& c) {

  int N = (int)(prof.size());
  std::vector<double> R = {0.0, 0.0, 0.0};

  double X[5];  // Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3),sigma(xi^4)
  for (int i = 0; i < 5; i++) {
    X[i] = 0.0;
    for (int j = 0; j < N; j++) {
      X[i] = X[i] + pow(prof[j].offset, i);
    }
  }

  double B[3][4];  // B is the Normal matrix(augmented) that will store the equations
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      B[i][j] = X[i + j];
    }
  }

  double Y[3];  // Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)
  for (int i = 0; i < 3; i++) {
    Y[i] = 0.0;
    for (int j = 0; j < N; j++) {
      Y[i] = Y[i] + pow(prof[j].offset, i) * prof[j].greyLevel;
    }
  }

  // Load the values of Y as the last column of B(Normal Matrix but augmented)
  for (int i = 0; i < 3; i++) B[i][3] = Y[i];

  // From now Gaussian Elimination starts
  for (int i = 0; i < 3; i++) {
    for (int k = i + 1; k < 3; k++) {
      if (B[i][i] < B[k][i]) {
        for (int j = 0; j < 4; j++) {
          std::swap(B[i][j], B[k][j]);
        }
      }
    }
  }

  // loop to perform the gauss elimination
  for (int i = 0; i < 2; i++) {
    for (int k = i + 1; k < 3; k++) {
      double t = B[k][i] / B[i][i];
      for (int j = 0; j < 4; j++) B[k][j] = B[k][j] - t * B[i][j];
    }
  }

  // Back-substitution
  for (int i = 2; i >= 0; i--) {
    R[i] = B[i][3];
    for (int j = 0; j < 3; j++)
      if (j != i) R[i] = R[i] - B[i][j] * R[j];
    R[i] = R[i] / B[i][i];
  }

  a = R[2];
  b = R[1];
  c = R[0];
}

//
// Trouve les centres des grains
//
std::vector<GrainSeed> findCenters(CImg<float>& img_dst_map, float rmin, int nmax = 1000) {
  CImg<float> img = img_dst_map;
  std::vector<GrainSeed> seeds;

  for (int n = 0; n < nmax; n++) {

    // TODO: il existe certainement une façon plus rapide de faire ça
    float max_dst = 0.0;
    int xm = 0, ym = 0, zm = 0;
    cimg_forXYZ(img, x, y, z) {
      float dst = img(x, y, z, 0);

      if (max_dst < dst) {
        max_dst = dst;
        xm = x;
        ym = y;
        zm = z;
      }
    }

    if (max_dst <= rmin) {
      break;
    } else {
      std::cout << "Sphere found at position (" << xm << ", " << ym << ", " << zm << "), radius = " << max_dst << "\n";
      draw_sphere<float>(img, xm, ym, zm, (int)floor(max_dst), 0.0);
      seeds.push_back(GrainSeed(xm, ym, zm, max_dst));
    }
  }
  std::cout << seeds.size() << " spheres found\n";

  return seeds;
}

//
// Recherche les voisins (en faisant l'hypothèse que les grains sont sphériques)
//
std::vector<Neighbor> findNeighbors(std::vector<GrainSeed>& seeds, int dmin = 10) {

  std::vector<Neighbor> neighbors;

  for (size_t i = 0; i < seeds.size(); i++) {
    for (size_t j = i + 1; j < seeds.size(); j++) {
      float xi = (float)seeds[i].x;
      float yi = (float)seeds[i].y;
      float zi = (float)seeds[i].z;
      float Ri = seeds[i].dst;
      float xj = (float)seeds[j].x;
      float yj = (float)seeds[j].y;
      float zj = (float)seeds[j].z;
      float Rj = seeds[j].dst;

      float dx = xj - xi;
      float dy = yj - yi;
      float dz = zj - zi;
      float distance = sqrt(dx * dx + dy * dy + dz * dz) - Ri - Rj;
      if (distance < dmin) {
        Neighbor C;
        C.i = i;
        C.j = j;
        neighbors.push_back(C);
      }
    }
  }
  return neighbors;
}

//
// 26-ways floodfill using stack routines
//
// image        Image sur laquelle on rempli
// distMap.      ...
// x,y,z        Position du point de remplissage
// nx, ny, nz   Normal au plan du cylindre épais
// distMaxAxe   Demi-épaisseur du cylindre
// distMaxRad   Rayon maximum du cylindre (on arrête de remplir si la distance est nulle)
// newColor1    Couleur de remplisable du coté du début de n
// newColor2    Couleur de remplisable du coté de la fin de n
template <typename T>
bool floodFill26SPlanDistRestricted(CImg<T>& image, CImg<float>& distMap, int x, int y, int z, float nx, float ny,
                                    float nz, float distMaxAxe, float distMaxRad, T newColor1, T newColor2) {
  static int deltaX[] = {1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, -1, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1};
  static int deltaY[] = {0, 1, 0, -1, 1, 1, 1, 0, 0, 1, 1, -1, -1, 0, -1, 0, 1, -1, -1, -1, 0, 0, -1, -1, 1, 1};
  static int deltaZ[] = {0, 0, 1, 0, 0, -1, 1, -1, 1, -1, 1, 1, -1, 0, 0, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};

  T oldColor = image(x, y, z);

  if (newColor1 == oldColor || newColor2 == oldColor) {
    return false;
  }

  int width = image.width();
  int height = image.height();
  int depth = image.depth();

  // std::deque should faster than std::vector (better cache locality)
  std::deque<Coord> stack;
  std::deque<T> newColorStack;

  stack.push_back(Coord(x, y, z));
  newColorStack.push_back(newColor2);

  while (!stack.empty()) {
    Coord pos = stack.back();
    stack.pop_back();
    T color = newColorStack.back();
    newColorStack.pop_back();

    image(pos.x, pos.y, pos.z) = color;

    for (int i = 0; i < 26; i++) {
      int nextX = pos.x + deltaX[i];
      int nextY = pos.y + deltaY[i];
      int nextZ = pos.z + deltaZ[i];

      float px = (nextX - x);
      float py = (nextY - y);
      float pz = (nextZ - z);
      float distAxe = px * nx + py * ny + pz * nz;
      float rx = px - distAxe * nx;
      float ry = py - distAxe * ny;
      float rz = pz - distAxe * nz;
      float distRad2 = rx * rx + ry * ry + rz * rz;

      if (nextX < 0 || nextX >= width || nextY < 0 || nextY >= height || nextZ < 0 || nextZ >= depth ||
          fabs(distAxe) > distMaxAxe || distRad2 > distMaxRad * distMaxRad) {
        continue;
      }

      if (image(nextX, nextY, nextZ) == oldColor && distMap(nextX, nextY, nextZ) > 0.1) {
        if (distAxe >= 0) {
          newColorStack.push_back(newColor2);
        } else {
          newColorStack.push_back(newColor1);
        }

        stack.push_back(Coord(nextX, nextY, nextZ));
      }
    }
  }

  return true;
}

//
// Lit les valeurs des distances sur une ligne connectant les centres
// détermine la postion du contact et ajoute un disque au niveau du contact
//
template <typename T>
void addContactCap(CImg<T>& image, CImg<float>& distance_map, std::vector<GrainSeed>& seeds, Neighbor& N,
                   int contactCapThickness, int contactCapRadiusInc) {
  int i = N.i;
  int j = N.j;
  int x1 = seeds[i].x;
  int y1 = seeds[i].y;
  int z1 = seeds[i].z;
  int x2 = seeds[j].x;
  int y2 = seeds[j].y;
  int z2 = seeds[j].z;

  N.nx = x2 - x1;
  N.ny = y2 - y1;
  N.nz = z2 - z1;
  float norm_vec_n = sqrt(N.nx * N.nx + N.ny * N.ny + N.nz * N.nz);
  N.nx /= norm_vec_n;
  N.ny /= norm_vec_n;
  N.nz /= norm_vec_n;

  // Profil de valeurs des voxels le long de la ligne droite reliant les centres
  std::vector<ProfilePoint> voxelValues = getVoxelValuesAlongLine(distance_map, x1, y1, z1, x2, y2, z2);

  // Si on trouve du vide entre les grains, c'est qu'il n'y a pas de contact
  for (size_t i = 0; i < voxelValues.size(); i++) {
    if (voxelValues[i].greyLevel < 0.1) {
      std::cout << "No contact between " << i + 1 << " and " << j + 1 << '\n';
      N.state = 0;
      return;
    }
  }

  float a, b, c;
  polyregSecondOrder(voxelValues, a, b, c);
  float min_parabola_s = -b / (2.0 * a);
  float min_parabola_dst = -b * b / (4.0 * a) + c;
  if (min_parabola_dst < 0.1) {
    std::cout << "No contact between " << i + 1 << " and " << j + 1 << '\n';
    N.state = 0;
    return;
  }

  N.surfaceEstimation = M_PI * min_parabola_dst * min_parabola_dst;
  N.state = 1;

  float Rmax = min_parabola_dst + contactCapRadiusInc;

  float ux = x2 - x1;
  float uy = y2 - y1;
  float uz = z2 - z1;
  float norm = sqrt(ux * ux + uy * uy + uz * uz);
  ux /= norm;
  uy /= norm;
  uz /= norm;

  float xI = min_parabola_s * ux + x1;
  float yI = min_parabola_s * uy + y1;
  float zI = min_parabola_s * uz + z1;
  std::cout << "add Cap to contact between " << i + 1 << " and " << j + 1 << '\n';
  floodFill26SPlanDistRestricted<int>(image, distance_map, xI, yI, zI, ux, uy, uz, contactCapThickness * 0.5, Rmax,
                                      i + 1, j + 1);

  // std::cout << i << " " << j << "\n";
  // std::cout << "Surface contact : 2 PI R = " << M_PI * min_parabola_y * min_parabola_y
  //           << ", surface réelle = " << N.surfaceEstimation << "\n";
}

//
// 26-ways floodfill using stack routines. This version does nothing else than floodfilling
//
template <typename T>
int floodFill26S(CImg<T>& image, int x, int y, int z, float distMax, T newColor) {
  static int deltaX[] = {1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, -1, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1};
  static int deltaY[] = {0, 1, 0, -1, 1, 1, 1, 0, 0, 1, 1, -1, -1, 0, -1, 0, 1, -1, -1, -1, 0, 0, -1, -1, 1, 1};
  static int deltaZ[] = {0, 0, 1, 0, 0, -1, 1, -1, 1, -1, 1, 1, -1, 0, 0, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};

  T oldColor = image(x, y, z);

  if (newColor == oldColor) {
    return 0;
  }

  int width = image.width();
  int height = image.height();
  int depth = image.depth();
  float distMax2 = distMax * distMax;

  // std::deque should faster than std::vector (better cache locality)
  std::deque<Coord> stack;

  stack.push_back(Coord(x, y, z));
  int nbPix = 0;
  while (!stack.empty()) {
    Coord pos = stack.back();
    stack.pop_back();
    image(pos.x, pos.y, pos.z) = newColor;
    nbPix++;

    for (int i = 0; i < 26; i++) {
      int nextX = pos.x + deltaX[i];
      int nextY = pos.y + deltaY[i];
      int nextZ = pos.z + deltaZ[i];

      float dx = nextX - x;
      float dy = nextY - y;
      float dz = nextZ - z;
      float dist2 = dx * dx + dy * dy + dz * dz;

      if (nextX < 0 || nextX >= width || nextY < 0 || nextY >= height || nextZ < 0 || nextZ >= depth ||
          dist2 > distMax2) {
        continue;
      }

      if (image(nextX, nextY, nextZ) == oldColor) {
        stack.push_back(Coord(nextX, nextY, nextZ));
      }
    }
  }

  return nbPix;
}

//
//
//
template <typename T>
GrainVoxSurf getGrainVoxSurf(CImg<T>& image, int center_x, int center_y, int center_z, int radius) {
  static int deltaX[] = {1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, -1, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1};
  static int deltaY[] = {0, 1, 0, -1, 1, 1, 1, 0, 0, 1, 1, -1, -1, 0, -1, 0, 1, -1, -1, -1, 0, 0, -1, -1, 1, 1};
  static int deltaZ[] = {0, 0, 1, 0, 0, -1, 1, -1, 1, -1, 1, 1, -1, 0, 0, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};

  GrainVoxSurf GRAIN;
  T colorGrain = image(center_x, center_y, center_z);
  std::cout << "colorGrain = " << colorGrain<< '\n';

  // Définir les bornes de la sphère à dessiner
  const int min_x = std::max(1, center_x - radius + 1);
  const int max_x = std::min(image.width() - 2, center_x + radius - 1);
  const int min_y = std::max(1, center_y - radius + 1);
  const int max_y = std::min(image.height() - 2, center_y + radius - 1);
  const int min_z = std::max(1, center_z - radius + 1);
  const int max_z = std::min(image.depth() - 2, center_z + radius - 1);

  for (int x = min_x; x <= max_x; ++x) {
    for (int y = min_y; y <= max_y; ++y) {
      for (int z = min_z; z <= max_z; ++z) {

        for (int i = 0; i < 26; i++) {
          int nextX = x + deltaX[i];
          int nextY = y + deltaY[i];
          int nextZ = z + deltaZ[i];

          if (image(x, y, z) == colorGrain && image(nextX, nextY, nextZ) != colorGrain) {
            GRAIN.voxPos.push_back(Coord(x, y, z));
            GRAIN.color.push_back(image(nextX, nextY, nextZ));
            break;
          }
        }
        
      } // z
    } // y
  } // x

  return GRAIN;
}

//
// identification des grains en utilisant le floodfill
//
template <typename T>
const CImg<int> addLabels(const CImg<T>& image, std::vector<GrainSeed>& seeds, float radiusInc) {

  // Crée une image int pour avoir des labels qui vont au delà de 255
  CImg<int> labels(image.width(), image.height(), image.depth(), 1);
  cimg_forXYZ(image, x, y, z) { labels(x, y, z) = static_cast<int>(image(x, y, z)); }

  std::cout << "Add labels for " << seeds.size() << " particles\n";

  for (size_t s = 0; s < seeds.size(); s++) {
    std::cout << "Label particle " << s << " / " << seeds.size() << " ... ";
    int nbPix = floodFill26S<int>(labels, seeds[s].x, seeds[s].y, seeds[s].z, seeds[s].dst + radiusInc, s + 1);
    std::cout << nbPix << " voxels\n";
  }

  return labels;
}

template <typename T>
void clean(CImg<T>& image, T from, T to) {
  cimg_forXYZ(image, x, y, z) {
    if (image(x, y, z) == from) {
      image(x, y, z) = to;
    }
  }
}

//
// Export pour Paraview (.vti)
//
void exportToVTI(const CImg<unsigned char>& image, const char* filename) {
  std::cout << "Save " << filename << std::endl;
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
    return;
  }
  int W = image.width();
  int H = image.height();
  int D = image.depth();
  std::cout << W << " " << H << " " << D << "\n";

  // Écriture de l'en-tête du fichier VTI
  outfile << "<?xml version=\"1.0\"?>" << std::endl;
  outfile << "<VTKFile type=\"ImageData\" version=\"1.0\">" << std::endl;
  outfile << "<ImageData WholeExtent=\"0 " << W - 1 << " 0 " << H - 1 << " 0 " << D - 1 << "\" ";
  outfile << "Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
  outfile << "<Piece Extent=\"0 " << W - 1 << " 0 " << H - 1 << " 0 " << D - 1 << "\">" << std::endl;
  outfile << "<PointData Scalars=\"ImageScalars\">" << std::endl;
  // outfile << "<CellData Scalars=\"ImageScalars\">" << std::endl;
  outfile << "<DataArray type=\"UInt8\" Name=\"ImageScalars\" format=\"ascii\">" << std::endl;

  // Écriture des valeurs de pixel
  for (int z = 0; z < D; ++z) {
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        outfile << (int)image(x, y, z) << " ";
      }
      outfile << std::endl;
    }
  }

  // Fermeture des balises et du fichier
  outfile << "</DataArray>" << std::endl;
  outfile << "</PointData>" << std::endl;
  outfile << "</Piece>" << std::endl;
  outfile << "</ImageData>" << std::endl;
  outfile << "</VTKFile>" << std::endl;

  std::cout << "Done\n";
}

class Peanut {
 public:
  int eraseCylinder{0};

  // bords du cylindre interne
  float xminCyl{156};
  float xmaxCyl{847};
  float yminCyl{165};
  float ymaxCyl{856};
  float zminCyl{31};
  float zmaxCyl{364};

  float rmin{40};  // rayon minimum des grains
  int nmax{1000};  // nombre maximum de grains à chercher

  // distance de recherche des voisins
  int distance_max{10};

  int contactCapThickness{5};
  int contactCapRadiusInc{2};

  int radiusInc{5};

  std::string filenameTIF{"troisSpheres.tif"};

  int blur{0};
  float blurSigma{1.5};

  int threshold{0};
  float thresholdValue{50};

  int exportVTI{0};
  int exportSeed{0};
  int exportContacts{0};

  int show_after_normalisation{1};
  int show_distance_map{1};
  int show_marked_contacts{1};
  int show_labels{1};

  void printHeader() {

    std::cout << "          ,+.  \n";
    std::cout << "         ((|)) \n";
    std::cout << "          )|(. \n";
    std::cout << "         ((|)) \n";
    std::cout << "          `-'  \n\n";
    std::cout << "Peanut Voxel Analyzer\n\n";
  }

  void readCommands(const char* filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "@readCommands, Cannot read " << filename << '\n';
    }

    std::string token;
    file >> token;

    while (file.good()) {
      if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
        getline(file, token);  // ignore the rest of the current line
        file >> token;         // next token
        continue;
      } else if (token == "eraseCylinder") {
        file >> eraseCylinder;
      } else if (token == "xminCyl") {
        file >> xminCyl;
      } else if (token == "xmaxCyl") {
        file >> xmaxCyl;
      } else if (token == "yminCyl") {
        file >> yminCyl;
      } else if (token == "ymaxCyl") {
        file >> ymaxCyl;
      } else if (token == "zminCyl") {
        file >> zminCyl;
      } else if (token == "zmaxCyl") {
        file >> zmaxCyl;
      } else if (token == "rmin") {
        file >> rmin;
      } else if (token == "nmax") {
        file >> nmax;
      } else if (token == "distance_max") {
        file >> distance_max;
      } else if (token == "contactCapThickness") {
        file >> contactCapThickness;
      } else if (token == "contactCapRadiusInc") {
        file >> contactCapRadiusInc;
      } else if (token == "radiusInc") {
        file >> radiusInc;
      } else if (token == "filenameTIF") {
        file >> filenameTIF;
      } else if (token == "exportVTI") {
        file >> exportVTI;
      } else if (token == "exportSeed") {
        file >> exportSeed;
      } else if (token == "exportContacts") {
        file >> exportContacts;
      } else if (token == "blur") {
        file >> blur;
      } else if (token == "blurSigma") {
        file >> blurSigma;
      } else if (token == "threshold") {
        file >> threshold;
      } else if (token == "thresholdValue") {
        file >> thresholdValue;
      } else if (token == "show_after_normalisation") {
        file >> show_after_normalisation;
      } else if (token == "show_distance_map") {
        file >> show_distance_map;
      } else if (token == "show_marked_contacts") {
        file >> show_marked_contacts;
      } else if (token == "show_labels") {
        file >> show_labels;
      } else if (token == "GenerateSphereSample") {

        std::string filename;
        int width, height, depth;
        size_t nb;
        file >> filename >> width >> height >> depth;
        std::cout << filename << " " << width << " " << height << " " << depth << "\n";
        CImg<unsigned char> image(width, height, depth, 1, 0);
        file >> nb;
        for (size_t i = 0; i < nb; i++) {
          int x, y, z, r;
          file >> x >> y >> z >> r;
          std::cout << x << " " << y << " " << z << " " << r << "\n";
          draw_sphere<unsigned char>(image, x, y, z, r, 255);
        }
        std::cout << "save " << filename.c_str() << "\n";
        image.save_tiff(filename.c_str(), 8);

      } else if (token == "exit") {
        exit(0);
      } else {
        std::cerr << "Unknown token: " << token << std::endl;
        exit(0);
      }

      file >> token;
    }
  }

  // This is the procedure
  void run() {

    CImg<int> image(filenameTIF.c_str());
    image.resize(image.width(), image.height(), image.depth(), 1);  // we only need a single channel

    // d'éventuels
    if (blur == 1) {
      std::cout << "Blur\n" << std::flush;
      image.blur(blurSigma);
    }

    if (threshold == 1) {
      std::cout << "Threshold\n" << std::flush;
      image.threshold(thresholdValue);
    }

    std::cout << "Normalize in range 0-255\n" << std::flush;
    image.normalize(0, nmax);
    if (show_after_normalisation == 1) {
      image.display();
    }

    // pour enlever l'exterieur du cylindre
    if (eraseCylinder == 1) {
      int x0 = (xminCyl + xmaxCyl) * 0.5;
      int y0 = (yminCyl + ymaxCyl) * 0.5;
      int z0 = zminCyl;
      int x1 = x0;
      int y1 = y0;
      int z1 = zmaxCyl;
      int rad = ((xmaxCyl - xminCyl) * 0.5 + (ymaxCyl - yminCyl) * 0.5) * 0.5;
      eraseCylOutside<int>(image, x0, y0, z0, x1, y1, z1, rad, nmax - 1);
    }

    // Carte de distance
    std::cout << "Distance map (can be long)\n" << std::flush;
    CImg<float> distance_map = image.get_distance(0, 2);  // 0 = distance to closest 0, 2 = euclidian distance
    if (show_distance_map == 1) {
      distance_map.display();
    }

    std::cout << "Find centers\n" << std::flush;
    std::vector<GrainSeed> seeds = findCenters(distance_map, rmin, nmax);  // rmin = 15.0
    std::cout << "Number found centers = " << seeds.size() << '\n';

    std::cout << "Find neighbors\n" << std::flush;
    std::vector<Neighbor> neighbors = findNeighbors(seeds, distance_max);
    std::cout << "Number found neigbors = " << neighbors.size() << '\n';

    std::cout << "Mark contacts\n" << std::flush;
    for (size_t i = 0; i < neighbors.size(); i++) {
      addContactCap<int>(image, distance_map, seeds, neighbors[i], contactCapThickness, contactCapRadiusInc);
    }
    if (show_marked_contacts == 1) {
      image.display();
    }

    if (exportSeed == 1) {
      std::ofstream partFile("particles.txt");
      partFile << "# x y z r" << '\n';
      for (size_t i = 0; i < seeds.size(); i++) {
        partFile << seeds[i].x << ' ' << seeds[i].y << ' ' << seeds[i].z << ' ' << seeds[i].dst << '\n';
      }
    }

    if (exportContacts == 1) {
      std::ofstream ctcFile("contacts.txt");
      ctcFile << "# id_i id_j S_contact nx ny" << '\n';
      for (size_t i = 0; i < neighbors.size(); i++) {
        ctcFile << neighbors[i].i << ' ' << neighbors[i].j << ' ' << neighbors[i].surfaceEstimation << ' '
                << neighbors[i].nx << ' ' << neighbors[i].ny << ' ' << neighbors[i].nz << '\n';
      }
    }

    std::cout << "Labelize grains\n" << std::flush;
    CImg<int> labels;
    labels = addLabels(image, seeds, (float)contactCapRadiusInc);

    // nettoyage des couleurs nmax -> 0 pour transformer les grains trop petits en vide
    clean(labels, nmax, 0);

    image.assign();  // free memory
    if (show_labels == 1) {
      labels.display();
    }

    // construire les peaux de grain
    std::vector<GrainVoxSurf> grainSkins;
    for (size_t i = 0; i < seeds.size(); i++) {
      grainSkins.push_back(
          getGrainVoxSurf(labels, seeds[i].x, seeds[i].y, seeds[i].z, (int)std::ceil(seeds[i].dst) + 5));
    }

    // provisoire
    std::ofstream bidon("skins.txt");
    // bidon << "TITLE = \"Skins of particles\"\n";
    // bidon << "VARIABLES : \"x\", \"y\", \"z\", \"c\"\n";
    // bidon << "     F:POINT\n";
    bidon << grainSkins.size() << '\n';
    for (size_t g = 0; g < grainSkins.size(); g++) {
      bidon << grainSkins[g].voxPos.size() << '\n';
      for (size_t i = 0; i < grainSkins[g].voxPos.size(); i++) {
        bidon << grainSkins[g].voxPos[i].x << " " << grainSkins[g].voxPos[i].y << " " << grainSkins[g].voxPos[i].z
              << " " << grainSkins[g].color[i] << "\n";
      }
      bidon << '\n';
    }

    if (exportVTI == 1) {
      exportToVTI(labels, "output_image.vti");
    }
  }
};

int main(int argc, char const* argv[]) {

  Peanut peanut;

  if (argc < 2 || argc > 2) {
    std::cout << "usage: " << argv[0] << " command-file\n";
    return 0;
  } else {
    peanut.readCommands(argv[1]);
  }

  peanut.printHeader();
  peanut.run();

  return 0;
}
