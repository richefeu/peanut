#ifndef SEE_HPP
#define SEE_HPP

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#include <GL/freeglut.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>
#include <fstream>
#include <iostream>

#include "AABB.hpp"
#include "glTools.hpp"

struct Coord {
  int x, y, z;
  Coord(int t_x, int t_y, int t_z) : x(t_x), y(t_y), z(t_z) {}
};

struct GrainVoxSurf {
  std::vector<Coord> voxPos;
  std::vector<int> color; // TODO std::vector<std::vector<int>> color;
};

std::vector<GrainVoxSurf> grainSkins;
AABB aabb;


int main_window;

// flags
int show_background = 0;

int width = 800;
int height = 800;
float wh_ratio = (float)width / (float)height;

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int display_mode = 0;  // sample or slice rotation
int mouse_start[2];
float view_angle;
float znear;
float zfar;
GLfloat Rot_Matrix[16];
GLfloat max_length;

vec3r eye;
vec3r center;
vec3r up;

// Drawing functions

void clear_background();
void drawVoxels();


// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void buildMenu();
void printHelp();
vec3r rotatePoint(vec3r const& p, vec3r const& center, vec3r const& axis, double theta);
void adjust_clipping_plans();
void fit_view();
bool fileExists(const char* fileName);
bool try_to_readConf(int num);
int screenshot(const char* filename);



#endif /* end of include guard: SEE_HPP */
