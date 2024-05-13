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
#include "OBB.hpp"
#include "glTools.hpp"
#include "ColorTable.hpp"

float vertices[] = {
    -0.49f, -0.49f, -0.49f,
     0.49f, -0.49f, -0.49f,
     0.49f,  0.49f, -0.49f,
     0.49f,  0.49f, -0.49f,
    -0.49f,  0.49f, -0.49f,
    -0.49f, -0.49f, -0.49f,

    -0.49f, -0.49f,  0.49f,
     0.49f, -0.49f,  0.49f,
     0.49f,  0.49f,  0.49f,
     0.49f,  0.49f,  0.49f,
    -0.49f,  0.49f,  0.49f,
    -0.49f, -0.49f,  0.49f,

    -0.49f,  0.49f,  0.49f,
    -0.49f,  0.49f, -0.49f,
    -0.49f, -0.49f, -0.49f,
    -0.49f, -0.49f, -0.49f,
    -0.49f, -0.49f,  0.49f,
    -0.49f,  0.49f,  0.49f,

     0.49f,  0.49f,  0.49f,
     0.49f,  0.49f, -0.49f,
     0.49f, -0.49f, -0.49f,
     0.49f, -0.49f, -0.49f,
     0.49f, -0.49f,  0.49f,
     0.49f,  0.49f,  0.49f,

    -0.49f, -0.49f, -0.49f,
     0.49f, -0.49f, -0.49f,
     0.49f, -0.49f,  0.49f,
     0.49f, -0.49f,  0.49f,
    -0.49f, -0.49f,  0.49f,
    -0.49f, -0.49f, -0.49f,

    -0.49f,  0.49f, -0.49f,
     0.49f,  0.49f, -0.49f,
     0.49f,  0.49f,  0.49f,
     0.49f,  0.49f,  0.49f,
    -0.49f,  0.49f,  0.49f,
    -0.49f,  0.49f, -0.49f
};

float normals[] = {
    0.0f, 0.0f, -1.0f,  // front face
    0.0f, 0.0f, 1.0f,   // back face
    -1.0f, 0.0f, 0.0f,  // left face
    1.0f, 0.0f, 0.0f,   // right face
    0.0f, -1.0f, 0.0f,  // bottom face
    0.0f, 1.0f, 0.0f    // top face
};


#include "GrainVoxSurf.hpp"

std::vector<GrainVoxSurf> grainSkins;
AABB aabb;
OBB obb;

ColorTable CT;

int main_window;

// flags
int show_background = 0;
int show_grainSkins = 1;
int show_grainSkins_subSkin = 1;
int show_grainSkins_subContact = 1;
int show_grainSkins_subContactWall = 1;
int show_SurfaceNetwork = 1;

int width = 800;
int height = 800;
float wh_ratio = (float)width / (float)height;

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int mouse_start[2];
float view_angle;
float znear;
float zfar;

vec3r eye;
vec3r center;
vec3r up;

// Drawing functions
void drawVoxelCube(float x, float y, float z);
void drawVoxels();
void drawLimits();
void drawSurfaceNetwork();

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

#endif /* end of include guard: SEE_HPP */
