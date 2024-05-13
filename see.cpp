#include "see.hpp"

void printHelp() {
  using namespace std;
  // cout << "" << endl;
  cout << endl;
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

    case 'b':
      show_background = 1 - show_background;
      break;

    case 'h':
      printHelp();
      break;

    case 'p':
      show_grainSkins = 1 - show_grainSkins;
      break;

    case 'q':
      exit(0);
      break;

    case 's':
      show_SurfaceNetwork = 1 - show_SurfaceNetwork;
      break;

    case '1':
      show_grainSkins_subSkin = 1 - show_grainSkins_subSkin;
      break;
    case '2':
      show_grainSkins_subContact = 1 - show_grainSkins_subContact;
      break;
    case '3':
      show_grainSkins_subContactWall = 1 - show_grainSkins_subContactWall;
      break;

    case '=': {
      fit_view();
      adjust_clipping_plans();
    } break;
  };

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
          mouse_mode = PAN;
        else
          mouse_mode = ROTATION;
        break;
      case GLUT_MIDDLE_BUTTON:
        mouse_mode = ZOOM;
        break;
    }
  }
}

vec3r rotatePoint(vec3r const& p, vec3r const& center_, vec3r const& axis, double theta) {
  double const c = cos(theta), s = sin(theta);
  double const C = 1.0 - c;
  vec3r tmp = p - center_;
  return center + vec3r(tmp[0] * (axis[0] * axis[0] * C + c) + tmp[1] * (axis[0] * axis[1] * C - axis[2] * s) +
                            tmp[2] * (axis[0] * axis[2] * C + axis[1] * s),
                        tmp[0] * (axis[1] * axis[0] * C + axis[2] * s) + tmp[1] * (axis[1] * axis[1] * C + c) +
                            tmp[2] * (axis[1] * axis[2] * C - axis[0] * s),
                        tmp[0] * (axis[2] * axis[0] * C - axis[1] * s) +
                            tmp[1] * (axis[2] * axis[1] * C + axis[0] * s) + tmp[2] * (axis[2] * axis[2] * C + c));
}

void motion(int x, int y) {
  if (mouse_mode == NOTHING) return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;
  double length;
  vec3r axis;

  switch (mouse_mode) {

    case ROTATION:
      axis = (cross(up, center - eye));
      axis.normalize();
      eye = rotatePoint(eye, center, up, -dx * M_PI);
      eye = rotatePoint(eye, center, axis, dy * M_PI);
      up = (rotatePoint((center + up), center, axis, dy * M_PI) - center);
      up.normalize();
      break;

    case ZOOM:
      eye = center + (eye - center) * (dy + 1.0);
      break;

    case PAN:
      length = (eye - center).length() * tan(view_angle * M_PI / 360.0) * 2.0;
      axis = cross(up, center - eye);
      axis.normalize();
      center = center + axis * dx * length * 0.8;
      center = center + up * dy * length;
      break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  display();
}

void display() {
  // clear_background();
  glTools::clearBackground((bool)show_background);
  adjust_clipping_plans();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  drawLimits();
  if (show_grainSkins == 1) {
    drawVoxels();
  }
  if (show_SurfaceNetwork == 1) {
    drawSurfaceNetwork();
  }

  glFlush();
  glutSwapBuffers();
}

void adjust_clipping_plans() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  wh_ratio = (float)width / (float)height;
  float zf = (float)((eye - center).normalize());
  vec3r mx = aabb.max - aabb.min;
  GLfloat max_length = (GLfloat)norm(mx);
  znear = zf - 0.5f * max_length;
  float close_dst = 0.1f * zf;
  if (znear < close_dst) znear = close_dst;
  zfar = zf + 0.5f * max_length;
  gluPerspective(view_angle, wh_ratio, znear, zfar);
  glMatrixMode(GL_MODELVIEW);
}

void fit_view() {
  vec3r dir = (eye - center);
  vec3r diag = aabb.max - aabb.min;
  dir.normalize();
  center = 0.5 * (aabb.max + aabb.min);
  GLfloat d = 0.5f * (GLfloat)diag.length() / (GLfloat)(atan(view_angle * M_PI / 360.0));
  eye = center + d * dir;
}

void reshape(int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, width, height);

  adjust_clipping_plans();
  glutPostRedisplay();
}

void drawVoxelCube(float x, float y, float z) {
  glTranslatef(x, y, z);
  glDrawArrays(GL_QUADS, 0, 24);
  glTranslatef(-x, -y, -z);
}

void drawVoxels() {

  if (mouse_mode != NOTHING) return;

  glEnable(GL_LIGHTING);

  for (size_t g = 0; g < grainSkins.size(); g++) {
    for (size_t v = 0; v < grainSkins[g].voxPos.size(); v++) {
      if (grainSkins[g].color[v] == 0) {
        color4f col;
        CT.getColor4f(g, &col);
        glColor4f(col.r, col.g, col.b, 1.0f);
        if (show_grainSkins_subSkin == 1) {
          drawVoxelCube(grainSkins[g].voxPos[v].x, grainSkins[g].voxPos[v].y, grainSkins[g].voxPos[v].z);
        }
      } else if (grainSkins[g].color[v] < (int)grainSkins.size()) {
        glColor4f(1.0f, 0.0f, 0.0f, 0.5f);
        if (show_grainSkins_subContact == 1) {
          drawVoxelCube(grainSkins[g].voxPos[v].x, grainSkins[g].voxPos[v].y, grainSkins[g].voxPos[v].z);
        }
      } else {
        glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
        if (show_grainSkins_subContactWall == 1) {
          drawVoxelCube(grainSkins[g].voxPos[v].x, grainSkins[g].voxPos[v].y, grainSkins[g].voxPos[v].z);
        }
      }
    }
  }
}

void drawLimits() {
  glColor3f(1.0f, 0.0f, 0.0f);
  glShape::obb(obb);
}

void drawSurfaceNetwork() {

  glEnable(GL_LIGHTING);
  glColor3f(0.5f, 0.5f, 0.5f);

  for (size_t g = 0; g < grainSkins.size(); g++) {
    for (size_t c = 0; c < grainSkins[g].contactPos.size(); c++) {
      double diam = 0.5 * sqrt(grainSkins[g].contactArea[c]);
      vec3r arrow = grainSkins[g].contactPos[c] - grainSkins[g].baryCenter;
      glShape::tube(grainSkins[g].baryCenter, arrow, diam);
    }
  }
}

bool readSkins(const char* filename) {
  std::ifstream file(filename);
  size_t nbg;
  file >> nbg;
  vec3i cmin(50000, 50000, 50000);
  vec3i cmax(0, 0, 0);
  for (size_t g = 0; g < nbg; g++) {
    size_t nbv;
    file >> nbv;
    GrainVoxSurf GVS;
    for (size_t v = 0; v < nbv; v++) {
      int x, y, z, c;
      file >> x >> y >> z >> c;
      GVS.voxPos.push_back(vec3i(x, y, z));
      GVS.color.push_back(c);
      if (x < cmin.x) {
        cmin.x = x;
      }
      if (y < cmin.y) {
        cmin.y = y;
      }
      if (z < cmin.z) {
        cmin.z = z;
      }
      if (x > cmax.x) {
        cmax.x = x;
      }
      if (y > cmax.y) {
        cmax.y = y;
      }
      if (z > cmax.z) {
        cmax.z = z;
      }
    }
    grainSkins.push_back(GVS);
  }

  aabb.min.set((float)cmin.x, (float)cmin.y, (float)cmin.z);
  aabb.max.set((float)cmax.x, (float)cmax.y, (float)cmax.z);
  obb.center = 0.5 * (aabb.min + aabb.max);
  obb.extent.set(0.5 * (aabb.max.x - aabb.min.x), 0.5 * (aabb.max.y - aabb.min.y), 0.5 * (aabb.max.z - aabb.min.z));

  CT.setMinMax(0, (float)grainSkins.size());
  CT.setSize((int)grainSkins.size());
  CT.setTableID(20);  // 20 = RANDOM
  CT.Rebuild();

  for (size_t g = 0; g < grainSkins.size(); g++) {
    grainSkins[g].computeBaryCenter();
    grainSkins[g].computeContacts();
  }

  return true;
}

void menu(int num) {
  switch (num) {

    case 0:
      exit(0);
      break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  glutCreateMenu(menu);  // Main menu
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {

  if (argc == 2) {
    readSkins(argv[1]);
  } else {
    std::cout << "grrrrr\n";
    exit(0);
  }

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("PBC3Dbox VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  // glutSpecialFunc(processSpecialKeys);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  // ==== Init the visualizer
  center.set(0.0, 0.0, 0.0);  // where we look at
  eye.set(1.0, 1.0, 0.0);     // from where we look
  up.set(0.0, 0.0, 1.0);      // direction (normalized)

  mouse_mode = NOTHING;
  view_angle = 45.0;
  znear = 0.01f;
  zfar = 10.0f;

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, vertices);
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, 0, normals);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_COLOR_MATERIAL);

  // Create light components
  GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
  GLfloat diffuseLight[] = {0.8f, 0.8f, 0.8f, 1.0f};
  GLfloat specularLight[] = {0.5f, 0.5f, 0.5f, 1.0f};
  GLfloat positionLight0[] = {1000000.0f, 1000000.0f, 1000000.0f, 1.0f};
  GLfloat positionLight1[] = {-1000000.0f, -1000000.0f, -1000000.0f, 1.0f};

  // Assign created components to GL_LIGHT0
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, positionLight0);

  // Assign created components to GL_LIGHT1
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT1, GL_POSITION, positionLight1);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  // ==== Enter GLUT event processing cycle
  adjust_clipping_plans();
  fit_view();
  glutMainLoop();
  return 0;
}
