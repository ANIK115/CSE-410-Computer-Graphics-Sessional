#include <GL/glut.h>
#include <cmath>

// Input values
float distances[2] = {1, 1000};
float fovY = 80;
float aspectRatio = 1;
float eye[3] = {0, 0, distances[0]}; // Assuming looking down the positive z-axis
float lookAt[3] = {0, 2, 0}; // Looking at point (0, 2, 0)
float up[3] = {0, 1, 0};

int checkerboardSize = 1; // Width of each checkerboard cell

void drawCheckerboard() {
    int numSquares = 10; // Number of squares to draw

    for (int i = -numSquares; i <= numSquares; ++i) {
        for (int j = -numSquares; j <= numSquares; ++j) {
            glColor3f(((i + j) % 2 == 0) ? 1.0f : 0.0f, ((i + j) % 2 == 0) ? 1.0f : 0.0f, ((i + j) % 2 == 0) ? 1.0f : 0.0f);

            glBegin(GL_QUADS);
            glVertex3f(i * checkerboardSize, 0, j * checkerboardSize);
            glVertex3f((i + 1) * checkerboardSize, 0, j * checkerboardSize);
            glVertex3f((i + 1) * checkerboardSize, 0, (j + 1) * checkerboardSize);
            glVertex3f(i * checkerboardSize, 0, (j + 1) * checkerboardSize);
            glEnd();
        }
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set up perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovY, aspectRatio, distances[0], distances[1]);

    // Set up camera transformation
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2], lookAt[0], lookAt[1], lookAt[2], up[0], up[1], up[2]);

    drawCheckerboard();

    glFlush();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(768, 768);
    glutCreateWindow("Camera and Checkerboard");
    glEnable(GL_DEPTH_TEST);
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}
