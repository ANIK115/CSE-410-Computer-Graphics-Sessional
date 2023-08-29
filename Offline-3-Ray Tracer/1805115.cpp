#include<bits/stdc++.h>
#include <GL/glut.h>
using namespace std;
#include "bitmap_image.hpp"
#include "utils.h"
#include <windows.h>
#define pi (2*acos(0.0))


// Global variables
struct point pos;   // position of the eye
struct point l;     // look/forward direction
struct point r;     // right direction
struct point u;     // up direction
bool isAxes = true;
int rotation_angle = 0; //this angle is to rotate the object w.r.t its own axis
int countDecreaseX = 0;
int countDecreaseZ = 0;
int countIncreaseX = 0;
int countIncreaseZ = 0;

double nearDistance, farDistance, fovY, fovX, aspectRatio;
int levelOfRecursion;
int imageSize;
double screenWidth, screenHeight;

vector<Object*> objects;
vector<PointLight*> pointLights;
vector<SpotLight*> spotLights;


/*Used from demo code that was provided in moodle*/
/* Draw axes: X in Red, Y in Green and Z in Blue */
// draw axes
/*Used from demo code that was provided in moodle*/
/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes()
{
    glLineWidth(3);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0); // Red
    // X axis
    glVertex3f(100, 0, 0);
    glVertex3f(-100, 0, 0);

    glColor3f(0, 1, 0); // Green
    // Y axis
    glVertex3f(0, -100, 0);
    glVertex3f(0, 100, 0);

    glColor3f(0, 0, 1); // Blue
    // Z axis
    glVertex3f(0, 0, 100);
    glVertex3f(0, 0, -100);
    glEnd();
}

/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display()
{
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW); // To operate on Model-View matrix
    glLoadIdentity();           // Reset the model-view matrix


    // cout << "eyes: " << pos.x << ", " << pos.y << ", " << pos.z << endl;
    // cout << "look: " << l.x << ", " << l.y << ", " << l.z << endl;
    // cout << "up: " << u.x << ", " << u.y << ", " << u.z << endl;
    // cout << "right: " << r.x << ", " << r.y << ", " << r.z << endl;

    gluLookAt(pos.x,pos.y,pos.z,
              pos.x+l.x,pos.y+l.y,pos.z+l.z,
              u.x,u.y,u.z);

    // draw
    glRotated(rotation_angle, 0, 1, 0);
    if (isAxes)
        drawAxes();

    for(int i=0; i<objects.size(); i++)
    {
        objects[i]->draw();
    }

    glutSwapBuffers(); // Render now
}



void reshapeListener(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();            // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}


void cross_product(point &u, point &l, point &r)
{
    r.x=u.y*l.z - l.y*u.z; 
    r.y=l.x*u.z - u.x*l.z;
    r.z=u.x*l.y - l.x*u.y;
}
/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int xx,int yy){
    double rate = 0.01;
    double clockwise_rotation_angle = 10.0;
	switch(key){

		case '2':
			r.x = r.x*cos(rate)+l.x*sin(rate);
			r.y = r.y*cos(rate)+l.y*sin(rate);
			r.z = r.z*cos(rate)+l.z*sin(rate);

			l.x = l.x*cos(rate)-r.x*sin(rate);
			l.y = l.y*cos(rate)-r.y*sin(rate);
			l.z = l.z*cos(rate)-r.z*sin(rate);
			break;

        case '1':
			r.x = r.x*cos(-rate)+l.x*sin(-rate);
			r.y = r.y*cos(-rate)+l.y*sin(-rate);
			r.z = r.z*cos(-rate)+l.z*sin(-rate);

			l.x = l.x*cos(-rate)-r.x*sin(-rate);
			l.y = l.y*cos(-rate)-r.y*sin(-rate);
			l.z = l.z*cos(-rate)-r.z*sin(-rate);
			break;

        case '3':
			l.x = l.x*cos(rate)+u.x*sin(rate);
			l.y = l.y*cos(rate)+u.y*sin(rate);
			l.z = l.z*cos(rate)+u.z*sin(rate);

			u.x = u.x*cos(rate)-l.x*sin(rate);
			u.y = u.y*cos(rate)-l.y*sin(rate);
			u.z = u.z*cos(rate)-l.z*sin(rate);
			break;

        case '4':
			l.x = l.x*cos(-rate)+u.x*sin(-rate);
			l.y = l.y*cos(-rate)+u.y*sin(-rate);
			l.z = l.z*cos(-rate)+u.z*sin(-rate);

			u.x = u.x*cos(-rate)-l.x*sin(-rate);
			u.y = u.y*cos(-rate)-l.y*sin(-rate);
			u.z = u.z*cos(-rate)-l.z*sin(-rate);
			break;

        case '6':
			u.x = u.x*cos(rate)+r.x*sin(rate);
			u.y = u.y*cos(rate)+r.y*sin(rate);
			u.z = u.z*cos(rate)+r.z*sin(rate);

			r.x = r.x*cos(rate)-u.x*sin(rate);
			r.y = r.y*cos(rate)-u.y*sin(rate);
			r.z = r.z*cos(rate)-u.z*sin(rate);
			break;

        case '5':
			u.x = u.x*cos(-rate)+r.x*sin(-rate);
			u.y = u.y*cos(-rate)+r.y*sin(-rate);
			u.z = u.z*cos(-rate)+r.z*sin(-rate);

			r.x = r.x*cos(-rate)-u.x*sin(-rate);
			r.y = r.y*cos(-rate)-u.y*sin(-rate);
			r.z = r.z*cos(-rate)-u.z*sin(-rate);
			break;
        case 'd':

            rotation_angle += clockwise_rotation_angle;
            rotation_angle = rotation_angle% 360;
            break;

        case 'a':
            rotation_angle -= clockwise_rotation_angle;
            rotation_angle = rotation_angle% 360;          
            break;
            
            break;
		default:
			break;
	}
	glutPostRedisplay();
}



void specialKeyListener(int key, int x,int y)
{
	switch(key){
		case GLUT_KEY_UP:		//down arrow key
			pos.x+=l.x;
			pos.y+=l.y;
			pos.z+=l.z;

			break;
		case GLUT_KEY_DOWN:		// up arrow key
			pos.x-=l.x;
			pos.y-=l.y;
			pos.z-=l.z;
			break;

		case GLUT_KEY_RIGHT:
			pos.x+=r.x;
			pos.y+=r.y;
			pos.z+=r.z;
			break;
		case GLUT_KEY_LEFT :
			pos.x-=r.x;
			pos.y-=r.y;
			pos.z-=r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x+=u.x;
			pos.y+=u.y;
			pos.z+=u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
            pos.x-=u.x;
			pos.y-=u.y;
			pos.z-=u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
	glutPostRedisplay();
}


void loadData()
{
    ifstream in("scene.txt");
    in >> nearDistance >> farDistance >> fovY >> aspectRatio;
    fovX = fovY * aspectRatio;

    in >> levelOfRecursion >> imageSize;

    screenHeight = 2 * nearDistance * tan(fovY/2);
    screenWidth = 2 * nearDistance * tan(fovX/2);
    
    double floorWidth;
    in >> floorWidth;
    double floor_ambient, floor_diffuse, floor_reflection;
    in >> floor_ambient >> floor_diffuse >> floor_reflection;


    Object *floor;
    PointVector floor_ref(pos.x, 0, pos.z);
    floor = new Floor(floorWidth, floor_ref);
    floor->setColor(Color(0.5, 0.5, 0.5));
    double co_efficients[] = {floor_ambient, floor_diffuse, 0.0, floor_reflection};
    floor->setCoEfficients(co_efficients);
    floor->setShininess(1);
    objects.push_back(floor);

    int numberOfObjects;
    in >> numberOfObjects;

    for(int i=0; i<numberOfObjects; i++)
    {
        string objectType;
        in >> objectType;

        Object *object;

        if(objectType == "sphere")
        {
            object = new Sphere();
            in >> *((Sphere*)object);
        }else
        {
            cout << "Unknown object type: " << objectType << endl;
        }
        objects.push_back(object);
    }

    int numberOfNormalLights;
    in >> numberOfNormalLights;

    for(int i=0; i<numberOfNormalLights; i++)
    {
        PointLight *pointLight = new PointLight();
        in >> *pointLight;
        pointLights.push_back(pointLight);
    }

}
void initialization()
{

    pos.x = 120; pos.y = 25, pos.z = 0;
    l.x = 0-1; l.y = 0; l.z = 0;
    r.x = 0; r.y = 0; r.z = -1;
    u.x = 0; u.y = 1; u.z = 0;

    loadData();
    glClearColor(0,0,0,0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovY, aspectRatio, nearDistance, farDistance); //fovY,aspect,zNear,zFar
}

/* Initialize OpenGL Graphics */
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
    initialization();
}


/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char **argv)
{
    glutInit(&argc, argv);                                    // Initialize GLUT
    glutInitWindowSize(640, 640);                             // Set the window's initial width & height
    glutInitWindowPosition(50, 50);                           // Position the window's initial top-left corner
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color
    glutCreateWindow("Ray Tracing");                          // Create a window with the given title

    glutDisplayFunc(display);                                 // Register display callback handler for window re-paint
    // glutReshapeFunc(reshapeListener);                         // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);                       // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);                      // Register callback handler for special-key event
    initGL();                                                 // Our own OpenGL initialization
    glutMainLoop();                                           // Enter the event-processing loop

    objects.clear();
    objects.shrink_to_fit();
    return 0;
}
