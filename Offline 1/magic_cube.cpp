#include <windows.h> // for MS Windows
#include <GL/glut.h> // GLUT, include glu.h and gl.h
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

#define pi (2 * acos(0.0))

struct point
{
    GLfloat x, y, z;
    point() {}
    point(double x, double y, double z) : x(x), y(y), z(z) {}
    point operator +(point b)  {return point(x+b.x,y+b.y, z+b.z);}
    point operator -(point b)  {return point(x-b.x,y-b.y, z-b.z);}
    point operator *(double b)  {return point(x*b,y*b, z*b);}
    point operator *(point b)  {return point(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);}
};

// Global variables
struct point pos;   // position of the eye
struct point l;     // look/forward direction
struct point r;     // right direction
struct point u;     // up direction

/* Initialize OpenGL Graphics */
void initGL()
{
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
}

double distance(point a, point b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

point centre_of_mass(point a, point b, point c)
{
    return point((a.x+b.x+c.x)/3.0, (a.y+b.y+c.y)/3.0, (a.z+b.z+c.z)/3.0);
}

// Global variables
bool isAxes = true, isCylinder = true, isOctahedron = true, isSphere = true;
double scale = 1.0;
int columns = 0;
int rotation_angle = 0; //this angle is to rotate the object w.r.t its own axis
point triangel_A(1, 0, 0), triangel_B(0, 1, 0), triangel_C(0, 0, 1);
double triangle_side = distance(triangel_A, triangel_B);
point centre(0,0,0);
double height = triangle_side;
double common_radius = distance(centre, centre_of_mass(triangel_A, triangel_B, triangel_C));


/*Used from demo code that was provided in moodle*/
/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes()
{
    glLineWidth(3);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0); // Red
    // X axis
    glVertex3f(0, 0, 0);
    glVertex3f(2, 0, 0);

    glColor3f(0, 1, 0); // Green
    // Y axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 2, 0);

    glColor3f(0, 0, 1); // Blue
    // Z axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 2);
    glEnd();
}



void drawTriangle()
{
    glPushMatrix();

    // Translating the points towards the centre of the triangle within the same plane
    glTranslated(1 / 3.0 - 1 * scale / 3.0, 1 / 3.0 - 1 * scale / 3.0, 1 / 3.0 - 1 * scale / 3.0);

    glBegin(GL_TRIANGLES);
    glVertex3f(scale, 0, 0);
    glVertex3f(0, scale, 0);
    glVertex3f(0, 0, scale);
    glEnd();

    glPopMatrix();
}


void drawOctadehron()
{
    /*This part is for the uppder octahedron*/
    glPushMatrix();
    glColor3f(0, 1, 1);
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(1, 0.4, 1);
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(0, 1, 1);
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(1, 0.4, 1);
    drawTriangle();
    glPopMatrix();

    /*This part is for the lower octahedron*/
    glPushMatrix();
    glRotated(180, 1, 0, 0);
    glColor3f(0, 1, 1);
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(1, 0.4, 1);
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(0, 1, 1);
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(1, 0.4, 1);
    drawTriangle();
    glPopMatrix();
    glPopMatrix();
    glColor3f(1, 1, 1);
}

// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
std::vector<point> buildUnitPositiveX(int subdivision)
{
    double sphere_radius = common_radius * (1 - scale);

    const float DEG2RAD = acos(-1) / 180.0f;

    std::vector<point> vertices;
    float n1[3]; // normal of longitudinal plane rotating along Y-axis
    float n2[3]; // normal of latitudinal plane rotating along Z-axis
    float v[3];  // direction vector intersecting 2 planes, n1 x n2
    float a1;    // longitudinal angle along Y-axis
    float a2;    // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;
    columns = pointsPerRow;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for (unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for (unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float n_scale = sphere_radius / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            v[0] *= n_scale;
            v[1] *= n_scale;
            v[2] *= n_scale;

            point p;
            p.x = v[0];
            p.y = v[1];
            p.z = v[2];

            // add a vertex into array
            vertices.push_back(p);
        }
    }

    return vertices;
}

vector<vector<point>> build_2D_matrix(vector<point> mat)
{
    vector<vector<point>> matrix(columns);
    for (int i = 0; i < columns; i++)
    {
        matrix[i] = vector<point>(columns);
        for (int j = 0; j < columns; j++)
        {
            matrix[i][j] = mat[i * columns + j];
        }
    }
    return matrix;
}

void drawSphereSegment()
{
    std::vector<point> vertices = buildUnitPositiveX(5);
    vector<vector<point>> matrix = build_2D_matrix(vertices);
    glPushMatrix();
    glTranslated(scale, 0, 0);
    glBegin(GL_QUADS);
    /*Here the points were generated column wise. That means
    1   5   9   13
    2   6   10  14
    3   6   11  15
    4   8   12  16  
    
    To draw a quad, we need to take 4 points in clockwise, or anti-clockwise.
    Here I have followed clock-wise ordering. Therefore I need to take points 1, 5, 6 and 2
    */
    for (int i = 0; i < columns - 1; i++)
    {
        for (int j = 0; j < columns - 1; j++)
        {
            glVertex3d(matrix[i][j].x, matrix[i][j].y, matrix[i][j].z);
            glVertex3d(matrix[i][j + 1].x, matrix[i][j + 1].y, matrix[i][j + 1].z);
            glVertex3d(matrix[i + 1][j + 1].x, matrix[i + 1][j + 1].y, matrix[i + 1][j + 1].z);
            glVertex3d(matrix[i + 1][j].x, matrix[i + 1][j].y, matrix[i + 1][j].z);
        }
    }
    glEnd();
    glPopMatrix();
}

void drawAllSphereSegment()
{
    /*horizontal 4 parts of the sphere*/
    glColor3d(1, 0, 0);
    drawSphereSegment();
    glPushMatrix();
    glRotated(90, 0, 1, 0);
    glColor3d(0, 1, 0);
    drawSphereSegment();
    glRotated(90, 0, 1, 0);
    glColor3d(1, 0, 0);
    drawSphereSegment();
    glRotated(90, 0, 1, 0);
    glColor3d(0, 1, 0);
    drawSphereSegment();
    glPopMatrix();

    /*Upper part of the sphere*/
    glPushMatrix();
    glRotated(90, 0, 0, 1);
    glColor3d(0, 0, 1);
    drawSphereSegment();
    glPopMatrix();

    /*Lower part of the sphere*/
    glPushMatrix();
    glRotated(-90, 0, 0, 1);
    glColor3d(0, 0, 1);
    drawSphereSegment();
    glPopMatrix();
}

void drawCylinderSegment(double height, double radius, int segments)
{
    /*Here the radius is same as the sphere radius
    The sphere radius is generated as follows:
    The centre is (0,0,0) and the sphere touches the octahedron at (1/3, 1/3, 1/3)
    So radius is 1/sqrt(3)

    The height of the cylinder is the length of a side of the triangle which is sqrt(2)
    */
    radius *= (1 - scale);
    double d = radius / sin((109.47 / 2) * M_PI / 180.0);
    height = height * scale;

    struct point points[segments + 1];

    double angle = 70.5287794;
    glPushMatrix();
    glRotated(-45, 0, 1, 0);    //This rotation takes the cylinder to the edge of the octahedron in ZX plane
    glTranslated(1 / sqrt(2) - d, 0, 0);    //This shifts the centre of the cylinder into correct position
    glRotated(-angle / 2.0, 0, 0, 1);   //This takes lower half portion of the cylinder to cover the lower gap of octahedron

    /*This portion of the code is used from the demo code that was provided in moodle*/
    double tempx = radius, tempy = 0;
    double currx, curry;
    glBegin(GL_QUADS);
    for (int i = 1; i <= segments; i++)
    {
        double theta = (i * angle * M_PI / 180) / segments;
        currx = radius * cos(theta);
        curry = radius * sin(theta);

        glColor3f(1, 1, 0);
        glVertex3f(currx, curry, height / 2.0);
        glVertex3f(currx, curry, -height / 2.0);

        glVertex3f(tempx, tempy, -height / 2.0);
        glVertex3f(tempx, tempy, height / 2.0);

        tempx = currx;
        tempy = curry;
    }
    glEnd();
    glPopMatrix();
}

void drawAllCylinderSegments()
{
    int segments = 50;
    double radius = common_radius;
    double height = triangle_side;
    drawCylinderSegment(height, radius, segments);
    glPushMatrix();
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glPopMatrix();

    glPushMatrix();
    glRotated(90, 0, 0, 1);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glPopMatrix();

    glPushMatrix();
    glRotated(90, 1, 0, 0);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glRotated(90, 0, 1, 0);
    drawCylinderSegment(height, radius, segments);
    glPopMatrix();
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
    if (isSphere)
        drawAllSphereSegment();
    if (isOctahedron)
        drawOctadehron();
    if (isCylinder)
        drawAllCylinderSegments();

    glutSwapBuffers(); // Render now
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
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

void rotation_3D(point &p, point &respect, double angle)
{
    p = p*cos(angle)+(respect*p)*sin(angle);
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

        case '.':
        if(scale<1.0)
            scale+=0.01;
        break;

        case ',':
        if(scale>0.0)
            scale-=0.01;
        break;

        case 'b':
            isCylinder = !isCylinder;
            break;
        case 'n':
            isOctahedron = !isOctahedron;
            break;
        case 'm':
            isSphere = !isSphere;
            break;
        case 'r':
            pos.x=3;pos.y=3;pos.z=3;
            l.x=-1.0/sqrt(3);l.y=-1.0/sqrt(3);l.z=-1.0/sqrt(3);
            u.x=-1.0/sqrt(6);u.y=2.0/sqrt(6);u.z=-1.0/sqrt(6);
            /*r is generated with help of the cross product of l and u*/
            
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

/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char **argv)
{
    pos.x=3;pos.y=3;pos.z=3;
    l.x=-1.0/sqrt(3);l.y=-1.0/sqrt(3);l.z=-1.0/sqrt(3);
    u.x=-1.0/sqrt(6);u.y=2.0/sqrt(6);u.z=-1.0/sqrt(6);
    /*r is generated with help of the cross product of l and u*/
    cross_product(u, l, r);
    
    glutInit(&argc, argv);                                    // Initialize GLUT
    glutInitWindowSize(640, 640);                             // Set the window's initial width & height
    glutInitWindowPosition(50, 50);                           // Position the window's initial top-left corner
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color
    glutCreateWindow("OpenGL 3D Drawing");                    // Create a window with the given title
    glutDisplayFunc(display);                                 // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);                         // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);                       // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);                      // Register callback handler for special-key event
    initGL();                                                 // Our own OpenGL initialization
    glutMainLoop();                                           // Enter the event-processing loop
    return 0;
}
