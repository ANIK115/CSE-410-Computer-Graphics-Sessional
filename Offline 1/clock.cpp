#include <iostream>
#include <GL/glut.h>
#include <iostream>
#include <ctime>
#include <cmath>

/*declaring global variables*/

double clockX, clockY;
double clockOuterRadius = 0.36;
double clockInnerRadius = 0.3;
double hourAngle, minuteAngle, secondAngle;
float hourHand = 2.5, minuteHand = 1.5, secondHand = 1.2;
float hourWidth = 3 , minuteWidth = 2.5, secondWidth = 1.5;
float pendulumRadius = 0.05;
float maxPendulumAngle = 30;
float pendulumAngle;
float pendulumTime;
double length_of_pendulum = 0.3;

//STRUCTURE DEFINING POINT
typedef struct Point {
    double x , y;
}point;


void drawClock(double radius)
{
    glBegin(GL_LINE_LOOP);
    glColor3f(0.5f, 0.5f, 1.0f); // Light-blue

    clockX = 0;
    clockY = 0.3;

    for (float theta = 0; theta < 360; theta += 10)
    {
        float x = clockX + radius * cos(theta / 180 * M_PI);
        float y = clockY + radius * sin(theta / 180 * M_PI);
        glVertex2f(x, y);
    }
    glEnd();
}

void drawClockBox()
{
    glLineWidth(5);
    glPushMatrix();
    glTranslated(clockX, clockY, 0);
    glBegin(GL_LINE_LOOP);
    glColor3f(1,1,1);
    glVertex2d(-clockOuterRadius-0.1, clockOuterRadius+0.1);
    glVertex2d(-clockOuterRadius-0.1, -clockOuterRadius-0.1);
    glVertex2d(0, -(clockY+clockOuterRadius+pendulumRadius+0.1));
    glVertex2d(clockOuterRadius+0.1, -clockOuterRadius-0.1);
    glVertex2d(clockOuterRadius+0.1, clockOuterRadius+0.1);
    
    glEnd();
    glPopMatrix();
    glLineWidth(1);

}

void drawClockCentre()
{
    glPointSize(10);
    glBegin(GL_POINTS);
    glColor3f(1,1,1);
    clockX = 0;
    clockY = 0.3;
    glVertex2d(clockX, clockY);
    glEnd();
    glPointSize(1);
}

void drawPoint(point p) {
	glBegin(GL_POINTS);
		glVertex2f(p.x , p.y);
	glEnd();
}

void drawLine(point p1, point p2)
{
    glBegin(GL_LINES);
        glVertex2f(p1.x, p1.y);
        glVertex2f(p2.x, p2.y);
    glEnd();
}

void drawMarks(void) {
	point p1,p2,p3;
	int count = 0;
    glMatrixMode(GL_MODELVIEW);     // To operate on Model-View matrix
    glLoadIdentity();

    p1.x = 0;
    p1.y = clockInnerRadius ;

    p2.x = 0;
    p2.y = clockInnerRadius-0.06 ;

    p3.x = 0;
    p3.y = clockInnerRadius-0.03 ;

	for(double i=0 ; i<= 360 ; i+=6) {
        glPushMatrix();
        glTranslatef(0, 0.3, 0);
        glRotatef(i,0.0,0.0,1.0);
            if(count % 5 == 0)  //for the hour ticks
            {
                glLineWidth(2);
                drawLine(p1,p2);
            }else{            //for the minute ticks    
                glLineWidth(1);
                drawLine(p1,p3);
            }
        glPopMatrix();
        count++;
	}
    glLineWidth(1);
}

/*This function takes fraction which provides a length with respect to the clock radius and width is used for
clock hand width */
void drawClockHand(float fraction, float width)
{
    glLineWidth(width);
    glBegin(GL_LINES);

        glVertex2d(0,0);
        glVertex2d(0,clockInnerRadius/fraction);

    glEnd();
    glLineWidth(1);
}

void drawClockHands()
{
    time_t now = time(0); // get current date and time  
    tm* ltm = localtime(&now);  
    float sec = ltm->tm_sec;
    float min = ltm->tm_min;
    int hour = ltm->tm_hour;
    if(hour > 12)
    {
        hour = hour % 12;
    }
    //Here everything is converted to seconds first and then to degrees
    secondAngle = sec * 6;
    minuteAngle = (min*60 + sec) * 360/3600.0;  //the minute hand will complete 360 degree in 60 minutes or 3600 seconds
    hourAngle = (hour * 3600 + min * 60 + sec) * 360 / (3600*12.0); //the hour hand will complete 360 degree in 12 hours


    glMatrixMode(GL_MODELVIEW);     // To operate on Model-View matrix
    glLoadIdentity();

    glPushMatrix();

    glTranslatef(clockX, clockY, 0);    //Translating the origin to the centre of the clock
    glPushMatrix();
    glRotated(-secondAngle, 0, 0, 1);
    drawClockHand(secondHand, secondWidth);
    glPopMatrix();

    glPushMatrix();
    glRotated(-minuteAngle, 0, 0, 1);
    drawClockHand(minuteHand, minuteWidth);
    glPopMatrix();

    glPushMatrix();
    glRotated(-hourAngle, 0, 0, 1);
    drawClockHand(hourHand, hourWidth);
    glPopMatrix();

    glPopMatrix();
}

void drawPendulum()
{
    pendulumAngle = maxPendulumAngle * cos(M_PI*pendulumTime);
    pendulumTime += 0.02;
    glMatrixMode(GL_MODELVIEW);     // To operate on Model-View matrix
    glLoadIdentity();

    glColor3f(1,1,1);
    glPushMatrix();
    glTranslatef(clockX, clockY, 0);             //To shift the origin to the centre of the clock 

    glTranslatef(0, -clockOuterRadius, 0);  //To shift the pendulum below the clock
    glRotatef(pendulumAngle, 0, 0, 1);      //Rotating the pendulum

    glBegin(GL_POLYGON);
        for (float theta = 0; theta < 360; theta += 10)
        {
            float x = pendulumRadius/3 * cos(theta / 180 * M_PI);
            float y = pendulumRadius/3 * sin(theta / 180 * M_PI);
            glVertex2f(x, y);
        }
    glEnd();

    glBegin(GL_LINES);
        glVertex2f(0,0);
        glVertex2f(0,-length_of_pendulum);
    glEnd();


    //drawing bob
    glTranslatef(0, -length_of_pendulum, 0);    //shifting the bob to the end of the string
    glBegin(GL_POLYGON);
    glColor3d(1,1,1);
    glRotatef(pendulumAngle, 0, 0, 1);
        for (float theta = 0; theta < 360; theta += 10)
        {
            float x = pendulumRadius * cos(theta / 180 * M_PI);
            float y = pendulumRadius * sin(theta / 180 * M_PI);
            glVertex2f(x, y);
        }

    glEnd();
    glPopMatrix();
}


void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    drawClockCentre();
    drawClockBox();
    drawClock(clockInnerRadius);
    drawClock(clockOuterRadius);
    drawMarks();
    drawClockHands();
    drawPendulum();
    glFlush();
}

void timerListener(int value)
{
    glutPostRedisplay();
    glutTimerFunc(1000, timerListener, 0);
}

void timerPendulum(int value)
{
    glutPostRedisplay();
    glutTimerFunc(10, timerPendulum, 0);
}


int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(100, 50);
    glutCreateWindow("Analog Clock");
    glutDisplayFunc(display);
    glutTimerFunc(0, timerListener, 0);
    glutTimerFunc(0, timerPendulum, 0);
    glutMainLoop();
    return 0;
}