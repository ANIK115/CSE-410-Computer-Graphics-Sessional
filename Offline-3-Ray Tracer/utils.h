#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include <GL/glut.h>
#include <cmath>
using namespace std;

#define PI acos(-1.0)*1.0

//create 4 Enums for different types of co-efficients and initialize them
enum
{
    AMBIENT,
    DIFFUSE,
    SPECULAR,
    REFLECTION
};

//create 4 Enums for different types of objects and initialize them
enum
{
    SPHERE,
    PYRAMID,
    CUBE,
    FLOOR
};

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

class PointVector
{

public:
    int dimension;
    vector<double> point_vector;

    PointVector()
    {
        dimension = 4;
        point_vector.resize(dimension);
        for (int i = 0; i < dimension - 1; i++)
        {
            point_vector[i] = 0.0;
        }
        point_vector[dimension - 1] = 1.0; // value of n
    }

    PointVector(double x, double y, double z)
    {
        dimension = 4;
        point_vector.resize(dimension);
        point_vector[0] = x;
        point_vector[1] = y;
        point_vector[2] = z;
        point_vector[3] = 1.0;
    }

    /*Copy Constructor*/
    PointVector(const PointVector &p)
    {
        this->dimension = p.dimension;
        this->point_vector.resize(dimension);
        for (int i = 0; i < dimension; i++)
        {
            this->point_vector[i] = p.point_vector[i];
        }
    }

    /*Operator overloading*/
    PointVector operator+(PointVector const& p)
    {
        PointVector result;
        for (int i = 0; i < dimension - 1; i++)
        {
            result.point_vector[i] = this->point_vector[i] + p.point_vector[i];
        }
        return result;
    }

    PointVector operator-(PointVector const &p)
    {
        PointVector result;
        for (int i = 0; i < dimension - 1; i++)
        {
            result.point_vector[i] = this->point_vector[i] - p.point_vector[i];
        }
        return result;
    }

    double operator*(PointVector const& p)
    {
        double result = 0.0;
        for (int i = 0; i < dimension - 1; i++)
        {
            result += this->point_vector[i] * p.point_vector[i];
        }
        return result;
    }

    PointVector operator*(double scalar)
    {
        PointVector result;
        for (int i = 0; i < dimension - 1; i++)
        {
            result.point_vector[i] = this->point_vector[i] * scalar;
        }
        return result;
    }


    void makeHomogeneous()
    {
        for(int i=0; i<dimension; i++)
        {
            point_vector[i] /= point_vector[dimension-1];
        }
    }

    double magnitude()
    {
        double sum = 0.0;
        for (int i = 0; i < dimension - 1; i++)
        {
            sum += point_vector[i] * point_vector[i];
        }
        return sqrt(sum);
    }

    void normalizePoints()
    {
        double magnitude = this->magnitude();
        for (int i = 0; i < dimension - 1; i++)
        {
            point_vector[i] /= magnitude;
        }
    }

};

istream& operator>>(istream &in, PointVector &p)
{
    for(int i=0; i<p.dimension-1; i++)
    {
        in >> p.point_vector[i];
    }
    return in;
}

ofstream& operator<<(ofstream &out, PointVector &p)
{
    out << fixed << setprecision(7);
    for(int i=0; i<p.dimension-1; i++)
    {
        out << p.point_vector[i];
        if (i != p.dimension - 2)
            out << " ";
    }
    return out;
}

PointVector cross_product(PointVector p, PointVector q)
{
    PointVector result;
    result.point_vector[0] = p.point_vector[1] * q.point_vector[2] - p.point_vector[2] * q.point_vector[1];
    result.point_vector[1] = p.point_vector[2] * q.point_vector[0] - p.point_vector[0] * q.point_vector[2];
    result.point_vector[2] = p.point_vector[0] * q.point_vector[1] - p.point_vector[1] * q.point_vector[0];
    return result;
}


class Color
{
public:
    double red, green, blue;

    Color(double red, double green, double blue)
    {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }

    Color()
    {
        red = 0.0;
        green = 0.0;
        blue = 0.0;
    }

    //override * to multiply with a scalar
    Color operator*(double scalar)
    {
        Color result;
        result.red = this->red * scalar;
        result.green = this->green * scalar;
        result.blue = this->blue * scalar;
        return result;
    }
    
};

class PointLight
{
public:
    PointVector position;
    Color color;

    PointLight(PointVector position, Color color)
    {
        this->position = position;
        this->color = color;
    }

    PointLight()
    {
        position = PointVector();
        color = Color();
    }

    void draw()
    {
        glPointSize(5);
        glBegin(GL_POINTS);
        glColor3f(color.red, color.green, color.blue);
        glVertex3f(position.point_vector[0], position.point_vector[1], position.point_vector[2]);
        glEnd();
    }

    friend istream& operator>>(istream &in, PointLight &pl)
    {
        in >> pl.position;
        in >> pl.color.red >> pl.color.green >> pl.color.blue;
        return in;
    }

};

class SpotLight
{
public:
    PointLight point_light;
    PointVector direction;
    double cut_off_angle;

    SpotLight(PointLight point_light, PointVector direction, double cut_off_angle)
    {
        this->point_light = point_light;
        this->direction = direction;
        this->cut_off_angle = cut_off_angle;
    }

    SpotLight()
    {
        point_light = PointLight();
        direction = PointVector();
        cut_off_angle = 0.0;
    }

    void draw()
    {
        //draw a cone shape
        glPushMatrix();
        glTranslatef(point_light.position.point_vector[0], point_light.position.point_vector[1], point_light.position.point_vector[2]);
        glRotatef(90, 1, 0, 0);
        glColor3f(point_light.color.red, point_light.color.green, point_light.color.blue);
        glutSolidCone(0.5, 1, 20, 20);
        glPopMatrix();
    }

    friend istream& operator>>(istream& in, SpotLight &spot_light)
    {
        in >> spot_light.point_light.position;
        in >> spot_light.point_light.color.red >> spot_light.point_light.color.green >> spot_light.point_light.color.blue;
        in >> spot_light.direction;
        in >> spot_light.cut_off_angle;
        return in;
    }

};

class Ray 
{
public:
    PointVector origin;
    PointVector direction;

    Ray(PointVector origin, PointVector direction)
    {
        this->origin = origin;
        this->direction = direction;
    }

    Ray()
    {
        origin = PointVector();
        direction = PointVector();
    }
        friend ostream& operator<<(ostream& out, Ray &ray)
    {
        out << "Origin: " << "(" << ray.origin.point_vector[0] << ", " << ray.origin.point_vector[1] << ", " << ray.origin.point_vector[2] << ")" << endl;
        out << "Direction: " << "(" << ray.direction.point_vector[0] << ", " << ray.direction.point_vector[1] << ", " << ray.direction.point_vector[2] << ")" << endl;
        return out;
    }
};

class Object;
extern vector <PointLight*> point_lights;
extern vector <Object*> objects;
extern vector <SpotLight*> spot_lights;
class Object
{
public:
    PointVector reference_point;
    double height, width, length;
    Color color;
    double co_efficients[4];
    int shininess;
    int type; // 1 for sphere, 2 for pyramid, 3 for cube, 4 for floor

    Object()
    {
        reference_point = PointVector();
        height = 0.0;
        width = 0.0;
        length = 0.0;
        color = Color();
        for (int i = 0; i < 4; i++)
        {
            co_efficients[i] = 0.0;
        }
        shininess = 0;
    }

    int getType()
    {
        return type;
    }
    
    
    void setColor(Color color)
    {
        this->color = color;
    }
    void setShininess(int shininess)
    {
        this->shininess = shininess;
    }
    void setCoEfficients(double co_efficients[])
    {
        for (int i = 0; i < 4; i++)
        {
            this->co_efficients[i] = co_efficients[i];
        }
    }
    virtual void draw() = 0;
    virtual double find_intersection(Ray ray) = 0;
    virtual Ray get_normal(PointVector intersection_point, Ray incident_ray) = 0;

    double intersect(Ray ray, Color &color, int level)
    {
        double t = find_intersection(ray);
        //No intersection
        if(t < 0)
            return -1;
        //When level is 0, the purpose of the intersect() method is to determine the nearest object only. No color computation is required then
        if(level == 0) return t;

        //level > 0 and we need to do some color computation
        PointVector intersection_point = ray.origin + ray.direction * t;

        //color of the object
        Color object_color = getColor(intersection_point);

        //ambient color
        Color ambient_color = object_color * co_efficients[AMBIENT];

        //calculating color for all the point lights
        for(int i=0; i<point_lights.size(); i++)
        {
            Ray incidentLightRay = Ray(point_lights[i]->position, intersection_point - point_lights[i]->position);
            incidentLightRay.direction.normalizePoints();

            if((intersection_point - incidentLightRay.origin).magnitude() < 1e-5)
            {
                continue;
            }

            //checking if the point is in shadow
            bool isInShadow = false;
            for(int j=0; j<objects.size(); j++)
            {
                double t = objects[j]->find_intersection(incidentLightRay);
                if(t > 0 && t+ 1e-5 < (intersection_point - incidentLightRay.origin).magnitude())
                {
                    isInShadow = true;
                    break;
                }
            }

            if(isInShadow)
                continue;
            
            //diffuse color

            Ray normal = get_normal(intersection_point, incidentLightRay);

            double lambert_value = max(0.0, normal.direction * incidentLightRay.direction);

            return 0;

        }

        
    }

    virtual Color getColor(PointVector intersection_point)
    {
        return color;
    }

};

class Floor : public Object
{
public:
    int number_of_tiles;
    int changeStartX, changeFinishX, changeStartZ, changeFinishZ;
    PointVector new_reference_point;
    Floor()
    {
        number_of_tiles = 100;
        type = FLOOR;
    }

    Floor(int shellWidth, PointVector reference_point)
    {
        this->reference_point = reference_point;
        this->width = shellWidth;
        number_of_tiles = 100;
        type = FLOOR;
        changeStartX = changeStartZ = changeFinishX = changeFinishZ = 0;
    }

    void decreaseChangeX()
    {
        changeStartX -= width;
        changeFinishX -= width;
    }
    void increaseChangeX()
    {
        changeStartX += width;
        changeFinishX += width;
    }

    void decreaseChangeZ()
    {
        changeStartZ -= width;
        changeFinishZ -= width;
    }

    void increaseChangeZ()
    {
        changeStartZ += width;
        changeFinishZ += width;
    }




    virtual void draw()
    {
        int startX = reference_point.point_vector[0]-1000+changeStartX;
        int finishX = reference_point.point_vector[0]+1000+changeFinishX;
        number_of_tiles = (finishX - startX) / width;
        // int startY = reference_point.point_vector[1]-1000;
        // int finishY = reference_point.point_vector[1]+1000;
        int startZ = reference_point.point_vector[2]-1000+changeStartZ;
        int finishZ = reference_point.point_vector[2]+1000+changeFinishZ;
        for(int i=0; i<number_of_tiles; i++)
        {
            for(int j=0; j<number_of_tiles; j++)
            {
                if((i+j)%2 == 0)
                    glColor3f(1.0, 1.0, 1.0);
                else
                    glColor3f(0.0, 0.0, 0.0);
                
                glBegin(GL_QUADS);
                {
                    //draw in XY plane
                    glVertex3f(startX + i*width, reference_point.point_vector[1], startZ + j*width);

                    glVertex3f(startX + i*width, reference_point.point_vector[1], startZ + (j+1)*width);

                    glVertex3f(startX + (i+1)*width, reference_point.point_vector[1], startZ + (j+1)*width);

                    glVertex3f(startX + (i+1)*width, reference_point.point_vector[1], startZ + j*width);
                }
                glEnd();
                
            }
        }
    }

    virtual double find_intersection(Ray ray)
    {
        double t = -1.0 * ray.origin.point_vector[1] / ray.direction.point_vector[1];
        return t;
    }

    virtual Ray get_normal(PointVector intersection_point, Ray incident_ray)
    {
        PointVector normal = PointVector(0, 1, 0);
        normal.normalizePoints();
        Ray normal_ray = Ray(intersection_point, normal);
        return normal_ray;
    }

    virtual double intersect(Ray ray, Color &color, int level)
    {
        double t = find_intersection(ray);
        if (t > 0)
        {
            if (level == 0)
            {
                color = this->color;
            }
            return t;
        }
        return -1;
    }

    virtual Color getColor(PointVector intersection_point)
    {
        int x = (intersection_point.point_vector[0] - reference_point.point_vector[0]) / width;
        int z = (intersection_point.point_vector[2] - reference_point.point_vector[2]) / width;
        if ((x + z) % 2 == 0)
            return Color(1.0, 1.0, 1.0);
        else
            return Color(0.0, 0.0, 0.0);
    }
};


class Sphere : public Object 
{
public:
    Sphere()
    {
        reference_point = PointVector();
        length = 50.0;
        type = SPHERE;
    }

    Sphere(PointVector center, double radius)
    {
        reference_point = center;
        length = radius;
        type = SPHERE;
    }

    void drawSphere(double radius, int stacks, int slices) 
    {
        struct point points[stacks+1][slices+1];
        for (int j = 0; j <= stacks; j++) {
            double phi = -M_PI / 2.0 + j * M_PI / stacks;
            double r = radius * cos(phi);
            double h = radius * sin(phi);
            for (int i = 0; i < slices+1; i++) {
                double theta = i * 2.0 * M_PI / slices;
                points[j][i].x = r * cos(theta);
                points[j][i].y = r * sin(theta);
                points[j][i].z = h;
            }
        }

        glPushMatrix();
        glTranslatef(reference_point.point_vector[0], reference_point.point_vector[1], reference_point.point_vector[2]);
        glColor3f(color.red, color.green, color.blue);

        glBegin(GL_QUADS);
            for (int j = 0; j < stacks; j++) {
                for (int i = 0; i < slices; i++) {
                    GLfloat c = (2+cos((i+j) * 2.0 * M_PI / slices)) / 3;
                    glVertex3f(points[j][i].x, points[j][i].y, points[j][i].z);
                    glVertex3f(points[j][i+1].x, points[j][i+1].y, points[j][i+1].z);

                    glVertex3f(points[j+1][i+1].x, points[j+1][i+1].y, points[j+1][i+1].z);
                    glVertex3f(points[j+1][i].x, points[j+1][i].y, points[j+1][i].z);
                }
            }
        glEnd();
        glPopMatrix();
    }

    virtual void draw()
    {
        drawSphere(length, 15, 50);
    }

    virtual double find_intersection(Ray ray)
    {
        //formula used from lecture slide and schaums book
        ray.direction.normalizePoints();    //just to be sure
        double a = 1.0;
        double b = 2.0 * (ray.direction * (ray.origin - reference_point));
        double c = (ray.origin - reference_point) * (ray.origin - reference_point) - length * length;
        double d = b * b - 4 * a * c;
        if (d < 0)
            return -1;
        else
        {
            double t1 = (-b + sqrt(d)) / (2.0 * a);
            double t2 = (-b - sqrt(d)) / (2.0 * a);
            if (t1 < 0 && t2 < 0)
                return -1;
            else if (t1 < 0)
                return t2;
            else if (t2 < 0)
                return t1;
            else
                return min(t1, t2);
        }
    }

    virtual Ray get_normal(PointVector intersection_point, Ray incident_ray)
    {
        return Ray();
    }

    friend std::istream& operator>>(std::istream& in, Sphere& s)
    {
    in >> s.reference_point >> s.length; // center and radius
    in >> s.color.red >> s.color.green >> s.color.blue; // color
    for(int i = 0; i < 4; i++) in >> s.co_efficients[i];
    in >> s.shininess;
    return in;
    }

};
