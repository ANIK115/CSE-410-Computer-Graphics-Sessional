#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include <GL/glut.h>
#include <cmath>
#include <math.h>
using namespace std;

#define PI acos(-1.0)*1.0
#define M_PI		3.14159265358979323846

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

    PointVector operator-()
    {
        PointVector result;
        for (int i = 0; i < dimension - 1; i++)
        {
            result.point_vector[i] = -this->point_vector[i];
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


class Matrix
{

public:
    int dimension;
    vector<vector<double>> matrix;

    Matrix()
    {
        dimension = 4;
        matrix.resize(dimension);
        for (int i = 0; i < dimension; i++)
            matrix[i].resize(dimension);
    }

    Matrix(int dim)
    {
        this->dimension = dim;
        matrix.resize(dimension);
        for (int i = 0; i < dimension; i++)
            matrix[i].resize(dimension);
    }

    void createIdentityMatrix()
    {
        for (int i = 0; i < dimension; i++)
        {
            matrix[i][i] = 1.0;
        }
    }

    /*Overloading * for matrix multiplication*/
    Matrix operator*(Matrix const &B)
    {
        Matrix result(this->dimension);
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < dimension; k++)
                {
                    sum += this->matrix[i][k] * B.matrix[k][j];
                }
                result.matrix[i][j] = sum;
            }
        }
        return result;
    }

    /*Overloading * for n*n matrix and n point-vector*/
    PointVector operator*(PointVector const &p)
    {
        PointVector transformed_point;
        for(int i=0; i<dimension; i++)
        {
            transformed_point.point_vector[i] = 0;
            for(int j=0; j<dimension; j++)
            {
                transformed_point.point_vector[i] += matrix[i][j] * p.point_vector[j];
            }
        }
        transformed_point.makeHomogeneous();
        return transformed_point;
    }
    double determinant()
    {

        double det = 0.0;
        for(int i=0; i<dimension; i++)
        {
            double v1 = 1.0, v2 = 1.0;
            for(int j=0; j<dimension; j++)
            {
                v1 *= matrix[(i+j)%dimension][j];
                v2 *= matrix[(i+j)%dimension][dimension-j-1];
            }
            det += v1 - v2;
        }
        return det;

    }
};


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

    //override * to multiply with another color
    Color operator*(Color const &c)
    {
        Color result;
        result.red = this->red * c.red;
        result.green = this->green * c.green;
        result.blue = this->blue * c.blue;
        return result;
    }

    Color operator+(Color const &c)
    {
        Color result;
        result.red = this->red + c.red;
        result.green = this->green + c.green;
        result.blue = this->blue + c.blue;
        return result;
    }
    
};

class PointLight
{
public:
    PointVector position;
    Color color;
    double falloff;

    PointLight(PointVector position, Color color)
    {
        this->position = position;
        this->color = color;
    }

    PointLight()
    {
        position = PointVector();
        color = Color(1.0, 1.0, 1.0);
    }

    void draw()
    {
        glPointSize(30);
        glBegin(GL_POINTS);
        glColor3f(color.red, color.green, color.blue);
        glVertex3f(position.point_vector[0], position.point_vector[1], position.point_vector[2]);
        glEnd();
    }

    friend istream& operator>>(istream &in, PointLight &pl)
    {
        in >> pl.position;
        in >> pl.falloff;
        pl.color = Color(1.0, 1.0, 1.0);
        return in;
    }

};

class SpotLight
{
public:
    PointLight point_light;
    PointVector direction;
    double cut_off_angle;
    double falloff;
    Color color;

    SpotLight(PointLight point_light, PointVector direction, double cut_off_angle)
    {
        this->point_light = point_light;
        this->direction = direction;
        this->cut_off_angle = cut_off_angle;
        this->color = point_light.color;
    }

    SpotLight()
    {
        point_light = PointLight();
        direction = PointVector();
        cut_off_angle = 0.0;
        color = Color(1.0, 1.0, 1.0);
    }

    void drawCone(double height, double radius, int segments, Color color) 
    {
        double tempx = radius, tempy = 0;
        double currx, curry;
        glBegin(GL_TRIANGLES);
            for (int i = 1; i <= segments; i++) {
                double theta = i * 2.0 * M_PI / segments;
                currx = radius * cos(theta);
                curry = radius * sin(theta);

                GLfloat c = (2+cos(theta))/3;
                glColor3f(color.red , color.green, color.blue);
                glVertex3f(0, 0, height/2);
                glVertex3f(currx, curry, -height/2);
                glVertex3f(tempx, tempy, -height/2);

                tempx = currx;
                tempy = curry;
            }
        glEnd();
    }

    void draw()
    {
        //draw a cone shape
        glPushMatrix();
        glTranslatef(point_light.position.point_vector[0], point_light.position.point_vector[1], point_light.position.point_vector[2]);

        PointVector dir = direction-point_light.position;
        dir.normalizePoints();

        double angleX = acos(dir.point_vector[0]);
        double angleY = acos(dir.point_vector[1]);
        double angleZ = acos(dir.point_vector[2]);

      
        glRotatef(angleY*180.0/M_PI, 0.0, 1.0, 0.0);  
        glRotatef(-angleX*180.0/M_PI, 1.0, 0.0, 0.0);
        glRotatef(-angleZ*180.0/M_PI, 0.0, 0.0, 1.0);

        glColor3f(color.red, color.green, color.blue);
        //set size of the solid cone
        drawCone(25, 10, 20, color);
        glPopMatrix();
    }

    friend istream& operator>>(istream& in, SpotLight &spot_light)
    {
        in >> spot_light.point_light.position;
        in >> spot_light.falloff;
        in >> spot_light.direction;
        in >> spot_light.cut_off_angle;
        spot_light.color = Color(1.0, 1.0, 1.0);
        spot_light.direction = spot_light.direction - spot_light.point_light.position;
        //convert to radian
        // spot_light.cut_off_angle = spot_light.cut_off_angle * M_PI / 180.0;
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
        this->direction.normalizePoints();
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
extern int depthOfRecursion;
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
        ray.direction.normalizePoints();
        double t = find_intersection(ray);
        //No intersection
        if(t < 0)
            return -1;
        //When level is 0, the purpose of the intersect() method is to determine the nearest object only. No color computation is required then
        if(level == 0) 
        {
            // color = getColor(ray.origin + ray.direction * t);
            return t;
        }

        //level > 0 and we need to do some color computation
        PointVector intersection_point = ray.origin + ray.direction * t;

        //color of the object
        Color object_color = getColor(intersection_point);


        //ambient color
        color = object_color * co_efficients[AMBIENT];
        // cout << "Ambient Color: " << color.red << " " << color.green << " " << color.blue << endl;



        //calculating color for all the point lights
        double diffuse_color =0, specular_color=0;
        for(int i=0; i<point_lights.size(); i++)
        {

            //scaling factor = exp(-falloff * distance^2)
            double distance = (intersection_point - ray.origin).magnitude();
            double scaling_factor = exp(-(point_lights[0]->falloff * distance * distance));
            Ray incidentLightRay = Ray(point_lights[i]->position, intersection_point - point_lights[i]->position);
            incidentLightRay.direction.normalizePoints();

            if((intersection_point - incidentLightRay.origin).magnitude() <= 0)
            {
                continue;
            }

            //check if the ray and the light ray are in the opposite direction
            // Ray lineRay = Ray(ray.origin, point_lights[i]->position - ray.origin);

            Ray lineRay = Ray(ray.origin, point_lights[i]->position - ray.origin);

            double intersection_t = find_intersection(lineRay);

            if(intersection_t > 0)
            { 
                double camera_lightSource_distance = (point_lights[i]->position - ray.origin).magnitude();
                PointVector camera_object_intersection = ray.origin + ray.direction * intersection_t;
                double camera_object_distance = (camera_object_intersection - ray.origin).magnitude();
                if(camera_lightSource_distance > camera_object_distance)
                    continue;
            }

            //checking if the point is in shadow
            bool isInShadow = false;
            for(int j=0; j<objects.size(); j++)
            {
                double t1 = objects[j]->find_intersection(incidentLightRay);
                if(t1 > 0 && t1+ 1e-5 < (intersection_point - incidentLightRay.origin).magnitude())
                {
                    isInShadow = true;
                    break;
                }
            }

            if(isInShadow)
                continue;
            
            //diffuse color

            Ray normal = get_normal(intersection_point, incidentLightRay);

            normal.direction.normalizePoints();
            incidentLightRay.direction.normalizePoints();

            if(this->type == FLOOR)
            {
                normal.direction = normal.direction * -1.0;
            }
            double lambert_value = max(0.0, -normal.direction * incidentLightRay.direction);

            //reflected ray
            Ray reflectedRay = Ray(intersection_point, incidentLightRay.direction - normal.direction * 2.0 * (incidentLightRay.direction * normal.direction));

            reflectedRay.direction.normalizePoints();

            //diffuse color
            // Color diffuse_color = object_color * co_efficients[DIFFUSE] * point_lights[i]->color * lambert_value * scaling_factor;
            diffuse_color = lambert_value * scaling_factor;

            double phong_value = max(0.0, -reflectedRay.direction * ray.direction);
            
            phong_value = pow(phong_value, shininess);

            // color = color + object_color* (point_lights[i]->color * co_efficients[2] * phong_value * scaling_factor);
            specular_color = phong_value * scaling_factor; 
            color = color + object_color * point_lights[i]->color * (diffuse_color * co_efficients[DIFFUSE] + specular_color * co_efficients[SPECULAR]);

        }

        diffuse_color = 0;
        specular_color = 0;

        // //calculating color for all the spot lights
        for(int i=0; i<spot_lights.size(); i++)
        {
            //scaling factor = exp(-falloff * distance^2)
            double distance = (intersection_point - ray.origin).magnitude();
            double scaling_factor = exp(-spot_lights[i]->falloff * distance * distance);
            spot_lights[i]->direction.normalizePoints();
            Ray incidentLightRay = Ray(spot_lights[i]->point_light.position, intersection_point - spot_lights[i]->point_light.position);

            incidentLightRay.direction.normalizePoints();

            double angle = acos(incidentLightRay.direction * spot_lights[i]->direction);

            //convert to degrees
            angle = angle * 180.0 / M_PI;


            if(fabs(angle) < spot_lights[i]->cut_off_angle) {

                // cout << "Angle: " << fabs(angle) << endl;

                Ray normal = get_normal(intersection_point, incidentLightRay);
                normal.direction.normalizePoints();

                if((intersection_point - incidentLightRay.origin).magnitude() < 1e-5)
                {
                    continue;
                }

                //check if the ray and the light ray are in the opposite direction
                Ray lineRay = Ray(ray.origin, point_lights[i]->position - ray.origin);

                double intersection_t = find_intersection(lineRay);

                if(intersection_t > 0)
                { 
                    double camera_lightSource_distance = (point_lights[i]->position - ray.origin).magnitude();
                    PointVector camera_object_intersection = ray.origin + ray.direction * intersection_t;
                    double camera_object_distance = (camera_object_intersection - ray.origin).magnitude();
                    if(camera_lightSource_distance > camera_object_distance)
                        continue;
                }

                bool isInShadow = false;
                for(int j=0; j<objects.size(); j++)
                {
                    double t1 = objects[j]->find_intersection(incidentLightRay);
                    if(t1 > 0 && t1+ 1e-5 < (intersection_point - incidentLightRay.origin).magnitude())
                    {
                        isInShadow = true;
                        break;
                    }
                }

                if(isInShadow)
                    continue;

                // cout << "Calculating color for spot light: " << i << endl;
                //diffuse color
                if(this->type == FLOOR)
                {
                    normal.direction = normal.direction * -1.0;
                }
                double lambert_value = max(0.0, -normal.direction * incidentLightRay.direction); 
                // cout << "Lambert Value: " << lambert_value << endl;


                diffuse_color = lambert_value * scaling_factor;

                //reflected ray
                Ray reflectedRay = Ray(intersection_point, incidentLightRay.direction - normal.direction * 2.0 * (incidentLightRay.direction * normal.direction));

                reflectedRay.direction.normalizePoints();

                //specular color
                double phong_value = max(0.0, -reflectedRay.direction * ray.direction);
                // cout << "Phong Value: " << phong_value << endl;
                if(shininess == 0)
                    phong_value = phong_value;
                else
                    phong_value = pow(phong_value, shininess);

                // color = color + object_color * (spot_lights[i]->color * co_efficients[SPECULAR] * phong_value * scaling_factor);
                specular_color = phong_value * scaling_factor;
                // cout << "Color: " << color.red << " " << color.green << " " << color.blue << endl;
                color = color + object_color * spot_lights[i]->color * (diffuse_color * co_efficients[DIFFUSE] + specular_color * co_efficients[SPECULAR]);
            }
               
        }

        // recursive reflection

        if(level > 0)
        {
            Ray reflectedRay = Ray(intersection_point, ray.direction - get_normal(intersection_point, ray).direction * 2.0 * (ray.direction * get_normal(intersection_point, ray).direction));
            reflectedRay.direction.normalizePoints();
            
            reflectedRay.origin = reflectedRay.origin + reflectedRay.direction * 1e-5;

            int nearest = -1;
            double min_t = 1e9, t=-1;
            for(int i=0; i<objects.size(); i++)
            {
                t = objects[i]->find_intersection(reflectedRay);
                if(t > 0 && t < min_t)
                {
                    min_t = t;
                    nearest = i;
                }
            }

            if(nearest != -1)
            {
                Color reflectedColor(0, 0, 0);
                double t = objects[nearest]->intersect(reflectedRay, reflectedColor, level-1);
                color = color + reflectedColor * co_efficients[REFLECTION];
            }
        }

        return t;
    }

    virtual Color getColor(PointVector intersection_point)
    {
        return color;
    }

};

class Floor : public Object
{
public:
    int number_of_tiles, number_of_tiles_in_X, number_of_tiles_in_Z, number_of_tiles_in_Y;
    int startX, finishX, startZ, finishZ, diffX, diffZ, startY, finishY, diffY;
    PointVector new_reference_point;
    Floor()
    {
        number_of_tiles = 100;
        type = FLOOR;
        number_of_tiles_in_X = number_of_tiles_in_Z = number_of_tiles_in_Y = 100;
        diffX = diffZ = diffY = 0;
    }

    Floor(int shellWidth, PointVector reference_point)
    {
        this->reference_point = reference_point;
        this->width = shellWidth;
        number_of_tiles = number_of_tiles_in_X = number_of_tiles_in_Z = number_of_tiles_in_Y = 100;
        type = FLOOR;
        diffX = diffZ = diffY = 0;
    }

    void setReferencePoint(PointVector reference_point)
    {
        this->reference_point = reference_point;
    }

    void setNewReferencePoint(PointVector new_reference_point)
    {
        this->new_reference_point = new_reference_point;
    }


    virtual void draw()
    {
        number_of_tiles = 60;
        diffX = (new_reference_point.point_vector[0] - reference_point.point_vector[0]) / width;
        diffY = (new_reference_point.point_vector[1] - reference_point.point_vector[1]) / width;

        startX = reference_point.point_vector[0]-(number_of_tiles/2)*width ;
        
        startY = reference_point.point_vector[1]-(number_of_tiles/2)*width ;

        if(diffX %2 == 0)
        {
            startX += diffX * width;
        }else
        {
            startX += (diffX-1) * width;
        }

        if(diffY %2 == 0)
        {
            startY += diffY * width;
        }else
        {
            startY += (diffY-1) * width;
        }

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
                    glVertex3f(startX + i*width,startY + j*width,  reference_point.point_vector[2]);

                    glVertex3f(startX + i*width, startY + (j+1)*width, reference_point.point_vector[2]);

                    glVertex3f(startX + (i+1)*width, startY + (j+1)*width, reference_point.point_vector[2]);

                    glVertex3f(startX + (i+1)*width, startY + j*width, reference_point.point_vector[2]);
                }
                glEnd();
                
            }
        }
    }

    virtual double find_intersection(Ray ray)
    {
        PointVector normal = PointVector(0.0, 0.0, 1.0);
        ray.direction.normalizePoints();
        double dotProduct = normal * ray.direction;

        if(round(dotProduct * 100) == 0)
            return -1;
        double t = -(normal * ray.origin) / dotProduct;    

        PointVector intersection_point = ray.origin + ray.direction * t + ray.direction * 2e-5;

        if(intersection_point.point_vector[0] <= startX || intersection_point.point_vector[0] >= startX + number_of_tiles * width)
            return -1;
        if(intersection_point.point_vector[1] <= startY || intersection_point.point_vector[1] >= startY + number_of_tiles * width)
            return -1;    
        return t;
    }

    virtual Ray get_normal(PointVector intersection_point, Ray incident_ray)
    {
        if(incident_ray.direction.point_vector[2] > 0)
        {
            return Ray(intersection_point, PointVector(0.0, 0.0, 1.0));
        }
        else
        {
            return Ray(intersection_point, PointVector(0.0, 0.0, -1.0));
        }
    }


    virtual Color getColor(PointVector intersection_point)
    {
        int x = (intersection_point.point_vector[0]- startX) / width;
        int y = (intersection_point.point_vector[1]- startY) / width;
        
        if((x+y)%2 == 0)
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
        if (d <= 0)
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
        return Ray(intersection_point, intersection_point - reference_point);
    }

    friend std::istream& operator>>(std::istream& in, Sphere& s)
    {
        in >> s.reference_point >> s.length; // center and radius
        in >> s.color.red >> s.color.green >> s.color.blue; // color
        in >> s.co_efficients[AMBIENT] >> s.co_efficients[DIFFUSE] >> s.co_efficients[SPECULAR] >> s.co_efficients[REFLECTION]; // co-efficients
        in >> s.shininess;
        return in;
    }

    virtual Color getColor(PointVector intersection_point)
    {
        return color;
    }

};

class Triangle : public Object
{
public:
    PointVector a, b, c;

    Triangle()
    {
        a = PointVector();
        b = PointVector();
        c = PointVector();
        type = PYRAMID;
    }

    Triangle(PointVector a, PointVector b, PointVector c, string type)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        type == "pyramid" ? this->type = PYRAMID : this->type = CUBE;
    }

    void drawTriangle()
    {
        glColor3f(color.red, color.green, color.blue);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.point_vector[0], a.point_vector[1], a.point_vector[2]);
            glVertex3f(b.point_vector[0], b.point_vector[1], b.point_vector[2]);
            glVertex3f(c.point_vector[0], c.point_vector[1], c.point_vector[2]);
        }
        glEnd();
    }

    virtual void draw()
    {
        drawTriangle();
    }

    virtual Ray get_normal(PointVector intersection_point, Ray incident_ray)
    {
        PointVector normal = cross_product(b - a, c - a);
        normal.normalizePoints();
        
        incident_ray.direction * normal < 0 ? normal = normal * -1 : normal = normal;
        return Ray(intersection_point, normal);
    }

    virtual double find_intersection(Ray ray)
    {
        ray.direction.normalizePoints();
        Matrix beta(3);
        Matrix gamma(3);
        Matrix t(3);
        Matrix A(3);

        A.matrix[0][0] = a.point_vector[0] - b.point_vector[0];
        A.matrix[0][1] = a.point_vector[0] - c.point_vector[0];
        A.matrix[0][2] = ray.direction.point_vector[0];
        A.matrix[1][0] = a.point_vector[1] - b.point_vector[1];
        A.matrix[1][1] = a.point_vector[1] - c.point_vector[1];
        A.matrix[1][2] = ray.direction.point_vector[1];
        A.matrix[2][0] = a.point_vector[2] - b.point_vector[2];
        A.matrix[2][1] = a.point_vector[2] - c.point_vector[2];
        A.matrix[2][2] = ray.direction.point_vector[2];


        beta.matrix[0][0] = a.point_vector[0] - ray.origin.point_vector[0];
        beta.matrix[0][1] = a.point_vector[0] - c.point_vector[0];
        beta.matrix[0][2] = ray.direction.point_vector[0];
        beta.matrix[1][0] = a.point_vector[1] - ray.origin.point_vector[1];
        beta.matrix[1][1] = a.point_vector[1] - c.point_vector[1];
        beta.matrix[1][2] = ray.direction.point_vector[1];
        beta.matrix[2][0] = a.point_vector[2] - ray.origin.point_vector[2];
        beta.matrix[2][1] = a.point_vector[2] - c.point_vector[2];
        beta.matrix[2][2] = ray.direction.point_vector[2];


        gamma.matrix[0][0] = a.point_vector[0] - b.point_vector[0];
        gamma.matrix[0][1] = a.point_vector[0] - ray.origin.point_vector[0];
        gamma.matrix[0][2] = ray.direction.point_vector[0];
        gamma.matrix[1][0] = a.point_vector[1] - b.point_vector[1];
        gamma.matrix[1][1] = a.point_vector[1] - ray.origin.point_vector[1];
        gamma.matrix[1][2] = ray.direction.point_vector[1];
        gamma.matrix[2][0] = a.point_vector[2] - b.point_vector[2];
        gamma.matrix[2][1] = a.point_vector[2] - ray.origin.point_vector[2];
        gamma.matrix[2][2] = ray.direction.point_vector[2];

        t.matrix[0][0] = a.point_vector[0] - b.point_vector[0];
        t.matrix[0][1] = a.point_vector[0] - c.point_vector[0];
        t.matrix[0][2] = a.point_vector[0] - ray.origin.point_vector[0];
        t.matrix[1][0] = a.point_vector[1] - b.point_vector[1];
        t.matrix[1][1] = a.point_vector[1] - c.point_vector[1];
        t.matrix[1][2] = a.point_vector[1] - ray.origin.point_vector[1];
        t.matrix[2][0] = a.point_vector[2] - b.point_vector[2];
        t.matrix[2][1] = a.point_vector[2] - c.point_vector[2];
        t.matrix[2][2] = a.point_vector[2] - ray.origin.point_vector[2];

        double detA = A.determinant();
        double betaValue = beta.determinant() / detA;
        double gammaValue = gamma.determinant() / detA;
        double tValue = t.determinant() / detA;

        if (betaValue > 0 && gammaValue > 0 && betaValue + gammaValue < 1 && tValue > 0)
            return tValue;
        else
            return -1;
    }

};

class OptimizedSquare : public Object
{
public:
    PointVector a, b, c, d;
    OptimizedSquare()
    {
        a = PointVector();
        b = PointVector();
        c = PointVector();
        d = PointVector();
        type = CUBE;
    }


    OptimizedSquare(PointVector a, PointVector b, PointVector c, PointVector d, string type)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        type == "pyramid" ? this->type = PYRAMID : this->type = CUBE;
    }

    void virtual draw()
    {
        glColor3f(color.red, color.green, color.blue);
        glBegin(GL_QUADS);
        {
            glVertex3f(a.point_vector[0], a.point_vector[1], a.point_vector[2]);
            glVertex3f(b.point_vector[0], b.point_vector[1], b.point_vector[2]);
            glVertex3f(c.point_vector[0], c.point_vector[1], c.point_vector[2]);
            glVertex3f(d.point_vector[0], d.point_vector[1], d.point_vector[2]);
        }
        glEnd();
    }

    virtual Ray get_normal(PointVector intersection_point, Ray incident_ray)
    {
        PointVector normal = cross_product(b - a, c - a);
        normal.normalizePoints();
        
        // incident_ray.direction * normal < 0 ? normal = normal * -1 : normal = normal;
        return Ray(intersection_point, normal);
    }

    virtual double find_intersection(Ray ray)
    {
        PointVector normal = cross_product(b - a, c - a);
        normal.normalizePoints();

        double dotProduct = normal * ray.direction;
        if(fabs(dotProduct) < 1e-5)
            return -1;
        
        double t = PointVector(a - ray.origin) * normal / dotProduct;
        if(t < 0)
            return -1;
        
        PointVector intersection_point = ray.origin + ray.direction * t;

        PointVector ab = b - a;
        PointVector ad = d - a;

        PointVector ap = intersection_point - a;

        double dot_ab_ap = ab * ap;
        double dot_ab_ab = ab * ab;
        double dot_ad_ap = ad * ap;
        double dot_ad_ad = ad * ad;

        if(dot_ab_ap < 0 || dot_ab_ap > dot_ab_ab || dot_ad_ap < 0 || dot_ad_ap > dot_ad_ad)
            return -1;
        
        return t;
    }


};

class Pyramid : public Object
{
public:
    Triangle triangles[4];
    OptimizedSquare square;


    Pyramid()
    {
        triangles[0] = Triangle();
        triangles[1] = Triangle();
        triangles[2] = Triangle();
        triangles[3] = Triangle();
        square = OptimizedSquare();
        type = PYRAMID;
    }

    Pyramid(PointVector lowestPoint, double width, double height)
    {
        //The base of the pyramid is square and it is in XY plane.
        //The base is centered at basePoint
        this->width = width;
        this->height = height;
        type = PYRAMID;

        //calculate the corner points of the Pyramid
        PointVector basePoint = lowestPoint + PointVector(width / 2.0, width / 2.0, 0.0);
        PointVector a = PointVector(basePoint.point_vector[0] - width / 2.0, basePoint.point_vector[1] - width / 2.0, basePoint.point_vector[2]);
        PointVector b = PointVector(basePoint.point_vector[0] + width / 2.0, basePoint.point_vector[1] - width / 2.0, basePoint.point_vector[2]);
        PointVector c = PointVector(basePoint.point_vector[0] + width / 2.0, basePoint.point_vector[1] + width / 2.0, basePoint.point_vector[2]);
        PointVector d = PointVector(basePoint.point_vector[0] - width / 2.0, basePoint.point_vector[1] + width / 2.0, basePoint.point_vector[2]);

        PointVector e = PointVector(basePoint.point_vector[0], basePoint.point_vector[1], basePoint.point_vector[2] + height);

        triangles[0] = Triangle(a, b, e, "pyramid");
        triangles[1] = Triangle(b, c, e, "pyramid");
        triangles[2] = Triangle(c, d, e, "pyramid");
        triangles[3] = Triangle(d, a, e, "pyramid");
        square = OptimizedSquare(a, b, c, d, "pyramid");
    }

    void setColor(Color color)
    {
        for (int i = 0; i < 4; i++)
        {
            triangles[i].setColor(color);
        }
        square.setColor(color);
        this->color = color;
    }

    void drawPyramid()
    {
        for (int i = 0; i < 4; i++)
        {
            triangles[i].drawTriangle();
        }
        // square.drawSquare();
    }

    virtual void draw()
    {
        drawPyramid();
    }

    virtual Ray get_normal(PointVector intersection_point, Ray ray)
    {
        PointVector normal1 = triangles[0].get_normal(intersection_point, ray).direction;
        PointVector normal2 = triangles[1].get_normal(intersection_point, ray).direction;
        PointVector normal3 = triangles[2].get_normal(intersection_point, ray).direction;
        PointVector normal4 = triangles[3].get_normal(intersection_point, ray).direction;
        // PointVector normal5 = square.get_normal(intersection_point, ray).direction;

        if (abs((intersection_point - triangles[0].a) * normal1) < 1e-5)
        {
            ray.direction * normal1 < 0 ? normal1 = normal1 * -1 : normal1 = normal1;
            return Ray(intersection_point, normal1);
        }
        else if (abs((intersection_point - triangles[1].a) * normal2) < 1e-5)
        {
            ray.direction * normal2 < 0 ? normal2 = normal2 * -1 : normal2 = normal2;
            return Ray(intersection_point, normal2);
        }
        else if (abs((intersection_point - triangles[2].a) * normal3) < 1e-5)
        {
            ray.direction * normal3 < 0 ? normal3 = normal3 * -1 : normal3 = normal3;
            return Ray(intersection_point, normal3);
        }
        else if (abs((intersection_point - triangles[3].a) * normal4) < 1e-5)
        {
            ray.direction * normal4 < 0 ? normal4 = normal4 * -1 : normal4 = normal4;
            return Ray(intersection_point, normal4);
        }
        else
        {
            // ray.direction * normal5 < 0 ? normal5 = normal5 * -1 : normal5 = normal5;
            // return Ray(intersection_point, normal5);
            return Ray();
        }
    }

    virtual double find_intersection(Ray ray)
    {
        double tValues[5];
        tValues[0] = triangles[0].find_intersection(ray);
        tValues[1] = triangles[1].find_intersection(ray);
        tValues[2] = triangles[2].find_intersection(ray);
        tValues[3] = triangles[3].find_intersection(ray);
        // tValues[4] = square.find_intersection(ray);

        double minT = -1;
        for (int i = 0; i < 5; i++)
        {
            if (tValues[i] > 0)
            {
                if (minT < 0)
                    minT = tValues[i];
                else
                    minT = min(minT, tValues[i]);
            }
        }

        minT < 0 ? minT = -1 : minT = minT;
        return minT;
    }

    virtual Color getColor(PointVector intersection_point)
    {
        return this->color;
    }

};


class Cube : public Object
{
public:
    PointVector points[8];
    PointVector bottom_lower_left;
    Triangle triangles[12];
    OptimizedSquare squares[6];

    Cube()
    {

    }

    Cube(PointVector bottom_lower_left, double side)
    {
        this->length = side;
        this->bottom_lower_left = bottom_lower_left;
        type = CUBE;
        points[0] = bottom_lower_left;
        points[1] = PointVector(bottom_lower_left.point_vector[0] + side, bottom_lower_left.point_vector[1], bottom_lower_left.point_vector[2]);
        points[2] = PointVector(bottom_lower_left.point_vector[0] + side, bottom_lower_left.point_vector[1], bottom_lower_left.point_vector[2] + side);
        points[3] = PointVector(bottom_lower_left.point_vector[0], bottom_lower_left.point_vector[1], bottom_lower_left.point_vector[2] + side);   
        points[4] = PointVector(bottom_lower_left.point_vector[0], bottom_lower_left.point_vector[1] + side, bottom_lower_left.point_vector[2]);
        points[5] = PointVector(bottom_lower_left.point_vector[0] + side, bottom_lower_left.point_vector[1] + side, bottom_lower_left.point_vector[2]);
        points[6] = PointVector(bottom_lower_left.point_vector[0] + side, bottom_lower_left.point_vector[1] + side, bottom_lower_left.point_vector[2] + side);
        points[7] = PointVector(bottom_lower_left.point_vector[0], bottom_lower_left.point_vector[1] + side, bottom_lower_left.point_vector[2] + side);
    }

    void createSquares()
    {
        squares[0] = OptimizedSquare(points[0], points[1], points[2], points[3], "cube");
        squares[1] = OptimizedSquare(points[0], points[1], points[5], points[4], "cube");
        squares[2] = OptimizedSquare(points[1], points[2], points[6], points[5], "cube");
        squares[3] = OptimizedSquare(points[2], points[3], points[7], points[6], "cube");
        squares[4] = OptimizedSquare(points[3], points[0], points[4], points[7], "cube");
        squares[5] = OptimizedSquare(points[4], points[5], points[6], points[7], "cube");
    }

    void createTriangles()
    {
        for(int i=0; i<12; i++)
        {
            triangles[i] = Triangle();
            triangles[i].setColor(color);
            triangles[i].type = CUBE;
            triangles[i].shininess = shininess;
            triangles[i].setCoEfficients(co_efficients);
        }

        triangles[0].a = points[0];
        triangles[0].b = points[1];
        triangles[0].c = points[2];

        triangles[1].a = points[0];
        triangles[1].b = points[2];
        triangles[1].c = points[3];

        triangles[2].a = points[0];
        triangles[2].b = points[1];
        triangles[2].c = points[4];

        triangles[3].a = points[1];
        triangles[3].b = points[4];
        triangles[3].c = points[5];

        triangles[4].a = points[1];
        triangles[4].b = points[2];
        triangles[4].c = points[5];

        triangles[5].a = points[2];
        triangles[5].b = points[5];
        triangles[5].c = points[6];

        triangles[6].a = points[3];
        triangles[6].b = points[2];
        triangles[6].c = points[6];

        triangles[7].a = points[3];
        triangles[7].b = points[6];
        triangles[7].c = points[7];

        triangles[8].a = points[0];
        triangles[8].b = points[3];
        triangles[8].c = points[7];

        triangles[9].a = points[0];
        triangles[9].b = points[4];
        triangles[9].c = points[7];

        triangles[10].a = points[4];
        triangles[10].b = points[5];
        triangles[10].c = points[6];

        triangles[11].a = points[4];
        triangles[11].b = points[6];
        triangles[11].c = points[7];
    }

    virtual void draw()
    {

    }
    virtual double find_intersection(Ray ray)
    {
        return -1;
    }

    virtual Ray get_normal(PointVector intersection_point, Ray incident_ray)
    {
        return Ray();
    }
};