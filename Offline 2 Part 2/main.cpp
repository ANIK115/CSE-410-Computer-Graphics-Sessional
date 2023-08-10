#include "transformations.h"
#include "bitmap_image.hpp"

void set_color(bitmap_image &image, int screen_width, int screen_height, int r, int g, int b)
{
    for(int i=0; i<screen_width; i++)
    {
        for(int j=0; j<screen_height; j++)
        {
            image.set_pixel(i, j, r, g, b);
        }
    }
}

int screen_width, screen_height;
double min_x, min_y, max_x, max_y;
double left_limit = -1, right_limit = 1, bottom_limit = -1, top_limit = 1;
double z_max = 1.0;
double z_min = -1.0;

void set_min_max(Triangle &t)
{
    min_x = min(t.points[0].point_vector[0], min(t.points[1].point_vector[0], t.points[2].point_vector[0]));
    min_y = min(t.points[0].point_vector[1], min(t.points[1].point_vector[1], t.points[2].point_vector[1]));
    max_x = max(t.points[0].point_vector[0], max(t.points[1].point_vector[0], t.points[2].point_vector[0]));
    max_y = max(t.points[0].point_vector[1], max(t.points[1].point_vector[1], t.points[2].point_vector[1]));
}

void clip_triangle(Triangle &t, double left_x, double right_x, double bottom_y, double top_y)
{
    min_x = max(left_x, min_x);
    min_y = max(bottom_y, min_y);
    max_x = min(right_x, max_x);
    max_y = min(top_y, max_y);
}

double clip_x(double x)
{
    x < min_x ? x = min_x : x > max_x ? x = max_x : x = x;
    return x;
}


void calculate_z_buffer(vector<Triangle> &triangles, vector<vector<double>> &z_buffer, bitmap_image &image)
{
    double dx,dy;
    dx = (right_limit - left_limit) / screen_width;
    dy = (top_limit - bottom_limit) / screen_height;

    double top_y, bottom_y, left_x, right_x;
    top_y = top_limit - (dy/2.0);
    bottom_y = bottom_limit + (dy/2.0);
    left_x = left_limit + (dx/2.0);
    right_x = right_limit - (dx/2.0);

    for(int i=0; i<triangles.size(); i++)
    {
        set_min_max(triangles[i]);
        clip_triangle(triangles[i], left_x, right_x, bottom_y, top_y);
    
        // cout << "Min x: " << min_x << endl;
        // cout << "Min y: " << min_y << endl;
        // cout << "Max x: " << max_x << endl;
        // cout << "Max y: " << max_y << endl;
        // cout << endl;

        Triangle triangle = triangles[i];
        // cout << "Triangle " << i+1 << endl;
        int start_y = round((top_y-min_y)/dy);
        int end_y = round((top_y-max_y)/dy);
        double top_scanline = top_y;

        for(int j=start_y; j>= end_y; j--)
        {
            double scan_y = top_scanline - (j * dy);
            vector<double> intersections_x(2), intersections_z(2), clipped_x(2);
            int index = 0;

            for(int k=0; k<3; k++)
            {
                /*If the two points are parallel to the scan_y line*/
                if(triangle.points[k].point_vector[1] == triangle.points[(k+1)%3].point_vector[1])
                {
                    continue;
                }

                /*If the two points are not parallel to the scan_y line*/
                double x1, y1, z1, x2, y2, z2, x, z;
                x1 = triangle.points[k].point_vector[0];
                y1 = triangle.points[k].point_vector[1];
                z1 = triangle.points[k].point_vector[2];
                x2 = triangle.points[(k+1)%3].point_vector[0];
                y2 = triangle.points[(k+1)%3].point_vector[1];
                z2 = triangle.points[(k+1)%3].point_vector[2];

                if(scan_y >= min(y1, y2) && scan_y <= max(y1, y2))
                {
                    x = x1 + ((scan_y - y1) * (x2 - x1)) / (y2 - y1);
                    z = z1 + ((scan_y - y1) * (z2 - z1)) / (y2 - y1);
                    intersections_x[index] = x;
                    intersections_z[index] = z;
                    clipped_x[index] = clip_x(x);
                    index++;
                }
            }

            intersections_z[0] = intersections_z[1]- (intersections_z[1] - intersections_z[0]) * (intersections_x[1]-clipped_x[0]) / (intersections_x[1] - intersections_x[0]); 
            intersections_z[1] = intersections_z[1]- (intersections_z[1] - intersections_z[0]) * (intersections_x[1]-clipped_x[1]) / (intersections_x[1] - intersections_x[0]); 
            double xa = clipped_x[0];
            double xb = clipped_x[1];
            double za = intersections_z[0];
            double zb = intersections_z[1];

            if(xb < xa)
            {
                swap(xa, xb);
                swap(za, zb);
            }

            int start = round((xa - left_x) / dx);
            int end = round((xb - left_x) / dx);

            for(int k=start; k<=end; k++)
            {
                double x = left_x + (k * dx);
                double z = zb- (zb - za) * (xb-x) / (xb - xa);
                if(z >= z_min && z <= z_max && z < z_buffer[j][k])
                {
                    z_buffer[j][k] = z;
                    image.set_pixel(k, j, triangle.rgb_color[0], triangle.rgb_color[1], triangle.rgb_color[2]);
                }
            }
        }

    }

    ofstream fout("z_buffer.txt");
    for(int i=0; i<screen_height; i++)
    {
        for(int j=0; j<screen_width; j++)
        {
            if(z_buffer[i][j] < z_max)
            {
                fout << setprecision(6) << fixed << z_buffer[i][j] << "\t";
            }
        }
        fout << endl;
    }
    fout.close();
    image.save_image("out.bmp");

    z_buffer.clear();
    z_buffer.shrink_to_fit();
}


int main()
{
    PointVector eye, look, up;
    double fovY, aspectRatio, near, far;

    ifstream fin("scene.txt");
    ofstream fout("stage1.txt");

    fin >> eye >> look >> up;
    fin >> fovY >> aspectRatio >> near >> far;

    Matrix M(4);
    stack<Matrix> S;
    M.createIdentityMatrix();
    S.push(M);

    int number_of_triangles = 0;

    while(true)
    {
        string command;
        fin >> command;
        if(command == "triangle")
        {
//            cout << "Inside triangle command\n";
            number_of_triangles++;
            PointVector a, b, c;
            fin >> a >> b >> c;

            a = S.top() * a;
            b = S.top() * b;
            c = S.top() * c;

            fout << a << endl;
            fout << b << endl;
            fout << c << endl;
            fout << endl;
        }else if(command == "translate")
        {
//            cout << "Inside translate command\n";
            PointVector translate;
            fin >> translate;

            Matrix T(4);
            T.createTranslationMatrix(translate);

            M = S.top();
            M = M * T;
            S.pop();
            S.push(M);
        }else if(command == "scale")
        {
//            cout << "Inside scale command\n";
            PointVector scale;
            fin >> scale;

            Matrix Sc(4);
//            cout << "Calling createScaleMatrix\n";
            Sc.createScaleMatrix(scale);

            M = S.top();
//            cout << "Calling Scale matrix multiplication\n";
            M = M * Sc;
            S.pop();
            S.push(M);
        }else if(command == "rotate")
        {
//            cout << "Inside rotate command\n";
            double angle;
            PointVector axis;
            fin >> angle >> axis;

            Matrix R(4);
            R.createRotationMatrix(angle, axis);

            M = S.top();
            M = M * R;
            S.pop();
            S.push(M);
        }else if(command == "push")
        {
            M = S.top();
            S.push(M);
//            cout << "Push command executed\n";
        }else if(command == "pop")
        {
//            cout << "Inside pop command\n";
            if(S.size() == 1)
            {
                cout << "Stack is empty" << endl;
                return 0;
            }
            S.pop();
        }else
        {
            break;
        }
    }
    fin.close();
    fout.close();


    fin.open("stage1.txt");
    fout.open("stage2.txt");

    Matrix ViewTransformMatrix(4);
    ViewTransformMatrix.createViewMatrix(eye, look, up);

    for(int i=0; i<number_of_triangles; i++)
    {
        PointVector a, b, c;
        fin >> a >> b >> c;

        a = ViewTransformMatrix * a;
        b = ViewTransformMatrix * b;
        c = ViewTransformMatrix * c;

        fout << a << endl;
        fout << b << endl;
        fout << c << endl;
        fout << endl;
    }

    fin.close();
    fout.close();

    fin.open("stage2.txt");
    fout.open("stage3.txt");

    Matrix ProjectionTransformMatrix(4);
    ProjectionTransformMatrix.createProjectionMatrix(fovY, aspectRatio, near, far);

    for(int i=0; i<number_of_triangles; i++)
    {
        PointVector a, b, c;
        fin >> a >> b >> c;

        a = ProjectionTransformMatrix * a;
        b = ProjectionTransformMatrix * b;
        c = ProjectionTransformMatrix * c;

        fout << a << endl;
        fout << b << endl;
        fout << c << endl;
        fout << endl;
    }

    fin.close();
    fout.close();

    fin.open("config.txt");
    fin >> screen_width >> screen_height;
    fin.close();

    bitmap_image image(screen_width, screen_height);
    set_color(image, screen_width, screen_height, 0,0,0);



    vector<vector<double>> z_buffer(screen_height, vector<double>(screen_width, z_max));

    fin.open("stage3.txt");
    vector<Triangle> triangles;
    for(int i=0; i<number_of_triangles; i++)
    {
        PointVector a, b, c;
        fin >> a >> b >> c;

        Triangle t(a,b,c);
        triangles.push_back(t);
    }

    calculate_z_buffer(triangles, z_buffer, image);
    fin.close();
    return 0;
}
