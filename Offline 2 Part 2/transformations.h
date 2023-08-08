#include <bits/stdc++.h>

using namespace std;

#define PI acos(-1.0)*1.0

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

    void createTranslationMatrix(PointVector &translation)
    {
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                if (i == j)
                    matrix[i][j] = 1.0;
                else
                    matrix[i][j] = 0.0;
            }
            matrix[i][dimension - 1] = translation.point_vector[i];
        }
    }

    void createScaleMatrix(PointVector &scale)
    {
        createIdentityMatrix();
        for (int i = 0; i < dimension-1; i++)
        {
            for (int j = 0; j < dimension-1; j++)
            {
                if (i == j)
                    matrix[i][j] = scale.point_vector[i];
                else
                    matrix[i][j] = 0.0;
            }
        }
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

    PointVector rodriguesFormula(PointVector &unitVec, PointVector &a, double angle)
    {
        angle = angle * PI / 180.0;
        PointVector c;
        c = unitVec * cos(angle) + (a * (unitVec * a)) * (1 - cos(angle)) + cross_product(a, unitVec) * sin(angle);
        return c;
    }

    void createRotationMatrix(double angle, PointVector &a)
    {
        a.normalizePoints();
        /*create 3 point vectors for i, j and k unit vector*/
        PointVector i(1.0, 0.0, 0.0);
        PointVector j(0.0, 1.0, 0.0);
        PointVector k(0.0, 0.0, 1.0);

        PointVector c1 = rodriguesFormula(i, a, angle);
        PointVector c2 = rodriguesFormula(j, a, angle);
        PointVector c3 = rodriguesFormula(k, a, angle);

        vector<PointVector> column_vector(3);
        column_vector[0] = c1;
        column_vector[1] = c2;
        column_vector[2] = c3;

        createIdentityMatrix();

        for(int i=0; i<dimension-1; i++)
        {
            for(int j=0; j<dimension-1; j++)
            {
                matrix[j][i] = column_vector[i].point_vector[j];
            }
        }
    }

    void createViewMatrix(PointVector &eye, PointVector &look, PointVector &up)
    {
        PointVector l = look - eye;
        l.normalizePoints();
        PointVector r = cross_product(l, up);
        r.normalizePoints();
        PointVector u = cross_product(r, l);
        u.normalizePoints();

        vector<PointVector> row_vector(3);

        Matrix T(4);
        T.createIdentityMatrix();
        PointVector e = eye * -1.0;
        T.createTranslationMatrix(e);

        Matrix R(4);
        R.createIdentityMatrix();
        for(int i=0; i<dimension-1; i++)
        {
            R.matrix[0][i] = r.point_vector[i];
            R.matrix[1][i] = u.point_vector[i];
            R.matrix[2][i] = -1 * l.point_vector[i];
        }

        Matrix V = R * T;
        this->matrix = V.matrix;
    }

    void createProjectionMatrix(double fovY, double aspectRatio, double near, double far)
    {
        double fovX = fovY * aspectRatio;
        double t = near * tan((fovY/2.0)*PI/180.0);
        double r = near * tan((fovX/2.0)*PI/180.0);

        createIdentityMatrix();
        matrix[0][0] = near/r;
        matrix[1][1] = near/t;
        matrix[2][2] = -(far+near)/(far-near);
        matrix[2][3] = -(2.0*far*near)/(far-near);
        matrix[3][2] = -1.0;
        matrix[3][3] = 0.0;
    }

};
