#include <bits/stdc++.h>
using namespace std;

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

    PointVector(const PointVector &p)
    {
        this->dimension = p.dimension;
        for (int i = 0; i < dimension; i++)
        {
            this->point_vector[i] = p.point_vector[i];
        }
    }

    /*Operator overloading*/
    PointVector operator+(PointVector p)
    {
        PointVector result;
        for (int i = 0; i < dimension - 1; i++)
        {
            result.point_vector[i] = this->point_vector[i] + p.point_vector[i];
        }
        return result;
    }

    PointVector operator-(PointVector p)
    {
        PointVector result;
        for (int i = 0; i < dimension - 1; i++)
        {
            result.point_vector[i] = this->point_vector[i] - p.point_vector[i];
        }
        return result;
    }

    PointVector operator*(PointVector p)
    {
        PointVector result;
        for (int i = 0; i < dimension - 1; i++)
        {
            result.point_vector[i] = this->point_vector[i] * p.point_vector[i];
        }
        return result;
    }

    void makeHomogeneous()
    {
        
    }
};

class Matrix
{

public:
    int dimension;
    vector<vector<double>> matrix;

    Matrix(int dim)
    {
        this->dimension = dim;
        matrix.resize(dimension);
        for (int i = 0; i < dimension; i++)
            matrix[i].resize(dimension);
    }

    void createTranslationMatrix(PointVector translation)
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

    void createScaleMatrix(PointVector scale)
    {
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                if (i == j)
                    matrix[i][j] = scale.point_vector[i];
                else
                    matrix[i][j] = 0.0;
            }
        }
    }

    /*Overloading * for matrix multiplication*/
    Matrix operator*(Matrix B)
    {
        Matrix result(this->dimension);
        for (int i = 0; i < dimension; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < dimension; j++)
            {
                for (int k = 0; k < dimension; k++)
                {
                    sum += this->matrix[i][k] * B.matrix[k][j];
                }
                result.matrix[i][j] = sum;
            }
        }
    }

    /*Overloading * for n*n matrix and n point-vector*/
    PointVector operator*(PointVector p)
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


    }
};