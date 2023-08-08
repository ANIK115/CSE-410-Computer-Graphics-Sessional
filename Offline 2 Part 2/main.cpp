#include "transformations.h"

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

    return 0;
}
