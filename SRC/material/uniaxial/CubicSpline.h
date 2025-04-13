#pragma once
class CubicSpline
{
    
public:
    CubicSpline();
    ~CubicSpline();

    void Fit(double* x, int xl, double* y, int yl);
    double Eval(double x);
    double EvalT(double x);
    /*
    double* FitAndEval(double* x, int xLength, double* y,int yLength, double* xs, int xsLength);
    void Fit(double* x, int xLength, double* y,int yLength);
    int GetNextXIndex(double x);
    double* Eval(double* x, int xLength);
    void Fit(double* x,int xLength, double* y,int yLength, double a0, double b0, double an, double bn);
    void Fit3(double* x,int xLength, double* y,int yLength, double a0, double b0, double an, double bn);
    void Fit4(double* x,int xLength, double* y,int yLength, double a0, double b0, double an, double bn);
    void Affiche();
    */

private:
    // N-1 spline coefficients for N points
    double* xs,* ys;
    int xsL,ysL;
    double* c1s,* c2s,* c3s;
    int c1sc,c2sc,c3sc,c1sL,c2sL,c3sL;
    int sc1sL,sc2sL,sc3sL;
    /*
    int xOrigLength;
    int yOrigLength;
    int _lastIndex;
    double* a;
    double* b;
    double* xOrig;
    double* yOrig;*/


};

