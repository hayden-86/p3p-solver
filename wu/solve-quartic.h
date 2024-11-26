/////////////////////
// solve-quartic.h //
/////////////////////

#ifndef solve_quartic_h
#define solve_quartic_h

#include <iostream>
#include <complex>
#include <cmath>

// const double PI = 3.141592653589793238462643383279502884;
// const double M_2PI = 2*PI;
const double M_2PI = 6.2831853071795864;
const double eps=1e-12;
const double EPSILON = 1e-30;
#define Might_high_precision 0
#define Time_tag 0
#include <wu/p3p_timers_wu.h>

unsigned int solveP3(double *x,double a,double b,double c) {

    double q  = (a*a - 3.0*b)/9.0; // q/3
	double r  = (a*(2.0*a*a-9.0*b) + 27.0*c)/54.0;  // r/2
    double r2 = r*r;
	double q3 = q*q*q;
	double A,B;

    if(r2<q3)
    {
        double q_sqrt = std::sqrt(q);
        double t=r/(q*q_sqrt);
        if( t<-1) t=-1;
        if( t> 1) t= 1;        
        t=std::acos(t);
        a/=3.0; 
        q=-2*q_sqrt;        
        x[0]=q*std::cos(t/3.0)-a;
#if Might_high_precision
        x[1]=q*std::cos((t+M_2PI)/3.0)-a;
        x[2]=q*std::cos((t-M_2PI)/3.0)-a;
#endif
        return 1;
    }   
    else
    {
        A =-pow(std::abs(r)+std::sqrt(r2-q3),1./3.0);
        if( r<0 ) A=-A;
        B = (0==A ? 0 : q/A);
        a/=3.0;
        x[0] =(A+B)-a;
#if Might_Low_precisio
        x[1] =-0.5*(A+B)-a;
        x[2] = 0.5*std::sqrt(3.)*(A-B);
        if(std::abs(x[2])<eps) { x[2]=x[1]; return 2; }
#endif
        return 1;
    }
}

void solve_quartic_inner(double a, double b, double c, double d, double retval[4],int& n_roots)
{
    double a3 = -b;
    double b3 = a * c - 4.0 * d;
    double c3 = -a * a * d - c * c + 4.0 * b * d;
    // time 19ns 

    double x3[3]={0,0,0};
    unsigned int iZeroes = solveP3(x3, a3, b3, c3);
    // time 70ns 

    double q1, q2, p1, p2, D, sqD, y;
    y = x3[0];
    D = y * y - 4.0 * d; // solve the quadratic equations for h1, h2.
#if Might_high_precision
    if (iZeroes != 1)
    {
        if (std::abs(x3[1]) > std::abs(y)) 
        {
            // static int i = 0;
            // std::cout<<"==: "<<i++<<std::endl;
            y = x3[1];
        }
        if (std::abs(x3[2]) > std::abs(y))
        {
            y = x3[2];
        }
    }
    if (std::abs(D) < EPSILON) // if D is close to 0
    {
        static int i = 0;
        std::cout<<"==: "<<i++<<std::endl;
        q1 = q2 = y * 0.5;
        
        D = a * a - 4.0 * (b - y); // solve the quadratic equations for p1,p2
        if (std::abs(D) < EPSILON)
        {
            p1 = p2 = a * 0.5;
        }
        else
        {
            sqD = sqrt(D);
            p1 = (a + sqD) * 0.5;
            p2 = (a - sqD) * 0.5;
        }
    }
    else
#endif
    {
        sqD = std::sqrt(D);
        double sqD_inv = 1/sqD;               
        q1 = (y + sqD) * 0.5;
        q2 = (y - sqD) * 0.5;
        p1 = (a * q1 - c) * sqD_inv; // Use Cramer's rule to solve for p1 and p2
        p2 = (c - a * q2) * sqD_inv;
    }
    // 83ns
    
    double D1 = p1 * p1 - 4.0 * q1; // solve the equation: x^2 + p1 * x + q1 = 0
    if (D1 >= 0.0)
    {
        sqD = std::sqrt(D1);
        retval[n_roots++] = (sqD - p1) * 0.5;
        retval[n_roots++] = (-sqD - p1) * 0.5;
    }
#if Might_high_precision
    else
    {
        if (std::abs(sqrt(-D1) * 0.5)<1e-10) 
        {
            retval[n_roots++] = -p1 * 0.5;
            retval[n_roots++] = -p1 * 0.5;
        }
    }
#endif        

    double D2 = p2 * p2 - 4.0 * q2; // solve the equation: x^2 + p2 * x + q2 = 0
    if (D2 >= 0.0)
    {
        sqD = std::sqrt(D2);
        retval[n_roots++] = (sqD - p2) * 0.5;
        retval[n_roots++] = (-sqD -p2) * 0.5;
    }
#if Might_high_precision
    else
    {
        if (std::abs(std::sqrt(-D2) * 0.5)<1e-10) 
        {
            retval[n_roots++] = -p2 * 0.5;
            retval[n_roots++] = -p2 * 0.5;
        }
    }
#endif
   
    // time 93ms 
}

void solve_quartic_Ferrari_Modified(double A3, double A2, double A1, double A0, double retval[4],int& n_roots) {
    double C = A3 / 4.0;
    double C_2 = C * C;
    double b2 = A2 - 6 * C_2;
    double b1 = A1 + (8 * C_2 - 2 * A2 ) * C;
    double b0 = A0 - A1 * C + (A2 - 3 * C_2) * C_2;

    double cubic_roots[3];
    double b2_b2_4_b0 = b2 * b2 / 4 - b0;
    unsigned int num_cubic_roots = solveP3(cubic_roots, b2, b2_b2_4_b0, -b1 * b1 / 8);
#if Might_high_precision
    double m = 0;
    for (unsigned int i = 0; i < num_cubic_roots; ++i) {
        if (cubic_roots[i] > 0) {
            m = cubic_roots[i];
            break;
        }
    }
#else
    double m = cubic_roots[0];
#endif  


    double Sigma = (b1 > 0) ? 1 : -1;
    double R = Sigma * std::sqrt((m + b2) * m + b2_b2_4_b0);
    double sqrt_m_2 = std::sqrt(m / 2.0);
    double negetive_m_b2_2 = -(m +b2) / 2.0;
    double negative_m_2_b2_2_R = negetive_m_b2_2 - R;
    double negative_m_2_b2_2_plus_R = negetive_m_b2_2 + R;
    if(negative_m_2_b2_2_R>=0)
    {
        double sqrt_negative_m_2_b2_2_R = std::sqrt(negative_m_2_b2_2_R);
        retval[n_roots++] = sqrt_m_2 - C + sqrt_negative_m_2_b2_2_R;
        retval[n_roots++] = sqrt_m_2 - C - sqrt_negative_m_2_b2_2_R;
    }
    if(negative_m_2_b2_2_plus_R>=0)
    {
        double sqrt_negative_m_2_b2_2_plus_R = std::sqrt(negative_m_2_b2_2_plus_R);
        retval[n_roots++] = -sqrt_m_2 - C + sqrt_negative_m_2_b2_2_plus_R;
        retval[n_roots++] = -sqrt_m_2 - C - sqrt_negative_m_2_b2_2_plus_R;        
    }
}

#endif // solve_quartic_h