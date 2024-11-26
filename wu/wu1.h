#pragma once
// #include <lambdatwist/p3p_timers.h>
#include <wu/p3p_timers_wu.h>
#include "solve-quartic.h"
#include <utils/cvl/matrix.h>
#include <cmath>
#include <iostream>
#include <complex>
#include <wu/refine_wu.h>
#include <iomanip>

#include"refine_wu.h"
#define Time_tag 0

using namespace std;

namespace cvl{

double find_xk_inner(double a, double b, double c) {
    if (a > 0){
        double delta = b * b - 4 * a * c;
        if (delta > 0) {
            return std::ceil((-b + std::sqrt(delta)) / (2 * a));
        } else if (delta == 0) {
            return std::ceil(-b / (2 * a) + 1);
        } else {
            return 0;
        }
    }
    else if (a < 0){
        double delta = b * b - 4 * a * c;
        if (delta > 0) {
            return (-b / (2 * a));
        }
        else if (delta == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    else{
        if (b > 0) {
            return std::ceil(c / b + 1);
        } else if (b < 0) {
            return std::floor(c / b - 1);
        } else {
            if (c > 0) {
                return 0;
            } else {
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
    }     
}


template<class T, int refinement_iterations=5>
int p3p_wu1( Vector3<T> y1,
            Vector3<T> y2,
            Vector3<T> y3,
            Vector3<T> x1,
            Vector3<T> x2,
            Vector3<T> x3,
            Vector<cvl::Matrix<T,3,3>,4>& Rs,
            Vector<Vector3<T>,4>& Ts
            ){
#if Time_tag
    TIC(lt1);
#endif
    // 1 real, 2 complex
    // y1 = Vector3<T>(-0.1477, 0.2224, 0.9637);
    // y2 = Vector3<T>(0.6786, -0.0315, 0.7338);
    // y3 = Vector3<T>(-0.5371, 0.3909, 0.7474);
    // x1 = Vector3<T>(10.4916, 6.0984, 0.6732);
    // x2 = Vector3<T>(67.7863, 1.4434, 91.6460);
    // x3 = Vector3<T>(56.3792, 31.6665, -27.3255);
    
    // // 3 real
    // y1 = Vector3<T>(-60.9753, 62.1944, 75.1468);
    // y2 = Vector3<T>(-25.2913, 30.4318, 91.2678);
    // y3 = Vector3<T>(-16.0361, 75.9026, 83.7972);
    // x1 = Vector3<T>(-36.008, 71.9073, 83.2485);
    // x2 = Vector3<T>(-24.0047, 89.9166, 37.7108);
    // x3 = Vector3<T>(-68.2862, 78.5985, 48.6707);

    y1.normalize();
    y2.normalize();
    y3.normalize();  

    T m12 = y1.dot(y2);
    T m13 = y1.dot(y3);
    T m23 = y2.dot(y3);

    Vector3<T> d12 = x1 - x2;
    Vector3<T> d13 = x1 - x3;
    Vector3<T> d23 = x2 - x3;
    Vector3<T> d12xd13(d12.cross(d13));
#if Time_tag
    TOC(lt1);
#endif
#if Time_tag
    TIC(lt2);
#endif
    T a12 = d12.squaredNorm();
    T a13 = d13.squaredNorm();
    T a23 = d23.squaredNorm();
    T a = a12/a23;
    T b = a13/a23;
 
    T b1 = -2.0*m12;
    T c1 = 1.0-a;
    T e1 = 2.0*a*m23;
    T f2 = 1.0-b;
    int valid = 0;
    T p1_0,p1_1;
    T p2_0,p2_1;
    T p3_0,p3_1;

#if 1
    p3_0 = -std::sqrt(a);
    p1_0 = -p3_0;
    p1_1 = 0;
    p2_0 = p1_0;

    T c1_inv = 1/c1;
    p2_1 = -(b1 * p1_0 + e1)*c1_inv;

    if(std::abs(p2_1)<0.05 || a==1)
    {
        double discriminant = b1 * b1 - 4 * c1; // (<0,ellipse; =0,parabola; >0,hyperbola)
        if (discriminant<0)
        {
            T cc = 4 * c1 * a + e1 * e1;
            T sqrt_delta1 = std::sqrt(cc);
            // p1_0 = 0;
            p2_0 = 0;
            T c1_ = 0.5*c1_inv;
            // p1_1 = c1_*(-e1+sqrt_delta1);
            p2_1 = c1_*(-e1-sqrt_delta1);
        }
        else if(discriminant>0)
        {
            T xk = p3_0 - 1;
            T bb1 = b1*xk+e1;
            T sqrt_delta1 = std::sqrt(bb1*bb1-4*c1*(xk*xk-a));
            T c1_ = 0.5*c1_inv;
            // p1_0 = xk;
            // p1_1 = c1_*(-bb1+sqrt_delta1);
            p2_0 = xk;
            p2_1 = c1_*(-bb1-sqrt_delta1);
        }
        else
        {
            T xk = find_xk_inner(-4 * c1 + b1 * b1, 2 * b1 * e1, 4 * c1 * a + e1 * e1);
            T bb1 = b1*xk+e1;
            T sqrt_delta1 = std::sqrt(bb1*bb1-4*c1*(xk*xk-a));
            T c1_ = 0.5*c1_inv;
            p1_0 = xk;
            p1_1 = c1_*(-bb1+sqrt_delta1);
            p2_0 = xk;
            p2_1 = c1_*(-bb1-sqrt_delta1);
        }
    }
#else
    p3_0 = 0.5*(-std::sqrt(-4 * f1));
    p3_1 = 0.0;
    if (std::abs(c1) > 1e-10)
    {
        T cc = -4 * c1 * f1 + e1 * e1;
        if (cc > 0) {
            T sqrt_delta1 = std::sqrt(cc);
            p1_0 = 0;
            T c1_ = 0.5/c1;
            p1_1 = c1_*(-e1+sqrt_delta1);
            p2_0 = 0;
            p2_1 = c1_*(-e1-sqrt_delta1);
        } else {
            T xk = find_xk_inner(-4 * c1 + b1 * b1, 2 * b1 * e1, cc);
            T bb1 = b1*xk+e1;
            T sqrt_delta1 = std::sqrt(bb1*bb1-4*c1*(xk*xk+f1));
            T c1_ = 0.5/c1;
            p1_0 = xk;
            p1_1 = c1_*(-bb1+sqrt_delta1);
            p2_0 = xk;
            p2_1 = c1_*(-bb1-sqrt_delta1);
        }
    }
#endif

    T b1_2 = 0.5*b1;
    T e1_2 = 0.5*e1;
    T aa1 = p1_0+b1_2*p1_1;
    T aa2 = b1_2*p1_0+c1*p1_1+e1_2;
    T aa3 = e1_2*p1_1-a;
    T bb1 = p2_0+b1_2*p2_1;
    T bb2 = b1_2*p2_0+c1*p2_1+e1_2;
    T bb3 = e1_2*p2_1-a;

    T P00 = aa2*bb3-aa3*bb2;
    T P01 = p1_0;
    T P02 = p2_0;
    T P10 = aa3*bb1-aa1*bb3;
    T P11 = p1_1;
    T P12 = p2_1;
    T P20 = aa1*bb2-aa2*bb1;

    T Pdet_inv;
    if (p1_0==0 && p2_0==0)
    {
        Pdet_inv = 1 / (P00 * (P11 - P12));
    }
    else
    {    
        Pdet_inv = 1 / (P00 * (P11 - P12)
              -P01 * (P10 - P12 * P20)
              +P02 * (P10 - P11 * P20));
    }
    T P_inv00 = (P11 - P12);
    T P_inv02 = (P01 * P12 - P02 * P11);
    T P_inv10 = (P12 * P20 - P10);
    T P_inv12 = (P02 * P10 - P00 * P12);
    T P_inv20 = (P10 - P11 * P20);
    T P_inv22 = (P00 * P11 - P01 * P10);

    T lambda0 = (P_inv00*p3_0+P_inv02) * Pdet_inv;
    T lambda1 = (P_inv10*p3_0+P_inv12) * Pdet_inv;
    T lambda2 = (P_inv20*p3_0+P_inv22) * Pdet_inv;

    T H00 = lambda0 * P00;
    T H01 = lambda1 * P01;
    T H02 = lambda2 * P02;
    T H10 = lambda0 * P10;
    T H11 = lambda1 * P11;
    T H12 = lambda2 * P12;
    T H20 = lambda0 * P20;
    T H21 = lambda1;
    T H22 = lambda2;

    T e2_2 = b*m23;
    T C2_00 = H00-H20*m13;
    T C2_01 = -H10*b+H20*e2_2;
    T C2_02 = H10*e2_2+H20*f2-H00*m13;
    T C2_10 = H01-H21*m13; 
    T C2_11 = -H11*b+H21*e2_2;
    T C2_12 = H11*e2_2+H21*f2-H01*m13;
    T C2_20 = H02-H22*m13; 
    T C2_21 = -H12*b+H22*e2_2;
    T C2_22 = H12*e2_2+H22*f2-H02*m13;

#if 1
    T coefficients[5];
    if (p1_0==0 && p2_0==0)
    {
        coefficients[0] = C2_11*H11+C2_12*H21;
        coefficients[1] = (C2_01*H11+C2_02*H21)+(C2_10*H00+C2_11*H10+C2_12*H20);
        coefficients[2] = (C2_00*H00+C2_01*H10+C2_02*H20)+(C2_11*H12+C2_12*H22)+(C2_21*H11+C2_22*H21);
        coefficients[3] = (C2_01*H12+C2_02*H22)+(C2_20*H00+C2_21*H10+C2_22*H20);
        coefficients[4] = C2_21*H12+C2_22*H22;
    }
    else
    {
        coefficients[0] = C2_10*H01+C2_11*H11+C2_12*H21;
        coefficients[1] = (C2_00*H01+C2_01*H11+C2_02*H21)+(C2_10*H00+C2_11*H10+C2_12*H20);
        coefficients[2] = (C2_00*H00+C2_01*H10+C2_02*H20)+(C2_10*H02+C2_11*H12+C2_12*H22)+(C2_20*H01+C2_21*H11+C2_22*H21);
        coefficients[3] = (C2_00*H02+C2_01*H12+C2_02*H22)+(C2_20*H00+C2_21*H10+C2_22*H20);
        coefficients[4] = C2_20*H02+C2_21*H12+C2_22*H22;
    }
#else
    T coefficients[5] = {
        C2_10*H01+C2_11*H11+C2_12*H21,
        (C2_00*H01+C2_01*H11+C2_02*H21)+(C2_10*H00+C2_11*H10+C2_12*H20),
        (C2_00*H00+C2_01*H10+C2_02*H20)+(C2_10*H02+C2_11*H12+C2_12*H22)+(C2_20*H01+C2_21*H11+C2_22*H21),
        (C2_00*H02+C2_01*H12+C2_02*H22)+(C2_20*H00+C2_21*H10+C2_22*H20),
        (C2_20*H02+C2_21*H12+C2_22*H22)
        };
#endif

    T roots[4]= {0, 0, 0, 0};
    int n_roots=0;
#if 0
    solve_quartic(coefficients, roots, n_roots); 
#else
    const double a_ = coefficients[0];
    const double a_inv_ = 1 / a_;
    const double b_ = coefficients[1] * a_inv_;
    const double c_ = coefficients[2] * a_inv_;
    const double d_ = coefficients[3] * a_inv_;
    const double e_ = coefficients[4] * a_inv_;
#if 1
    if(std::abs(a_)<10000){
        solve_quartic_inner(b_, c_, d_, e_, roots,n_roots);
        if(n_roots==0)
        {
            solve_quartic_Ferrari_Modified(b_, c_, d_, e_, roots,n_roots);
        }
    }
    else{      
    	solve_quartic_Ferrari_Modified(b_, c_, d_, e_, roots,n_roots);
    }
#else
    if(a_>-50000){
	    solve_quartic_inner(b_, c_, d_, e_, roots,n_roots);
        if(n_roots==0)
        {
            solve_quartic_Ferrari_Modified(b_, c_, d_, e_, roots,n_roots);
        }
    }
    else{
    	solve_quartic_Ferrari_Modified(b_, c_, d_, e_, roots,n_roots);
    }
#endif

#endif

#if 0
    Matrix<T,3,3> X(d12(0),d13(0),d12xd13(0),
                    d12(1),d13(1),d12xd13(1),
                    d12(2),d13(2),d12xd13(2));
    X=X.inverse();
#else
    double det_inv = 1 / (d12(0) * (d13(1) * d12xd13(2) - d12xd13(1) * d13(2)) -
                 d13(0) * (d12(1) * d12xd13(2) - d12xd13(1) * d12(2)) +
                 d12xd13(0) * (d12(1) * d13(2) - d13(1) * d12(2)));
    Matrix<T,3,3> X((d13(1) * d12xd13(2) - d12xd13(1) * d13(2)) * det_inv,
                (d12xd13(0) * d13(2) - d13(0) * d12xd13(2)) * det_inv,
                (d13(0) * d12xd13(1) - d12xd13(0) * d13(1)) * det_inv,
                (d12xd13(1) * d12(2) - d12(1) * d12xd13(2)) * det_inv,
                (d12(0) * d12xd13(2) - d12xd13(0) * d12(2)) * det_inv,
                (d12xd13(0) * d12(1) - d12(0) * d12xd13(1)) * det_inv,
                (d12(1) * d13(2) - d13(1) * d12(2)) * det_inv,
                (d13(0) * d12(2) - d12(0) * d13(2)) * det_inv,
                (d12(0) * d13(1) - d13(0) * d12(1)) * det_inv);
#endif

    Vector3<T> yd1;
    Vector3<T> yd2;
    Vector3<T> yd1xd2;     
  
    for(int i=0;i<n_roots;i++)
    {
        T tmp_sol0 = roots[i]; 
        T tmp_sol1 = tmp_sol0*tmp_sol0; 
        T sol0 = H00*tmp_sol0+H01*tmp_sol1+H02;
        T sol1 = H10*tmp_sol0+H11*tmp_sol1+H12;

#if 1
        if(sol0*sol1<=0)
        {
            continue;
        }
        T sol2 = (H20*tmp_sol0+H21*tmp_sol1+H22);
        if(sol0 * sol2>0)
        {
            T denom = (sol0 * sol0 + sol1 * sol1-2 * m12 * sol0 * sol1);
            if(denom<0)
                continue;
            // T sqrt_d3_square = std::sqrt(a12 * (sol2*sol2) / denom);
            // T d1 = sol0 * sqrt_d3_square/sol2;
            // T d2 = sol1 * sqrt_d3_square/sol2;
            T sqrt_d3_square;
            T d1 = 0;
            T d2 = 0;
            T sqrt_d3_square_tmp = std::sqrt(a12 / denom);
            if (sol2<0)
            {
                sqrt_d3_square = -sol2 * sqrt_d3_square_tmp;
                d1 = -sol0 * sqrt_d3_square_tmp;
                d2 = -sol1 * sqrt_d3_square_tmp;
            }
            else
            {
                sqrt_d3_square = sol2 * sqrt_d3_square_tmp;
                d1 = sol0 * sqrt_d3_square_tmp;
                d2 = sol1 * sqrt_d3_square_tmp;
            }

            gauss_newton_refineL_wu<T, refinement_iterations>(d1, d2, sqrt_d3_square,a12,a13,a23,m12,m13,m23);
            // Vector3<T> y1_d1 = y1*d1;
            yd1=y1*d1-y2*d2;
            yd2=y1*d1-y3*sqrt_d3_square;
            yd1xd2=yd1.cross(yd2);
            Matrix<T,3,3> Y(yd1(0),yd2(0),yd1xd2(0),
                            yd1(1),yd2(1),yd1xd2(1),
                            yd1(2),yd2(2),yd1xd2(2));
            Rs[valid]=Y*X;
            Ts[valid]=(y1*d1 - Rs[valid]*x1);
            valid++;
        }
#else
        if(sol0*sol1<0)
        {
            continue;
        }
        T sol2_inv = 1/(H20*tmp_sol0+H21*tmp_sol1+H22);
        sol0 *= sol2_inv;
        sol1 *= sol2_inv;
        if(sol0>0 && sol1>0)
        {
            T denom = sol0 * sol0 + sol1 * sol1-2 * m12 * sol0 * sol1;
            if(denom<0)
                continue;
            T d3_square = a12 / denom;
            T sqrt_d3_square = std::sqrt(d3_square);
            T d1 = sol0 * sqrt_d3_square;
            T d2 = sol1 * sqrt_d3_square;
            gauss_newton_refineL_wu<T, refinement_iterations>(d1, d2, sqrt_d3_square,a12,a13,a23,m12,m13,m23);
            // Vector3<T> y1_d1 = y1*d1;
            yd1=y1*d1-y2*d2;
            yd2=y1*d1-y3*sqrt_d3_square;
            yd1xd2=yd1.cross(yd2);
            Matrix<T,3,3> Y(yd1(0),yd2(0),yd1xd2(0),
                            yd1(1),yd2(1),yd1xd2(1),
                            yd1(2),yd2(2),yd1xd2(2));
            Rs[valid]=Y*X;
            Ts[valid]=(y1*d1 - Rs[valid]*x1);
            valid++;
        }
#endif
    }

#if Time_tag
    TOC(lt2);
#endif 



    return valid;
}

}
