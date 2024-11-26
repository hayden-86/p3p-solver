#pragma once

#include <ding/solve_cubic_ding.h>
#include <ding/refine_lambda_ding.h>
#include <lambdatwist/p3p_timers.h>

namespace cvl
{

    template <class T, int refinement_iterations = 5>
    int p3p_ding_null(Vector3<T> x0,
                Vector3<T> x1,
                Vector3<T> x2,
                Vector3<T> X0,
                Vector3<T> X1,
                Vector3<T> X2,
                Vector<cvl::Matrix<T, 3, 3>, 4> &Rs,
                Vector<Vector3<T>, 4> &Ts)
    {

        // normalize the length of ys, we could expect it, but lets not...
        TIC(lt1);
        x0.normalize();
        x1.normalize();
        x2.normalize();

        // implicit creation of Vector3<T> can be removed
        Vector3<T> X01 = X0 - X1;
        Vector3<T> X02 = X0 - X2;
        Vector3<T> X12 = X1 - X2;
        
        T a01 = X01.squaredNorm();
        T a02 = X02.squaredNorm();
        T a12 = X12.squaredNorm();

        if (a01 > a02) {
            if (a01 > a12) {
                std::swap(x0, x2);
                std::swap(X0, X2);
                std::swap(a01, a12);
                X01 = -X12;
                X02 = -X02;
            }
        } else if (a02 > a12) {
            std::swap(x0, x1);
            std::swap(X0, X1);
            std::swap(a02, a12);
            X01 = -X01;
            X02 = X12;
        }

        T a12d = 1.0 / a12;
        T a = a01 * a12d;
        T b = a02 * a12d;

        T m01 = x0.dot(x1);
        T m02 = x0.dot(x2);
        T m12 = x1.dot(x2);

        // Ugly parameters to simplify the calculation
        T m12sq = -m12 * m12 + 1.0;
        T m02sq = -1.0 + m02 * m02;
        T m01sq = -1.0 + m01 * m01;
        T ab = a * b;
        T bsq = b * b;
        T asq = a * a;
        T m013 = -2.0 + 2.0 * m01 * m02 * m12;
        T bsqm12sq = bsq * m12sq;
        T asqm12sq = asq * m12sq;
        T abm12sq = 2.0 * ab * m12sq;

        T k3_inv = 1.0 / (bsqm12sq + b * m02sq);
        T k2 = k3_inv * ((-1.0 + a) * m02sq + abm12sq + bsqm12sq + b * m013);
        T k1 = k3_inv * (asqm12sq + abm12sq + a * m013 + (-1.0 + b) * m01sq);
        T k0 = k3_inv * (asqm12sq + a * m01sq);

        
        T s;
        bool G = cubick_ding(k2, k1, k0, s);
        
        T A00 = -a + s * (1 - b);
        T A01 = -m02 * s;
        T A02 = a * m12 + b * m12 * s;
        T A11 = s + 1;
        T A12 = -m01;
        T A22 = -a - b * s + 1;

        Matrix<T, 2, 2> Aiv(A01, A02,
                          A11, A12);
        Vector2d Ab(-A00,-A01);
        Vector2d yz = Aiv.inverse()*Ab;
        Vector3<T> Bp;
        Bp(0) = 1; Bp(1) = yz(0); Bp(2) = yz(1);
        Bp.normalize();
        Bp = Bp*std::sqrt(A01*A01+A02*A02+A12*A12-A00*A11-A00*A22-A22*A11);

        Matrix<T, 3, 3> A(A00, A01, A02,
                          A01, A11, A12,
                          A02, A12, A22); 

        Matrix<T, 3, 3> Mp;
        Mp = Bp.crossMatrix();

        Matrix<T, 3, 3> Cc;
        Cc = A + Mp;

        Matrix<T, 3, 2> pq;

        pq(0, 0) = Cc(0, 0);
        pq(1, 0) = Cc(1, 0);
        pq(2, 0) = Cc(2, 0);
        pq(0, 1) = Cc(0, 0);
        pq(1, 1) = Cc(0, 1);
        pq(2, 1) = Cc(0, 2);

        T d0, d1, d2;
        Vector3<T> x02x02(X01.cross(X02));
        Matrix<T, 3, 3> XX(X01(0), X02(0), x02x02(0),
                           X01(1), X02(1), x02x02(1),
                           X01(2), X02(2), x02x02(2));
        XX = XX.inverse();

        Vector3<T> v1, v2;
        Vector3<T> v1v2;
        // Matrix<T, 3, 3> YY;


        int n_sols = 0;

        for (int i = 0; i < 2; ++i) {
            // [p0 p1 p2] * [1; x; y] = 0, or [p0 p1 p2] * [d2; d0; d1] = 0
            T p0 = pq(0, i);
            T p1 = pq(1, i);
            T p2 = pq(2, i);
            // here we run into trouble if p0 is zero,
            // so depending on which is larger, we solve for either d0 or d1
            // The case p0 = p1 = 0 is degenerate and can be ignored
            bool switch_12 = std::abs(p0) <= std::abs(p1);

            if (switch_12) {
                // eliminate d0
                T w0 = -p0 / p1;
                T w1 = -p2 / p1;
                T ca = 1.0 / (w1 * w1 - b);
                T cb = 2.0 * (b * m12 - m02 * w1 + w0 * w1) * ca;
                T cc = (w0 * w0 - 2 * m02 * w0 - b + 1.0) * ca;
                Vector2<T> taus; 
                if (!root2real(cb, cc, taus[0], taus[1]))
                    continue;
                for (T tau : taus) {
                    if (tau <= 0)
                        continue;
                    // positive only
                    d2 = std::sqrt(a12 / (tau * (tau - 2.0 * m12) + 1.0));
                    d1 = tau * d2;
                    d0 = (w0 * d2 + w1 * d1);
                    if (d0 < 0)
                        continue;

                    gauss_newton_refineL_ding<T, refinement_iterations>(d0, d1, d2, a01, a02, a12, m01, m02, m12);
                    v1 = d0 * x0 - d1 * x1;
                    v2 = d0 * x0 - d2 * x2;
                    v1v2 = v1.cross(v2);
                    Matrix<T, 3, 3> Y(v1(0), v2(0), v1v2(0),
                              v1(1), v2(1), v1v2(1),
                              v1(2), v2(2), v1v2(2));
                              
                    Rs[n_sols] = Y * XX;
                    Ts[n_sols] = d0 * x0 - Rs[n_sols] * X0;
                    ++n_sols;
                }
            } else {
                T w0 = -p1 / p0;
                T w1 = -p2 / p0;
                T ca = 1.0 / (-a * w1 * w1 + 2 * a * m12 * w1 - a + 1);
                T cb = 2 * (a * m12 * w0 - m01 - a * w0 * w1) * ca;
                T cc = (1 - a * w0 * w0) * ca;

                Vector2<T> taus;
                if (!root2real(cb, cc, taus[0], taus[1]))
                    continue;
                for (T tau : taus) {
                    if (tau <= 0)
                        continue;
                    d0 = std::sqrt(a01 / (tau * (tau - 2.0 * m01) + 1.0));
                    d1 = tau * d0;
                    d2 = w0 * d0 + w1 * d1;

                    if (d2 < 0)
                        continue;

                    gauss_newton_refineL_ding<T, refinement_iterations>(d0, d1, d2, a01, a02, a12, m01, m02, m12);
                    v1 = d0 * x0 - d1 * x1;
                    v2 = d0 * x0 - d2 * x2;

                    v1v2 = v1.cross(v2);
                    Matrix<T, 3, 3> Y(v1(0), v2(0), v1v2(0),
                              v1(1), v2(1), v1v2(1),
                              v1(2), v2(2), v1v2(2));

                    Rs[n_sols] = Y * XX;
                    Ts[n_sols] = d0 * x0 - Rs[n_sols] * X0;
                    ++n_sols;
                }
            }

            if (n_sols > 0 && G)
                break;
        }


        return n_sols;
    }
}
