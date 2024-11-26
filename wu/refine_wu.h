#pragma once
#include <iostream>
#include <utils/cvl/matrix.h>
namespace cvl{

template<class T, int iterations>
/**
 * @brief refineL
 * @param L
 * @param a12
 * @param a13
 * @param a23
 * @param b12
 * @param b13
 * @param b23
 *
 * Gauss-Newton Solver
 * For unknown reasons it always works for the correct solution, but not always for the other solutions!
 *
 */
void gauss_newton_refineL_wu(T & l1, T & l2, T & l3, 
                          T a12, T a13, T a23,
                          T b12, T b13, T b23 ){

    // const expr makes it easier for the compiler to unroll
    for(int i=0;i<iterations;++i){
        T l1xl1= l1 * l1;
        T l2xl2= l2 * l2;
        T l3xl3= l3 * l3;
        T r1 = (l1xl1 - 2.0 * l1 * l2 * b12 + l2xl2 - a12);
        T r2 = (l1xl1 - 2.0 * l1 * l3 * b13 + l3xl3 - a13);
        T r3 = (l2xl2 - 2.0 * l2 * l3 * b23 + l3xl3 - a23);
        if (std::abs(r1) + std::abs(r2) + std::abs(r3) < 1e-10)
            return;
        T x11 = l1 - l2 * b12;
        T x12 = l2 - l1 * b12;
        T x21 = l1 - l3 * b13;
        T x23 = l3 - l1 * b13;
        T x32 = l2 - l3 * b23;
        T x33 = l3 - l2 * b23;
        T detJ = 0.5 / (x11 * x23 * x32 + x12 * x21 * x33); // half minus inverse determinant
        // This uses the closed form of the inverse for the jacobean.
        // Due to the zero elements this actually becomes quite nice.
        l1 += (-x23 * x32 * r1 - x12 * x33 * r2 + x12 * x23 * r3) * detJ;
        l2 += (-x21 * x33 * r1 + x11 * x33 * r2 - x11 * x23 * r3) * detJ;
        l3 += (x21 * x32 * r1 - x11 * x32 * r2 - x12 * x21 * r3) * detJ;

    }
    // cout<<i<<endl;


}


}
