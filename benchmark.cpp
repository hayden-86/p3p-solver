#include "benchmark.h"
#include "lambdatwist/lambdatwist.p3p.h"
#include "ding/ding_direct.h"
#include "ding/ding_null.h"
#include "ding/ding_adjoint.h"
#include "kneip/kneip.h"
#include "kneip/nakano.h"
#include "utils/mlibtime.h"
#include <ke/ke.h>
#include <p3p_generator.h>
#include <ke/ke_utils.h>
#include <data.h>
#include <kneip/kneip_utils.h>
#include <kneip/nakano_utils.h>

#include "utils/string_helpers.h"

#include <fstream>
#include "new/new.p3p.h"
#include "new/new.p3p_direct.h"
#include "new/new.p3p_method1.h"
#include "new/new.p3p_method2.h"
#include "wu/wu1.h"
/*********
 *
 *
 *
 * Standardize input output format,
 * The speed comparison is primarily between two methods that output rotation matrixes. Therefore use that as basis.
 * The ke uses a silly format to maximize their speed...
 * convert to something resonable on output since it would always be followed by that anyways.
 *
 * Internally the methods should test the feasibility of the returned solutions, given the three points.
 * Ie dont return geometrically invalid solutions, nans, non rotations, transforms that move points behind cameras, or large reprojection errors
 * If they fail to do this it will always be done directly after.
 *
 * It would be very natural to use inheritance to create a nice modular structure for the benchmark. Doing so would hide the performance however,
 *  as vtable lookups take quite some time compared to the 300ns or so that the routines take to run... yes I tried it...
 *
 *
 *   Terms:
 * a valid solution is one considered valid by the solver, should always atleast test for geometric feasibility and nans
 * a correct solution is one which is both valid and satisfies the rotation matrix and reprojection criteria.
 *
 *
 * **/
#include <sstream>
#include <solver.h>

using std::cout;
using std::endl;
using namespace cvl;

#define TOTAL_ACCURACY 0


namespace mlib
{
    namespace klas
    {

        template <class T, uint Rows, uint Cols>
        std::string displayMatrix(cvl::Matrix<T, Rows, Cols> M, bool displayrowcol)
        {

            std::vector<std::string> headers, rownames;
            if (displayrowcol)
            {
                for (uint i = 0; i < Cols; ++i)
                    headers.push_back(mlib::toStr(i));
                for (uint i = 0; i < Rows; ++i)
                    rownames.push_back(mlib::toStr(i) + ": ");
            }
            std::vector<std::vector<T>> rows;
            for (uint r = 0; r < Rows; ++r)
            {
                std::vector<T> row;
                for (uint c = 0; c < Cols; ++c)
                    row.push_back(M(r, c));
                rows.push_back(row);
            }

            return mlib::displayTable(headers, rows, rownames);
        }

        template <class T>
        class KeSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {
                return kes::ke_p3p_fair(data, Rs, Ts);
            }
            std::string get_name() { return "Ke"; }
        };
        template <class T>
        class KneipSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {
                return kneip::kneip_p3p_fair(data, Rs, Ts);
            }
            std::string get_name() { return "Kneip"; }
        };

        template <class T>
        class nakanoSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {
                return nakano::nakano_p3p_fair(data, Rs, Ts);
            }
            std::string get_name() { return "nakano"; }
        };

        template <class T>
        class LambdaSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_lambdatwist(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "Lambda"; }
        };

        template <class T>
        class dingSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_ding(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "ding"; }
        };

        template <class T>
        class dingSolverd
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_ding_direct(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "direct"; }
        };

        template <class T>
        class dingSolver1
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_ding_null(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "null"; }
        };

        template <class T>
        class dingSolver2
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_ding_adjoint(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "ding_new"; }
        };

        template <class T>
        class NewSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "New"; }
        };

        template <class T>
        class NewSolverd
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new_direct(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "ding_old"; }
        };

        template <class T>
        class NewSolver1
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new_method1(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "method1_new"; }
        };

        template <class T>
        class NewSolver2
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new_method2(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "method2_new"; }
        };

        // verification of answers can be made slowly, so...
        template <class T>
        class WuSolver1
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {
                // wuhd: x0-(X,Y,Z) 3D points;   xr-(x,y,z) camera coordinates;   ys-(xs,ys) 2Dimage points
                Vector3<Vector3<T>> yss = data.xr; // camera (x,y,z), 3 groups
                Vector3<Vector3<T>> xss = data.x0; // world (X,Y,Z), 3 groups
                // wuhd: return - valid(the biggest number is 4), Rs(4), Ts(4)
                return p3p_wu1(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "ours"; }
        };

        template <class T, class Solver>
        P3PResult compute_accuracy(Solver S,
                                   std::vector<Data<T>> datas,
                                   T error_limit = 1e-6)
        {
            P3PResult res(S.get_name(), datas.size());

            std::ofstream outFile;

            for (Data<T> &data : datas)
            {
                cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                cvl::Vector<Vector<T, 3>, 4> Ts;

                int valid = S.solve(data, Rs, Ts);
                res.valid += valid; // valid according to solver.

                int duplicates = 0;
                int sols = data.good_solutions(Rs, Ts, valid, duplicates); // correct, with duplicates included

                if (sols == 0)
                    res.no_solution++;

                if (valid > sols)
                    res.incorrect_valid += (valid - sols);

                res.solutions += (sols - duplicates); // correct unique
                res.duplicates += duplicates;

                T error = data.min_error(Rs, Ts, valid);

                if (error < error_limit)
                {
                    res.errors.push_back(error);
                    res.ground_truth_in_set++;
                }

                // res.errors.push_back(error);
                // if (error < error_limit)
                //     res.ground_truth_in_set++;

            }
            outFile.close();
            return res;
        }

        template <class T>
        void test(const std::vector<Data<T>> &datas)
        {

            cout << "Beginning Test: " << endl;
            T error_limit = 1e-6;

            P3PResult wu1 = compute_accuracy<T, WuSolver1<T>>(WuSolver1<T>(), datas, error_limit);
            P3PResult newp3p_direct_new = compute_accuracy<T, NewSolverd<T>>(NewSolverd<T>(), datas, error_limit);
            // P3PResult newp3p_method1 = compute_accuracy<T, NewSolver1<T>>(NewSolver1<T>(), datas, error_limit);
            // P3PResult newp3p_method2 = compute_accuracy<T, NewSolver2<T>>(NewSolver2<T>(), datas, error_limit);
            // P3PResult dingp3p_direct = compute_accuracy<T, dingSolverd<T>>(dingSolverd<T>(), datas, error_limit);
            // P3PResult dingp3p_null = compute_accuracy<T, dingSolver1<T>>(dingSolver1<T>(), datas, error_limit);
            P3PResult dingp3p_adjoint = compute_accuracy<T, dingSolver2<T>>(dingSolver2<T>(), datas, error_limit);
            P3PResult lambda = compute_accuracy<T,LambdaSolver<T>>(LambdaSolver<T>(),datas,error_limit);
            P3PResult ke = compute_accuracy<T,KeSolver<T>>(KeSolver<T>(),datas,error_limit);
            P3PResult kneip= compute_accuracy<T,KneipSolver<T>>(KneipSolver<T>(),datas,error_limit);
            P3PResult nakano= compute_accuracy<T,nakanoSolver<T>>(nakanoSolver<T>(),datas,error_limit);

#if TOTAL_ACCURACY // total accuracy
            static double data_size=0;
            static int i=0;
            data_size += datas.size();
            i++;
            static long int ground_truth_in_set_wu1=0;
            static long int no_solution_wu1=0;
            static long int valid_wu1=0;
            static long int duplicates_wu1=0;
            static long int solutions_wu1=0;
            static long int incorrect_valid_wu1=0;
            ground_truth_in_set_wu1 += wu1.ground_truth_in_set; 
            no_solution_wu1 += wu1.no_solution; 
            valid_wu1 += wu1.valid; 
            duplicates_wu1 += wu1.duplicates; 
            solutions_wu1 += wu1.solutions; 
            incorrect_valid_wu1 += wu1.incorrect_valid; 

            static long int ground_truth_in_set_direct_new=0;
            static long int no_solution_direct_new=0;
            static long int valid_direct_new=0;
            static long int duplicates_direct_new=0;
            static long int solutions_direct_new=0;
            static long int incorrect_valid_direct_new=0;
            ground_truth_in_set_direct_new += newp3p_direct_new.ground_truth_in_set; 
            no_solution_direct_new += newp3p_direct_new.no_solution; 
            valid_direct_new += newp3p_direct_new.valid; 
            duplicates_direct_new += newp3p_direct_new.duplicates; 
            solutions_direct_new += newp3p_direct_new.solutions; 
            incorrect_valid_direct_new += newp3p_direct_new.incorrect_valid; 

            static long int ground_truth_in_set_adjoint=0;
            static long int no_solution_adjoint=0;
            static long int valid_adjoint=0;
            static long int duplicates_adjoint=0;
            static long int solutions_adjoint=0;
            static long int incorrect_valid_adjoint=0;
            ground_truth_in_set_adjoint += dingp3p_adjoint.ground_truth_in_set; 
            no_solution_adjoint += dingp3p_adjoint.no_solution; 
            valid_adjoint += dingp3p_adjoint.valid; 
            duplicates_adjoint += dingp3p_adjoint.duplicates; 
            solutions_adjoint += dingp3p_adjoint.solutions; 
            incorrect_valid_adjoint += dingp3p_adjoint.incorrect_valid; 

            static long int ground_truth_in_set_Lambda=0;
            static long int no_solution_Lambda=0;
            static long int valid_Lambda=0;
            static long int duplicates_Lambda=0;
            static long int solutions_Lambda=0;
            static long int incorrect_valid_Lambda=0;
            ground_truth_in_set_Lambda += lambda.ground_truth_in_set; 
            no_solution_Lambda += lambda.no_solution; 
            valid_Lambda += lambda.valid; 
            duplicates_Lambda += lambda.duplicates; 
            solutions_Lambda += lambda.solutions; 
            incorrect_valid_Lambda += lambda.incorrect_valid; 

            static long int ground_truth_in_set_Ke=0;
            static long int no_solution_Ke=0;
            static long int valid_Ke=0;
            static long int duplicates_Ke=0;
            static long int solutions_Ke=0;
            static long int incorrect_valid_Ke=0;
            ground_truth_in_set_Ke += ke.ground_truth_in_set; 
            no_solution_Ke += ke.no_solution; 
            valid_Ke += ke.valid; 
            duplicates_Ke += ke.duplicates; 
            solutions_Ke += ke.solutions; 
            incorrect_valid_Ke += ke.incorrect_valid; 

            static long int ground_truth_in_set_Kneip=0;
            static long int no_solution_Kneip=0;
            static long int valid_Kneip=0;
            static long int duplicates_Kneip=0;
            static long int solutions_Kneip=0;
            static long int incorrect_valid_Kneip=0;
            ground_truth_in_set_Kneip += kneip.ground_truth_in_set; 
            no_solution_Kneip += kneip.no_solution; 
            valid_Kneip += kneip.valid; 
            duplicates_Kneip += kneip.duplicates; 
            solutions_Kneip += kneip.solutions; 
            incorrect_valid_Kneip += kneip.incorrect_valid; 

            static long int ground_truth_in_set_nakano=0;
            static long int no_solution_nakano=0;
            static long int valid_nakano=0;
            static long int duplicates_nakano=0;
            static long int solutions_nakano=0;
            static long int incorrect_valid_nakano=0;
            ground_truth_in_set_nakano += nakano.ground_truth_in_set; 
            no_solution_nakano += nakano.no_solution; 
            valid_nakano += nakano.valid; 
            duplicates_nakano += nakano.duplicates; 
            solutions_nakano += nakano.solutions; 
            incorrect_valid_nakano += nakano.incorrect_valid; 
            if (i%10==0)
            {
                const int colWidth = 12;

                std::cout << "\n\n--------------------------------------------------------------------------------------------------------------------------------------------\n";
                std::cout << "------------------------------------------------------------------Wu Table------------------------------------------------------------------\n";
                std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------\n";
                std::cout << "The "<< i << "th iteration"<<std::endl;
                std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------\n";
                std::cout << std::left<< std::setw(48) << "Table"
                        << std::right
                        << std::setw(colWidth) << "wu1"
                        << std::setw(colWidth) << "direct_new"
                        << std::setw(colWidth) << "adjoint"
                        << std::setw(colWidth) << "Lambda"
                        << std::setw(colWidth) << "Ke"
                        << std::setw(colWidth) << "Kneip"
                        << std::setw(colWidth) << "nakano" 
                        << std::endl;
                std::cout << std::left<< std::setw(48) << "ground truth found"
                        << std::right
                        << std::setw(colWidth) << ground_truth_in_set_wu1 
                        << std::setw(colWidth) << ground_truth_in_set_direct_new
                        << std::setw(colWidth) << ground_truth_in_set_adjoint 
                        << std::setw(colWidth) << ground_truth_in_set_Lambda
                        << std::setw(colWidth) << ground_truth_in_set_Ke 
                        << std::setw(colWidth) << ground_truth_in_set_Kneip
                        << std::setw(colWidth) << ground_truth_in_set_nakano 
                        << std::endl;

                std::cout << std::left<< std::setw(48) << "any solution found"
                        << std::right
                        << std::setw(colWidth) << (data_size - no_solution_wu1) 
                        << std::setw(colWidth) << (data_size - no_solution_direct_new)
                        << std::setw(colWidth) << (data_size - no_solution_adjoint)
                        << std::setw(colWidth) << (data_size - no_solution_Lambda)
                        << std::setw(colWidth) << (data_size - no_solution_Ke) 
                        << std::setw(colWidth) << (data_size - no_solution_Kneip)
                        << std::setw(colWidth) << (data_size - no_solution_nakano)
                        << std::endl;

                std::cout << std::left<< std::setw(48) << "ground truth found ratio"
                        << std::right
                        << std::setw(colWidth) << ground_truth_in_set_wu1 / data_size
                        << std::setw(colWidth) << ground_truth_in_set_direct_new / data_size
                        << std::setw(colWidth) << ground_truth_in_set_adjoint / data_size
                        << std::setw(colWidth) << ground_truth_in_set_Lambda / data_size
                        << std::setw(colWidth) << ground_truth_in_set_Ke / data_size
                        << std::setw(colWidth) << ground_truth_in_set_Kneip / data_size
                        << std::setw(colWidth) << ground_truth_in_set_nakano / data_size
                        << std::endl;

                std::cout << std::left<< std::setw(48) << "no solution found"
                        << std::right
                        << std::setw(colWidth) << no_solution_wu1 
                        << std::setw(colWidth) << no_solution_direct_new
                        << std::setw(colWidth) << no_solution_adjoint 
                        << std::setw(colWidth) << no_solution_Lambda
                        << std::setw(colWidth) << no_solution_Ke 
                        << std::setw(colWidth) << no_solution_Kneip
                        << std::setw(colWidth) << no_solution_nakano 
                        << std::endl;

                std::cout << std::left<< std::setw(48) << "valid according to solver"
                        << std::right
                        << std::setw(colWidth) << valid_wu1 
                        << std::setw(colWidth) << valid_direct_new
                        << std::setw(colWidth) << valid_adjoint 
                        << std::setw(colWidth) << valid_Lambda
                        << std::setw(colWidth) << valid_Ke 
                        << std::setw(colWidth) << valid_Kneip
                        << std::setw(colWidth) << valid_nakano 
                        << std::endl;

                std::cout << std::left<< std::setw(48) << "duplicates"
                        << std::right
                        << std::setw(colWidth) << duplicates_wu1 
                        << std::setw(colWidth) << duplicates_direct_new
                        << std::setw(colWidth) << duplicates_adjoint 
                        << std::setw(colWidth) << duplicates_Lambda
                        << std::setw(colWidth) << duplicates_Ke 
                        << std::setw(colWidth) << duplicates_Kneip
                        << std::setw(colWidth) << duplicates_nakano 
                        << std::endl;

                std::cout << std::left<< std::setw(48) << "unique correct solutions"
                        << std::right
                        << std::setw(colWidth) << solutions_wu1 
                        << std::setw(colWidth) << solutions_direct_new
                        << std::setw(colWidth) << solutions_adjoint 
                        << std::setw(colWidth) << solutions_Lambda
                        << std::setw(colWidth) << solutions_Ke 
                        << std::setw(colWidth) << solutions_Kneip
                        << std::setw(colWidth) << solutions_nakano 
                        << std::endl;

                std::cout << std::left<< std::setw(48) << "incorrect solutions output by the solver"
                        << std::right
                        << std::setw(colWidth) << incorrect_valid_wu1 
                        << std::setw(colWidth) << incorrect_valid_direct_new
                        << std::setw(colWidth) << incorrect_valid_adjoint 
                        << std::setw(colWidth) << incorrect_valid_Lambda
                        << std::setw(colWidth) << incorrect_valid_Ke 
                        << std::setw(colWidth) << incorrect_valid_Kneip
                        << std::setw(colWidth) << incorrect_valid_nakano 
                        << std::endl;
                std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------\n";
                std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------\n\n\n";
            }
#endif

            // make the result table
            // std::vector<P3PResult> res = {wu1, newp3p_direct_new, newp3p_method1, newp3p_method2, dingp3p_direct, dingp3p_null, dingp3p_adjoint,lambda,ke,kneip,nakano};
            std::vector<P3PResult> res = {wu1, newp3p_direct_new, dingp3p_adjoint,lambda,ke,kneip,nakano};
            std::vector<std::string> headers;
            std::vector<std::string> columns;
            std::vector<std::vector<std::string>> rows;

            for (auto &r : res)
                headers.push_back(r.name);

            {
                // number of times ground truth was found
                std::vector<int> gt;
                for (auto &r : res)
                    gt.push_back(r.ground_truth_in_set);
                rows.push_back(toStrVec(gt));
                columns.push_back("ground truth found");
            }
            {
                // number of times a valid solution was found
                std::vector<int> gt;
                for (auto &r : res)
                    gt.push_back(datas.size() - r.no_solution);
                rows.push_back(toStrVec(gt));
                columns.push_back("any solution found");
            }
            {
                // ratio for ground truth found
                std::vector<double> gtratio;
                for (auto &r : res)
                    gtratio.push_back((double)r.ground_truth_in_set / datas.size());
                rows.push_back(toStrVec(gtratio));
                columns.push_back("ground truth found ratio");
            }

            {
                // number of no solution at all found
                std::vector<int> nosol;
                for (auto &r : res)
                    nosol.push_back(r.no_solution);
                rows.push_back(toStrVec(nosol));
                columns.push_back("no solution found");
            }
            {
                // number of valid solutions according to solver
                std::vector<int> solutions;
                for (auto &r : res)
                    solutions.push_back(r.valid);
                rows.push_back(toStrVec(solutions));
                columns.push_back("valid according to solver");
            }
            {
                // duplicates
                std::vector<int> duplicates;
                for (auto &r : res)
                    duplicates.push_back(r.duplicates);
                rows.push_back(toStrVec(duplicates));
                columns.push_back("duplicates");
                //
            }
            {
                // unique correct
                std::vector<int> usols;
                for (auto &r : res)
                    usols.push_back(r.solutions);
                rows.push_back(toStrVec(usols));
                columns.push_back("unique correct solutions");
                //
            }
            {
                // incorrect valid solutions
                std::vector<int> solutions;
                for (auto &r : res)
                    solutions.push_back(r.incorrect_valid);
                rows.push_back(toStrVec(solutions));
                columns.push_back("incorrect solutuons output by the solver");
            }

            cout << displayTable(headers, rows, columns, "Table: ") << endl;
            {
                std::vector<std::vector<float128>> errors;
                for (P3PResult r : res)
                    errors.push_back(r.errors);
                std::string str = display(errors, headers);
                cout << str<< endl;              
#if 0 // save errors
                static int i = 1;  
                std::ofstream outFile("../data"+ std::to_string(i)+".txt");
                i++;
                outFile << std::fixed << std::setprecision(15);
                for (const auto& vec : errors) {
                    for (float num : vec) {
                        outFile << num << " ";
                    }
                    outFile << std::endl;
                }
                outFile.close();
#endif                
            }

        }

        template <class T>
        void versus(const std::vector<Data<T>> &datas,
                    int repeat_versus = 5,
                    bool inner_timers = false)
        {

            cout << "Beginning Versus" << endl;

            Timer kneip("Kneip: ", datas.size());
            Timer nakano("nakano: ", datas.size());
            Timer lt("Lambda: ", datas.size());
            Timer kes("ke ", datas.size());
            // Timer nsd("direct: ", datas.size());
            // Timer ns1("m1: ", datas.size());
            Timer ns2("m2: ", datas.size());
            Timer nsd_new("direct_new: ", datas.size());
            // Timer ns1_new("m1_new: ", datas.size());
            // Timer ns2_new("m2_new: ", datas.size());
            Timer wu1("wu1: ", datas.size());

            Timer tot_lambda("Lambda Total");
            Timer tot_kes("ke Total");
            Timer tot_kneip("Kneip Total");
            Timer tot_nakano("nakano Total");
            // Timer tot_dingd("ding_direct Total");
            // Timer tot_ding1("ding_null Total");
            Timer tot_ding2("ding_new Total");
            Timer tot_newd("ding_old Total");
            // Timer tot_new1("m1_new Total");
            // Timer tot_new2("m2_new Total");
            Timer tot_wu1("ours1 Total");

            // dont opt away
            int sols = 0;
            // ke

            bool dokneip = true;
            for (int i = 0; i < repeat_versus; ++i)
            {

                // kneip
                 if(dokneip){

                    tot_kneip.tic();
                    for(uint i=0;i<datas.size();++i){
                        //if(inner_timers){if(i==100){                kneip.clear();            }kneip.tic();}

                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        sols+= kneip::kneip_p3p_fair(datas[i],Rs, Ts);

                        //    if(inner_timers)kneip.toc();
                    }

                    tot_kneip.toc();
                }
                // nakano
                if(dokneip){

                    tot_nakano.tic();
                    for(uint i=0;i<datas.size();++i){
                        //if(inner_timers){if(i==100){                kneip.clear();            }kneip.tic();}

                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        sols+= nakano::nakano_p3p_fair(datas[i],Rs, Ts);

                        //    if(inner_timers)kneip.toc();
                    }

                    tot_nakano.toc();
                }
                if(dokneip){
                    tot_kes.tic();

                    for(uint i=0;i<datas.size();++i){
                        //    if(inner_timers)if(i==100){            kes.clear();        }        kes.tic();
                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        sols+=kes::ke_p3p_fair(datas[i],Rs,Ts);
                        //      if(inner_timers)kes.toc();

                    }
                    tot_kes.toc();
                }
                /*
                {
                    tot_dingd.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_ding_direct(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_dingd.toc();
                }
                {
                    tot_ding1.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_ding_null(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_ding1.toc();
                }
                */
                {
                    tot_ding2.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_ding_adjoint(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_ding2.toc();
                }
                {
                    tot_lambda.tic();
                    for(uint i=0;i<datas.size();++i){
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        Data<T> data=datas[i];
                        Vector3<Vector<T,3>> yss=data.xr;
                        Vector3<Vector<T,3>> xss=data.x0;
                        sols+=p3p_lambdatwist(yss[0],yss[1],yss[2],xss[0],xss[1],xss[2],Rs,Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_lambda.toc();
                }

                {
                    tot_newd.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_new_direct(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_newd.toc();
                }
                /*
                {
                    tot_new1.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_new_method1(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_new1.toc();
                }
                {
                    tot_new2.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_new_method2(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_new2.toc();
                }
                */
                {
                    tot_wu1.tic();
                    for(uint i=0;i<datas.size();++i){
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        Data<T> data=datas[i];
                        Vector3<Vector<T,3>> yss=data.xr;
                        Vector3<Vector<T,3>> xss=data.x0;
                        sols+=p3p_wu1(yss[0],yss[1],yss[2],xss[0],xss[1],xss[2],Rs,Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_wu1.toc();
                }
                // std::vector<mlib::Timer> tots = {tot_lambda, tot_kes, tot_kneip, tot_nakano,tot_dingd,tot_ding1,tot_ding2,tot_newd,tot_new1,tot_new2, tot_wu1};
                std::vector<mlib::Timer> tots = {tot_lambda, tot_kes, tot_kneip, tot_nakano,tot_ding2,tot_newd, tot_wu1};

                // printtimers();
                printtimers_wu();
                std::vector<mlib::Timer> ts = {lt, kes, kneip, nakano, ns2,nsd_new, wu1};
                // if ((i+1)%10==0)
                cout << tots << endl
                     << "wrote the ts" << endl;
            }
            if (false)
            {

                std::ofstream out("timing.csv");

                {
                    auto ts = kneip.getTimes();
                    for (Time time : ts)
                    {

                        out << (int)(time.ns) << ", ";
                    }
                    out << "\n";
                }
                {
                    auto ts = lt.getTimes();
                    for (Time time : ts)
                    {

                        out << (int)(time.ns) << ", ";
                    }
                    out << "\n";
                }
                {
                    auto ts = kes.getTimes();
                    for (Time time : ts)
                    {

                        out << (int)(time.ns) << ", ";
                    }
                    out << "\n";
                }
                out << std::endl;
                out.close();
                // cout<<"means: "<<mean(kn)<<" "<<mean(ls)<<endl;
            }

            // cout<<mlib::display(samples,false)<<endl;
            cout << "Versus Done !" << endl;
        }

        template <class T>
        void test_type()
        {
            // generate data
            cout << "generating";
            cout.flush();
            mlib::Timer timer("generator");

            Generator gennie;

            int N = std::pow(10, 7);
            // int N = std::pow(10, 5);
            // int N = std::pow(10, 2);

            std::vector<Data<T>> datas;
            datas.reserve(N);
            timer.tic();
            for (int i = 0; i < N; ++i)
            {

                datas.push_back(gennie.next<T>());

                // if((i%1000)==999)        cout<<"i: "<<i<<endl;
            }
            timer.toc();

            cout << timer << " done " << datas.size() << endl;

            test(datas);
            // printtimers();

            // versus(datas,3);
            versus(datas,100);

            // testReal();
        }

        void test_special()
        {
            // generate data
            cout << "generating";
            cout.flush();
            mlib::Timer timer("generator");

            Generator gennie;

            int N = std::pow(10, 4);

            std::vector<Data<double>> datas;
            datas.reserve(N);
            timer.tic();
            for (int i = 0; i < N; ++i)
            {

                datas.push_back(gennie.special_case0<double>());

                // if((i%1000)==999)        cout<<"i: "<<i<<endl;
            }
            timer.toc();

            cout << timer << " done " << datas.size() << endl;

            test(datas);
            // printtimers();

            versus(datas);

            // testReal();
        }

        template <class T>
        void profile_lambda()
        {

            Generator gennie;

            int sols = 0;
            cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
            cvl::Vector<Vector<T, 3>, 4> Ts;
            Data<T> data = gennie.next<T>();

            for (int i = 0; i < 100000; ++i)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                sols += p3p_lambdatwist(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            cout << sols << endl;
        }

        void testAll()
        {

            test_type<double>();
            // test_special();
            // profile_lambda<double>();
        }
    } // end namespace klas
} // end namespace mlib
