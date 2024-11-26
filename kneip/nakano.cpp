#include "kneip/nakano.h"
#include <iostream>


#include <Eigen/Dense>
#include <iostream>

// Function for root polishing using Gauss-Newton refinement
void rootpolishing(const Eigen::Matrix<double, 6, 1>& f,
                   const Eigen::Matrix<double, 6, 1>& g,
                   Eigen::VectorXd& x,
                   Eigen::VectorXd& y,
                   int iterations) {
    for (int i = 0; i < iterations; ++i) {
        // Compute x^2, xy, y^2 using element-wise operations
        Eigen::ArrayXd x2 = x.array().square();
        Eigen::ArrayXd xy = x.array() * y.array();
        Eigen::ArrayXd y2 = y.array().square();

        // Compute fv and gv using element-wise operations
        Eigen::ArrayXd fv = f(0) * x2 + f(1) * xy + f(3) * x.array() + f(4) * y.array() + f(5);
        Eigen::ArrayXd gv = g(0) * x2 - y2 + g(3) * x.array() + g(4) * y.array() + g(5);

        // Find indices where both fv and gv are close to zero
        Eigen::Array<bool, Eigen::Dynamic, 1> ii = (fv.abs() < 1e-14) && (gv.abs() < 1e-14);
        if (ii.all()) {
            return;
        }
        // Compute partial derivatives dfdx, dfdy, dgdx, dgdy using element-wise operations
        Eigen::ArrayXd dfdx = 2 * f(0) * x.array() + f(1) * y.array() + f(3);
        Eigen::ArrayXd dfdy = f(1) * x.array() + f(4);
        Eigen::ArrayXd dgdx = 2 * g(0) * x.array() + g(3);
        Eigen::ArrayXd dgdy = -2 * y.array() + g(4);

        // Compute the inverse of the determinant of the Jacobian matrix
        Eigen::ArrayXd inv_detJ = 1.0 / (dfdx * dgdy - dfdy * dgdx);

        // Compute dx and dy
        Eigen::ArrayXd dx = (dgdy * fv - dfdy * gv) * inv_detJ;
        Eigen::ArrayXd dy = (-dgdx * fv + dfdx * gv) * inv_detJ;

        // // Set dx and dy to zero where the condition is met
        // dx = dx * (!ii).cast<double>();
        // dy = dy * (!ii).cast<double>();

        // Update x and y
        x = x.array() - dx;
        y = y.array() - dy;
    }
}


int P3Pnokano::computePoses2(const Eigen::Matrix3d &feature_vectors, const Eigen::Matrix3d &world_points,
                      Eigen::Matrix<Eigen::Matrix<double, 3, 4>, 4, 1> &solutions)
{
  // Extraction of world points
  Eigen::Vector3d X1 = world_points.col(0);
  Eigen::Vector3d X2 = world_points.col(1);
  Eigen::Vector3d X3 = world_points.col(2);

  Eigen::Vector3d X32 = X3 - X2;
  Eigen::Vector3d X31 = X3 - X1;
  Eigen::Vector3d X21 = X2 - X1;

  Eigen::Matrix3d m = feature_vectors;
  Eigen::Matrix3d X = world_points;

  if (X32.squaredNorm() > X31.squaredNorm() && X32.squaredNorm() > X21.squaredNorm())
  {
    m.col(0) = feature_vectors.col(1);
    m.col(1) = feature_vectors.col(2);
    m.col(2) = feature_vectors.col(0);
    X.col(0) = world_points.col(1);
    X.col(1) = world_points.col(2);
    X.col(2) = world_points.col(0);
  }
  else if (X31.squaredNorm() > X32.squaredNorm() && X31.squaredNorm() > X21.squaredNorm())
  {
    m.col(1) = feature_vectors.col(2);
    m.col(2) = feature_vectors.col(1);
    X.col(1) = world_points.col(2);
    X.col(2) = world_points.col(1);
  }

  X21 = X.col(1) - X.col(0);
  X31 = X.col(2) - X.col(0);

  Eigen::Vector3d nx = X21;
  nx = nx / nx.norm();
  Eigen::Vector3d nz = nx.cross(X31);
  nz = nz / nz.norm();
  Eigen::Vector3d ny = nz.cross(nx);
  Eigen::Matrix3d N;
  N.col(0) = nx;
  N.col(1) = ny;
  N.col(2) = nz;

  Eigen::Vector3d mc0 = m.col(0);
  Eigen::Vector3d mc1 = m.col(1);
  Eigen::Vector3d mc2 = m.col(2);
  

  double a = nx.transpose() * (X21);
  double b = nx.transpose() * (X31);
  double c = ny.transpose() * (X31);

  double M12 = mc0.transpose() * mc1;
  double M13 = mc0.transpose() * mc2;
  double M23 = mc1.transpose() * mc2;
  double p = b / a;
  double q = (b * b + c * c) / (a * a);

  Eigen::Matrix<double, 6, 1> f;
  f << p, -M23, 0, -M12 * (2.0 * p - 1.0), M13, p - 1.0;

  Eigen::Matrix<double, 6, 1> g;
  g << q, 0, -1.0, -2.0 * M12 * q, 2.0 * M13, q - 1.0;

  Eigen::Matrix<double, 5, 1> h;
  h(0) = -f(0) * f(0) + g(0) * f(1) * f(1);
  h(1) = f(1) * f(1) * g(3) - 2.0 * f(0) * f(3) - 2.0 * f(0) * f(1) * f(4) + 2.0 * f(1) * f(4) * g(0);
  h(2) = f(4) * f(4) * g(0) - 2.0 * f(0) * f(4) * f(4) - 2.0 * f(0) * f(5) + f(1) * f(1) * g(5) - f(3) * f(3) - 2.0 * f(1) * f(3) * f(4) + 2.0 * f(1) * f(4) * g(3);
  h(3) = f(4) * f(4) * g(3) - 2.0 * f(3) * f(4) * f(4) - 2.0 * f(3) * f(5) - 2.0 * f(1) * f(4) * f(5) + 2.0 * f(1) * f(4) * g(5);
  h(4) = -2 * f(4) * f(4) * f(5) + g(5) * f(4) * f(4) - f(5) * f(5);

  Eigen::Matrix<double, 4, 1> xxx;

  int nsols = P3Pnokano::solveQuartic2(h, xxx);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(nsols, 1);
  for (int i = 0; i < nsols; ++i)
  {
    x(i) = xxx(i);
  }

  Eigen::VectorXd y = Eigen::VectorXd::Zero(nsols, 1);
  for (int i = 0; i < nsols; ++i)
  {
    y(i) = -((f(0) * x(i) + f(3)) * x(i) + f(5)) / (f(4) + f(1) * x(i));
  }
  
  // f << 0.7397, -0.9603, 0, -0.4749, 0.9861, -0.2603;
  // g << 0.8074, 0, -1.0000, -1.5997, 1.9723, -0.1926;
  // x = Eigen::VectorXd::Zero(2, 1);
  // y = Eigen::VectorXd::Zero(2, 1);
  // x << 1.9142, 0.6667;
  // y << 1.8084, 0.7174;
  rootpolishing(f,g,x,y,5);



  Eigen::Matrix3d A;
  A.col(0) = -mc0;
  A.col(1) = mc1;
  A.col(2) <<0,0,0;

  Eigen::Matrix3d B;
  B.col(0) = -mc0;
  B.col(1) << 0,0,0;
  B.col(2) = mc2;

  Eigen::Matrix3d C = B - p * A;

  for (int i = 0; i < nsols; ++i)
  {
    Eigen::Vector3d lambda;
    lambda << 1, x(i), y(i);
    double s = (A * lambda).norm() / a;
    Eigen::Vector3d d = lambda / s;
    Eigen::Vector3d r1 = (A * d) / a;
    Eigen::Vector3d r2 = (C * d) / c;
    Eigen::Vector3d r3 = r1.cross(r2);
    Eigen::Matrix3d Rc;
    Rc.col(0) = r1;
    Rc.col(1) = r2;
    Rc.col(2) = r3;
    Eigen::Vector3d tc = d(0) * mc0;

    Eigen::Matrix3d R = Rc * N.transpose();
    Eigen::Matrix<double, 3, 4> solution;
    solution.block<3, 3>(0, 0) = R;
    solution.col(3) = tc - R * X.col(0);

    solutions(i) = solution;
  }

  return 0;
}

int P3Pnokano::solveQuartic2(const Eigen::Matrix<double, 5, 1> &factors, Eigen::Matrix<double, 4, 1> &real_roots)
{
  double A = factors(0);
  double B = factors(1);
  double C = factors(2);
  double D = factors(3);
  double E = factors(4);

  double A_pw2 = A * A;
  double B_pw2 = B * B;
  double A_pw3 = A_pw2 * A;
  double B_pw3 = B_pw2 * B;
  double A_pw4 = A_pw3 * A;
  double B_pw4 = B_pw3 * B;

  double alpha = -3 * B_pw2 / (8 * A_pw2) + C / A;
  double beta = B_pw3 / (8 * A_pw3) - B * C / (2 * A_pw2) + D / A;
  double gamma = -3 * B_pw4 / (256 * A_pw4) + B_pw2 * C / (16 * A_pw3) - B * D / (4 * A_pw2) + E / A;

  double alpha_pw2 = alpha * alpha;
  double alpha_pw3 = alpha_pw2 * alpha;

  std::complex<double> P(-alpha_pw2 / 12 - gamma, 0);
  std::complex<double> Q(-alpha_pw3 / 108 + alpha * gamma / 3 - pow(beta, 2) / 8, 0);
  std::complex<double> R = -Q / 2.0 + sqrt(pow(Q, 2.0) / 4.0 + pow(P, 3.0) / 27.0);

  std::complex<double> U = pow(R, (1.0 / 3.0));
  std::complex<double> y;

  if (U.real() == 0)
    y = -5.0 * alpha / 6.0 - pow(Q, (1.0 / 3.0));
  else
    y = -5.0 * alpha / 6.0 - P / (3.0 * U) + U;

  std::complex<double> w = sqrt(alpha + 2.0 * y);

  std::complex<double> temp;

  int i = 0;
  real_roots << 0, 0, 0, 0;
  temp = -B / (4.0 * A) + 0.5 * (w + sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
  if (abs(temp.imag()) < 1e-8 && temp.real() > 0)
  {
    real_roots(i) = temp.real();
    i++;
  }
  temp = -B / (4.0 * A) + 0.5 * (w - sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
  if (abs(temp.imag()) < 1e-8 && temp.real() > 0)
  {
    real_roots(i) = temp.real();
    i++;
  }
  temp = -B / (4.0 * A) + 0.5 * (-w + sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
  if (abs(temp.imag()) < 1e-8 && temp.real() > 0)
  {
    real_roots(i) = temp.real();
    i++;
  }
  temp = -B / (4.0 * A) + 0.5 * (-w - sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
  if (abs(temp.imag()) < 1e-8 && temp.real() > 0)
  {
    real_roots(i) = temp.real();
    i++;
  }

  return i;
}
