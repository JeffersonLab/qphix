#include "mesfield.h"

namespace QPhiX
{
//! Calculates the antihermitian field strength tensor  iF(mu,nu)
/*
 * \ingroup glue
 *
 *    F(mu,nu) =  (1/4) sum_p (1/2) [ U_p(x) - U^dag_p(x) ]
 *
 *  where
 *    U_1 = u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)
 *    U_2 = u(x,nu)*u_dag(x-mu+nu,mu)*u_dag(x-mu,nu)*u(x-mu,mu)
 *    U_3 = u_dag(x-mu,mu)*u_dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
 *    U_4 = u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u_dag(x,mu)

 *  NOTE: I am using  i*F  of the usual F defined by UKQCD, Heatlie et.al.

 * Arguments:

 *  \param f   field strength tensor f(mu,nu) (Write)
 *  \param u   gauge field (Read)
 */
template <typename U>
void mesFieldT(QDP::multi1d<U> &f, const QDP::multi1d<U> &u)
{
  f.resize(QDP::Nd * (QDP::Nd - 1) / 2);

  U tmp_0;
  U tmp_1;
  U tmp_2;
  U tmp_3;
  U tmp_4;

  QDP::Real fact = 0.125;

  int offset = 0;

  for (int mu = 0; mu < QDP::Nd - 1; ++mu) {
    for (int nu = mu + 1; nu < QDP::Nd; ++nu) {
      tmp_3 = QDP::shift(u[nu], FORWARD, mu);
      tmp_4 = QDP::shift(u[mu], FORWARD, nu);
      tmp_0 = u[nu] * tmp_4;
      tmp_1 = u[mu] * tmp_3;

      f[offset] = tmp_1 * QDP::adj(tmp_0);

      tmp_2 = QDP::adj(tmp_0) * tmp_1;
      tmp_1 = QDP::shift(tmp_2, BACKWARD, nu);
      f[offset] += QDP::shift(tmp_1, BACKWARD, mu);
      tmp_1 = tmp_4 * QDP::adj(tmp_3);
      tmp_0 = QDP::adj(u[nu]) * u[mu];

      f[offset] += QDP::shift(tmp_0 * QDP::adj(tmp_1), BACKWARD, nu);
      f[offset] += QDP::shift(QDP::adj(tmp_1) * tmp_0, BACKWARD, mu);

      tmp_0 = QDP::adj(f[offset]);
      f[offset] -= tmp_0;
      f[offset] *= fact;

      ++offset;
    }
  }
}

void mesField(QDP::multi1d<QDP::LatticeColorMatrixF> &f,
              const QDP::multi1d<QDP::LatticeColorMatrixF> &u)
{
  mesFieldT(f, u);
}

void mesField(QDP::multi1d<QDP::LatticeColorMatrixD> &f,
              const QDP::multi1d<QDP::LatticeColorMatrixD> &u)
{
  mesFieldT(f, u);
}
}
