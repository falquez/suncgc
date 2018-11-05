/*
 *  Copyright (C) 2018 Carlos Falquez (falquez@nuberisim.de)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ClebschGordan.cpp"
#include <SUNCGC/symbol3j.h>

namespace SUNCG::Symbol3J {

  template <int N>
  std::vector<std::array<int, 3>> list_decompositions(const int r1, const int r2) {
    std::vector<std::array<int, 3>> R;

    clebsch::weight R1(N, r1);
    clebsch::weight R2(N, r2);
    clebsch::decomposition decomp(R1, R2);

    for (int i = 0; i < decomp.size(); ++i)
      R.push_back({R1.index(), R2.index(), conjugate_representation<N>(decomp(i).index())});

    return R;
  }

  template <int N>
  int representation_dimension(const int r_idx) {
    clebsch::weight R(N, r_idx);

    return R.dimension();
  }

  template <>
  int conjugate_representation<2>(const int r_idx) {
    return r_idx;
  }

  template <>
  int conjugate_representation<3>(const int r_idx) {

    clebsch::weight r(3, r_idx);
    clebsch::weight rc(3);

    rc(1) = r(1);
    rc(2) = r(1) - r(2);
    rc(3) = 0;

    return rc.index();
  }

  template <>
  double representation_casimir2<2>(const int r_idx) {
    clebsch::weight r(2, r_idx);
    int j2 = r(1);

    double c2 = static_cast<double>(j2 * (j2 + 2)) / 4;

    return c2;
  }
  template <>
  double representation_casimir2<3>(const int r_idx) {
    clebsch::weight r(3, r_idx);
    int p = r(1) - r(2);
    int q = r(2);

    double c2 = static_cast<double>(p * p + p * q + q * q + 3 * p + 3 * q) / 3;

    return c2;
  }

  template <>
  double weight_phase<2>(const int r_idx, const int w_idx) {
    clebsch::weight R(2, r_idx);
    clebsch::pattern M(R, w_idx);

    double m = double(2 * M(1, 1) - M(1, 2)) / 2;

    return m;
  }

  template <>
  double weight_phase<3>(const int r_idx, const int w_idx) {
    clebsch::weight R(3, r_idx);
    clebsch::pattern M(R, w_idx);

    double I = double(M(1, 2) - M(2, 2)) / 2;
    double Iz = M(1, 1) - M(2, 2) - I;
    double Y = 2 * (M(1, 1) - Iz - double(M(1, 3) + M(2, 3)) / 3);

    double phase = (Iz + 1.5 * Y);

    return phase;
  }

  template <>
  std::vector<double> representation_weights<2>(const int r_idx, const int w_idx) {
    int w = 2 * w_idx - r_idx;
    return std::vector{0.5 * r_idx, 0.5 * w};
  }

  template <>
  std::vector<double> representation_weights<3>(const int r_idx, const int w_idx) {
    clebsch::weight R(3, r_idx);
    clebsch::pattern M(R, w_idx);

    double I = double(M(1, 2) - M(2, 2)) / 2;
    double Iz = M(1, 1) - M(2, 2) - I;
    double Y = 2 * (M(1, 1) - Iz - double(M(1, 3) + M(2, 3)) / 3);

    return std::vector{I, Iz, Y};
  }

  /*
    |j  m> =  (2j,0)(j+m) = (w_12,0)(w_11)
    |j -m> -> (2j,0)(j-m) = (w_12,0)(w_12-w_11)
  */
  template <>
  int conjugate_weight<2>(const int r_idx, const int w_idx) {

    clebsch::weight r(2, r_idx);
    clebsch::pattern w(r, w_idx);
    clebsch::pattern wbar(2);

    wbar(1, 2) = w(1, 2);
    wbar(2, 2) = 0;
    wbar(1, 1) = w(1, 2) - w(1, 1);

    return wbar.index();
  }

  /*
    |(p,q) I Iz Y> -> |(q,p) I -Iz -Y>
    (p+q,q,0)(I + Y/2 + (p + 2q)/3, -I + Y/2 + (p + 2q)/3 )(Iz + Y/2 + (p + 2q)/3)
    ->  (p+q,p,0)(I - Y/2 + (q + 2p)/3, -I + -Y/2 + (q + 2p)/3 )(-Iz + -Y/2 + (q + 2p)/3)
  */

  template <>
  int conjugate_weight<3>(const int r_idx, const int w_idx) {

    clebsch::weight r(3, r_idx);
    clebsch::pattern w(r, w_idx);
    clebsch::pattern wbar(3);

    wbar(1, 3) = w(1, 3);
    wbar(2, 3) = w(1, 3) - w(2, 3);
    wbar(3, 3) = 0;
    wbar(1, 2) = w(1, 3) - w(2, 2);
    wbar(2, 2) = w(1, 3) - w(1, 2);
    wbar(1, 1) = w(1, 3) - w(1, 1);

    return wbar.index();
  }

  template <>
  std::string representation_to_string<3>(const int r) {
    std::string s;
    int dimR = representation_dimension<3>(r);
    int Rbar = conjugate_representation<3>(r);
    s = std::to_string(dimR == 1 ? 0 : dimR);
    s = (r > Rbar) ? s + "*" : s + " ";
    return s;
  }

  template <int N>
  Symbol3J<N>::Symbol3J(const std::array<int, 3> &R) : R{R} {

    clebsch::weight S1(N, R[0]);
    clebsch::weight S2(N, R[1]);
    clebsch::weight S3(N, R[2]);
    clebsch::weight S3bar(N, conjugate_representation<N>(R[2]));

    clebsch::decomposition decomp(S1, S2);

    bool present = false;
    for (int i = 0; i < decomp.size(); ++i)
      if (decomp(i) == S3bar && (present = true))
        break;

    if (!present)
      throw NotInDecomposition();

    int dim1 = S1.dimension();
    int dim2 = S2.dimension();
    int dim3 = S3.dimension();

    double hw_phase1 = weight_phase<N>(S1.index(), S1.dimension() - 1);
    double hw_phase2 = weight_phase<N>(S2.index(), S2.dimension() - 1);
    double hw_phase3 = weight_phase<N>(S3.index(), S3.dimension() - 1);

    P = {hw_phase1, hw_phase2, hw_phase3};

    // TODO: multiplicity
    int alpha = 0;
    const clebsch::coefficients C(S3bar, S1, S2);
    for (int m3 = 0; m3 < dim3; ++m3) {
      for (int m1 = 0; m1 < dim1; ++m1) {
        for (int m2 = 0; m2 < dim2; ++m2) {

          int m3bar = conjugate_weight<N>(S3.index(), m3);
          double x = double(C(m1, m2, alpha, m3bar));

          double w_phase = weight_phase<N>(S3.index(), m3);
          double w_phase2 = weight_phase<N>(S3bar.index(), m3bar);

          double phase = hw_phase1 - hw_phase2 - w_phase; //( IzH + 3*YH/2) - ( Iz + 3*Y/2);
          // double phase = hw_phase3 - w_phase;//( IzH + 3*YH/2) - ( Iz + 3*Y/2);

          int iphase = std::abs(std::round(phase));
          double signp = iphase % 2 ? -1.0 : 1.0;
          double f = signp / sqrt((double)dim3);

          // std::cout << "(" << m3 << " -> " << phase << "," << iphase << ")";

          double cg = x * f;

          if (fabs(x) > clebsch::EPS) {
            elements.emplace(std::array{m1, m2, m3}, Coefficient<N>(R, {m1, m2, m3}, cg));
            // elements.insert( std::pair<std::array<int, 3>, Coefficient> { {m1,m2,m3}, Coefficient(N, R, {m1,m2,m3},
            // cg) });
          }
        }
      }
    }
  }

  template <int N>
  double Symbol3J<N>::Symbol3J::operator[](const std::array<int, 3> &W) const {
    auto search = elements.find(W);
    if (search != elements.end()) {
      return search->second.value();
    } else {
      return 0.0;
    }
  }
  template <>
  Symbol6J<3>::Symbol6J(const std::array<int, 6> &R) : R{R}, value{0}, valid{false} {
    std::vector<int> J(7);
    std::vector<int> Jbar(7);
    std::vector<int> dimJ(7);
    std::vector<double> phase_h(7);
    for (int i = 0; i < 6; i++) {
      J[i + 1] = R[i];
      Jbar[i + 1] = conjugate_representation<3>(R[i]);
      dimJ[i + 1] = representation_dimension<3>(R[i]);
      phase_h[i + 1] = weight_phase<3>(J[i + 1], dimJ[i + 1] - 1);
    }

    try {
      Symbol3J<3> a1({J[1], J[2], J[3]});
      Symbol3J<3> a2({J[5], J[1], Jbar[6]});
      Symbol3J<3> a3({J[6], J[2], Jbar[4]});
      Symbol3J<3> a4({J[4], J[3], Jbar[5]});
      std::vector<int> m(7);
      for (const auto &[c, e] : a1.coefficients()) {
        m[1] = c[0];
        m[2] = c[1];
        m[3] = c[2];
        for (m[4] = 0; m[4] < dimJ[4]; m[4]++) {
          for (m[5] = 0; m[5] < dimJ[5]; m[5]++) {
            for (m[6] = 0; m[6] < dimJ[6]; m[6]++) {

              int m4_bar = conjugate_weight<3>(J[4], m[4]);
              int m5_bar = conjugate_weight<3>(J[5], m[5]);
              int m6_bar = conjugate_weight<3>(J[6], m[6]);

              double m4_phase = weight_phase<3>(J[4], m[4]);
              double m5_phase = weight_phase<3>(J[5], m[5]);
              double m6_phase = weight_phase<3>(J[6], m[6]);

              // double phase = h1_phase-w1_phase+ h2_phase-w2_phase+ h3_phase-w3_phase;//( IzH + 3*YH/2) - ( Iz +
              // 3*Y/2);
              double phase = phase_h[1] + phase_h[2] + phase_h[3] + phase_h[4] - m4_phase + phase_h[5] - m5_phase +
                             phase_h[6] - m6_phase;

              int iphase = std::abs(std::round(phase));
              double signp = iphase % 2 ? -1.0 : 1.0;

              double c1 = a1[{m[1], m[2], m[3]}];
              double c2 = a2[{m[5], m[1], m6_bar}];
              double c3 = a3[{m[6], m[2], m4_bar}];
              double c4 = a4[{m[4], m[3], m5_bar}];
              value += signp * c1 * c2 * c3 * c4;
            }
          }
        }
      }
      valid = true;
    } catch (const std::exception &e) {
    }
  }

  template <>
  std::ostream &operator<<(std::ostream &out, const Coefficient<3> &c) {
    auto weights = c.weights();

    out << "(" << weights[0] << "," << weights[1] << "," << weights[2] << ")";
    out << "(" << weights[3] << "," << weights[4] << "," << weights[5] << ")=";

    out << c.value(); // << std::endl;

    return out;
  }

  template <>
  std::ostream &operator<<(std::ostream &out, const Coefficient<2> &c) {
    std::vector<std::string> rs(3), ws(3);
    auto reps = c.representations();
    for (int i = 0; i < 3; i++) {
      rs[i] = reps[i] % 2 ? std::to_string(reps[i]) + "/2" : " " + std::to_string(reps[i] / 2) + " ";
    }

    auto weights = c.weights();
    for (int i = 0; i < 3; i++) {
      int w = 2 * weights[i] - reps[i];
      ws[i] = w % 2 ? std::to_string(w) + "/2" : " " + std::to_string(w / 2) + " ";
      ws[i] = w < 0 ? ws[i] : " " + ws[i];
    }

    out << "[(" << rs[0] << "," << rs[1] << "," << rs[2] << ")";
    out << " (" << ws[0] << "," << ws[1] << "," << ws[2] << ")]";
    out << "=" << c.value(); // << std::endl;

    return out;
  }

  template <>
  std::ostream &operator<<(std::ostream &out, const Symbol3J<2> &s) {
    std::vector<std::string> rs(3);
    auto reps = s.representations();
    for (int i = 0; i < 3; i++) {
      rs[i] = reps[i] % 2 ? std::to_string(reps[i]) + "/2" : " " + std::to_string(reps[i] / 2) + " ";
    }
    auto ps = s.phases();

    out << "Symbol3J<2>(" << rs[0] << "," << rs[1] << "," << rs[2] << ")";
    out << "\t(" << ps[0] << "," << ps[1] << "," << ps[2] << ")" << std::endl;

    for (const auto &[c, e] : s.coefficients()) {
      std::vector<std::string> ws(3);

      auto weights = e.weights();
      for (int i = 0; i < 3; i++) {
        int w = 2 * weights[i] - reps[i];
        ws[i] = w % 2 ? std::to_string(w) + "/2" : " " + std::to_string(w / 2) + " ";
        ws[i] = w < 0 ? ws[i] : " " + ws[i];
      }

      out << "\t(" << ws[0] << "," << ws[1] << "," << ws[2] << ")";
      out << "\t=" << e.value() << std::endl;
    }

    return out;
  }

  template <>
  std::ostream &operator<<(std::ostream &out, const Symbol3J<3> &s) {
    std::vector<std::string> rs(3);
    auto reps = s.representations();
    for (int i = 0; i < 3; i++) {
      int dimR = representation_dimension<3>(reps[i]);
      int Rbar = conjugate_representation<3>(reps[i]);
      rs[i] = std::to_string(dimR == 1 ? 0 : dimR);
      rs[i] = (reps[i] > Rbar) ? rs[i] + "*" : rs[i] + " ";
    }
    auto ps = s.phases();

    out << "Symbol3J<3>(" << rs[0] << "," << rs[1] << "," << rs[2] << ")";
    out << "\t(" << ps[0] << "," << ps[1] << "," << ps[2] << ") ->" << ps[0] + ps[1] + ps[2] << std::endl;

    for (const auto &[c, e] : s.coefficients()) {
      std::vector<std::string> ws(3);
      auto weights = e.weights();
      for (int i = 0; i < 3; i++) {
        std::vector<double> w = representation_weights<3>(reps[i], weights[i]);
        ws[i] += "(Idx=" + std::to_string(weights[i] + 1);
        ws[i] += ",Iz=" + std::to_string(w[1]);
        ws[i] += ",Y=" + std::to_string(w[2]) + ")";
        // int w = 2*weights[i] -reps[i];
        // ws[i] = w%2 ?  std::to_string(w) + "/2" : " " + std::to_string(w/2) + " ";
        // ws[i] = w < 0 ? ws[i] : " " + ws[i];
      }

      out << "\t[" << ws[0] << "," << ws[1] << "," << ws[2] << "]";
      out << "\t=" << e.value() << std::endl;
    }

    return out;
  }

  template <>
  std::ostream &operator<<(std::ostream &out, const Symbol6J<3> &s) {
    std::vector<std::string> rs(6);
    for (int i = 0; i < 6; i++) {
      int dimR = representation_dimension<3>(s.R[i]);
      int Rbar = conjugate_representation<3>(s.R[i]);
      rs[i] = std::to_string(dimR == 1 ? 0 : dimR);
      rs[i] = (s.R[i] > Rbar) ? rs[i] + "*" : rs[i] + " ";
    }
    out << "Symbol6J(" << rs[0] << "," << rs[1] << "," << rs[2] << "," << rs[3] << "," << rs[4] << "," << rs[5] << ")";
    out << " = " << s.value << std::endl;
    return out;
  }

  template <>
  int test_phases<3>(const int r_idx) {
    int r_bar = conjugate_representation<3>(r_idx);
    clebsch::weight S1(3, r_idx);
    clebsch::weight S1bar(3, r_bar);
    clebsch::weight S0(3, 0);

    int dim1 = S1.dimension();
    std::string r_str = std::to_string(dim1);
    std::string r_str1 = (r_bar > r_idx) ? r_str + " " : r_str + "*";
    std::string r_str2 = (r_bar > r_idx) ? r_str + "*" : r_str + " ";

    std::vector<double> highest_weight = representation_weights<3>(S1.index(), S1.dimension() - 1);

    std::cout << "Symbol<3>( " << r_str1 << " " << r_str2 << " 0 )";
    std::cout << "HW=(Iz=" + std::to_string(highest_weight[1]);
    std::cout << ",Y=" + std::to_string(highest_weight[2]) + ")" << std::endl;

    // TODO: multiplicity
    int alpha = 0;
    const clebsch::coefficients C(S0, S1, S1bar);
    // for (int m3 = 0; m3 < dim3; ++m3) {
    for (int m1 = 0; m1 < dim1; ++m1) {

      int m3 = 0;
      int m1_bar = conjugate_weight<3>(S1.index(), m1);
      double x = double(C(m1, m1_bar, alpha, m3));

      std::vector<std::string> ws(2);
      std::vector<double> w1, w2;
      w1 = representation_weights<3>(r_idx, m1);
      ws[0] += "(Iz=" + std::to_string(w1[1]);
      ws[0] += ",Y=" + std::to_string(w1[2]) + ")";
      w2 = representation_weights<3>(r_bar, m1_bar);
      ws[1] += "(Iz=" + std::to_string(w2[1]);
      ws[1] += ",Y=" + std::to_string(w2[2]) + ")";

      double phase0 = highest_weight[1] - w1[1] + 0.5 * (highest_weight[2] - w1[2]);
      double phase = weight_phase<3>(r_idx, dim1 - 1) - weight_phase<3>(r_idx, m1);
      std::cout << "\t[" << ws[0] << "," << ws[1] << "]";
      std::cout << "\t=" << x * std::sqrt(dim1) << ", p0=" << phase0 << ", p=" << phase << std::endl;
    }

    return 0;
  }

  template <>
  int test_series<3>(const int F1) {

    std::cout << "Testing with " << representation_to_string<3>(F1) << std::endl;

    std::vector<std::array<int, 3>> r1 = list_decompositions<3>(F1, 1);
    std::vector<std::array<int, 3>> r2 = list_decompositions<3>(F1, 2);

    for (int i1 = 0; i1 < r1.size(); i1++) {
      for (int i2 = 0; i2 < r2.size(); i2++) {
        const int E1 = r1[i1][2];
        const int E1bar = conjugate_representation<3>(E1);
        const int Fbar = r2[i2][2];
        const int F = conjugate_representation<3>(Fbar);
        const int dimE1 = representation_dimension<3>(E1);

        Symbol6J<3> s6j({F1, 2, Fbar, 1, E1, 2});

        try {
          SUNCG::Symbol3J::Symbol3J<3> a1({E1, F1, 1});
          SUNCG::Symbol3J::Symbol3J<3> a2({2, 2, 2});
          SUNCG::Symbol3J::Symbol3J<3> a3({1, Fbar, E1bar});
          SUNCG::Symbol3J::Symbol3J<3> a0({F1, 2, Fbar});

          std::cout << "(E1" << representation_to_string<3>(E1) << ",F1=" << representation_to_string<3>(F1) << ", 3)";
          std::cout << "(3*,3*,3*)";
          std::cout << "(3,Fbar=" << representation_to_string<3>(Fbar)
                    << ",E1bar=" << representation_to_string<3>(E1bar) << ")";
          std::cout << "\t[E1=" << E1 << ",F1=" << F1 << ",Fbar=" << Fbar << "]" << std::endl;

          std::cout << a0 << std::endl;
          for (const auto &[c, e] : a0.coefficients()) {
            double r = 0;
            // std::cout << "(" << c[0] << "," << c[1] << "," << c[2] << ")" << std::endl;

            for (int m1 = 0; m1 < dimE1; m1++) {
              int m1_bar = conjugate_weight<3>(E1, m1);
              for (int q1 = 0; q1 < 3; q1++) {
                int q1_bar = conjugate_weight<3>(1, q1);
                for (int a = 0; a < 3; a++) {
                  int a_bar = conjugate_weight<3>(1, a);

                  double h1_phase = weight_phase<3>(E1, dimE1 - 1);
                  double w1_phase = weight_phase<3>(E1, m1);
                  double h2_phase = weight_phase<3>(2, 2);
                  double w2_phase = weight_phase<3>(2, q1_bar);
                  double h3_phase = weight_phase<3>(1, 2);
                  double w3_phase = weight_phase<3>(1, a);

                  double phase = h1_phase - w1_phase + h2_phase - w2_phase + h3_phase -
                                 w3_phase; //( IzH + 3*YH/2) - ( Iz + 3*Y/2);
                  // double phase = h3_phase-w3_phase;//( IzH + 3*YH/2) - ( Iz + 3*Y/2);
                  // double phase = hw_phase3 - w_phase;//( IzH + 3*YH/2) - ( Iz + 3*Y/2);

                  int iphase = std::abs(std::round(phase));
                  double signp = iphase % 2 ? -1.0 : 1.0;
                  double c1 = a1[{m1, c[0], q1}];
                  double c2 = a2[{q1_bar, c[1], a_bar}];
                  double c3 = a3[{a, c[2], m1_bar}];
                  r += signp * c1 * c2 * c3;
                  // if(c1*c2*c3 > 1E-9){
                  // std::cout << "\t (" << m1+1 << "," << q1+1 << "," << a+1 << ") -> ";
                  // std::cout << "c1=" << c1 << ", c2=" << c2 << ", c3=" << c3;
                  // std::cout << " p=" << signp << " -> " << c1*c2*c3 << " -> "<<r << std::endl;

                  //}
                }
              }
            }

            std::cout << "\t\t r=" << r << ", e=" << e.value() << " e/r=" << e.value() * s6j.value / r;
            std::cout << std::endl;
          }
        } catch (const std::exception &e) {
          std::cout << "error" << std::endl;
          continue;
        }
      }
    }

    return 0;
  }

} // namespace SUNCG::Symbol3J

template class SUNCG::Symbol3J::Coefficient<2>;
template class SUNCG::Symbol3J::Coefficient<3>;

template class SUNCG::Symbol3J::Symbol3J<2>;
template class SUNCG::Symbol3J::Symbol3J<3>;

template int SUNCG::Symbol3J::representation_dimension<2>(const int);
template int SUNCG::Symbol3J::representation_dimension<3>(const int);

template std::vector<std::array<int, 3>> SUNCG::Symbol3J::list_decompositions<2>(const int, const int);
template std::vector<std::array<int, 3>> SUNCG::Symbol3J::list_decompositions<3>(const int, const int);
