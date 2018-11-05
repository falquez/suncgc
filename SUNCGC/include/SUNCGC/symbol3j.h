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

#ifndef _SUNCGC_SYMBOL3J_H
#define _SUNCGC_SYMBOL3J_H

#include <array>
#include <iostream>
#include <map>
#include <vector>

namespace SUNCG::Symbol3J {
  struct NotInDecomposition : public std::exception {
    const char *what() const throw() { return "Not in decomposition"; }
  };

  template <int N>
  std::string representation_to_string(const int r);

  template <int N>
  std::vector<std::array<int, 3>> list_decompositions(int maxR1, int maxR2);

  template <int N>
  int conjugate_representation(const int r_idx);

  template <int N>
  double weight_phase(const int r_idx, const int w_idx);

  template <int N>
  int conjugate_weight(const int r_idx, const int w_idx);

  template <int N>
  int representation_dimension(const int r_idx);

  template <int N>
  std::vector<double> representation_weights(const int r_idx, const int w_idx);

  template <int N>
  double representation_casimir2(const int r_idx);
  /*
    template <int N>
    int highest_weight( const int r_idx );

    */
  template <int N>
  class Symbol6J {
  public:
    std::array<int, 6> R;
    double value;
    bool valid;
    Symbol6J(const std::array<int, 6> &R);
  };

  template <int N>
  class Coefficient {

    std::array<int, 3> R;
    std::array<int, 3> W;

    double coeff;
    bool zero;

  public:
    Coefficient(const std::array<int, 3> &rep, const std::array<int, 3> &weight, const double c,
                const bool zero = false)
        : R{rep}, W{weight}, coeff{c}, zero{zero} {}

    double value() const { return coeff; }

    const std::array<int, 3> representations() const { return std::array{R[0], R[1], R[2]}; }
    const std::array<int, 3> weights() const { return std::array{W[0], W[1], W[2]}; }

    bool isZero() const { return zero; }
  };

  template <int N>
  class Symbol3J {
    std::array<int, 3> R;
    std::array<double, 3> P;

    std::map<std::array<int, 3>, Coefficient<N>> elements;

  public:
    Symbol3J(const std::array<int, 3> &R);
    Coefficient<N> operator()(const std::array<int, 3> &W) const;
    double operator[](const std::array<int, 3> &W) const;

    std::map<std::array<int, 3>, Coefficient<N>> coefficients() const { return elements; }

    const std::array<int, 3> representations() const { return std::array{R[0], R[1], R[2]}; }

    const std::array<double, 3> phases() const { return std::array{P[0], P[1], P[2]}; }
  };

  template <int N>
  int test_phases(const int r_idx);

  template <int N>
  int test_series(const int F1);

  template <int N>
  std::ostream &operator<<(std::ostream &out, const Coefficient<N> &c);

  template <int N>
  std::ostream &operator<<(std::ostream &out, const Symbol3J<N> &c);

  template <int N>
  std::ostream &operator<<(std::ostream &out, const Symbol6J<N> &c);

} // namespace SUNCG::Symbol3J

#endif
