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

#include <SUNCGC/hilbertspace.h>

namespace SUNCG::Space {

  template <int N>
  SUN_HilbertSpace<N>::SUN_HilbertSpace(const int maxR) : maxR{maxR} {
    std::vector<std::array<int, 3>> decomp;
    for (int r1 = 0; r1 < maxR; r1++) {
      for (int r2 = 0; r2 < N; r2++) {
        auto R = Symbol3J::list_decompositions<N>(r1, r2);
        decomp.insert(decomp.end(), R.begin(), R.end());
      }
    }

    for (auto const &r : decomp) {
      if (r[2] < maxR) {
        basis.push_back(Singlet<N>(r, false));
        if (r[1] == 0)
          basis.push_back(Singlet<N>(r, true));
      }
    }

    unsigned int i = 0;
    for (auto const &e : basis) {
      auto charge = e.charges();
      std::array<int, 3> c{charge[0], (e.qbar() ? -1 : charge[1]), charge[2]};
      index.emplace(c, i);
      i++;
    }
  }

  template <int N>
  unsigned int SUN_HilbertSpace<N>::operator()(const Singlet<N> &vec) const {
    unsigned int res = 0;
    auto charge = vec.charges();
    std::array<int, 3> c{charge[0], (vec.qbar() ? -1 : charge[1]), charge[2]};
    try {
      res = index.at(c);
    } catch (const std::out_of_range &e) {
      res = basis.size() + 1;
    }
    return res;
  }

  template <>
  std::ostream &operator<<<2>(std::ostream &out, const Singlet<2> &s) {
    auto [R, Q, L] = s.charges();
    auto [cR, cQ, cL] = s.casimir();

    std::string sR = R % 2 ? std::to_string(R) + "/2" : "  " + std::to_string(R / 2);
    std::string sQ = Q % 2 ? std::to_string(Q) + "/2" : "  " + std::to_string(Q / 2);
    std::string sL = L % 2 ? std::to_string(L) + "/2" : "  " + std::to_string(L / 2);

    sQ = s.qbar() ? sQ + "*" : sQ + " ";

    std::string str("[ ");
    // str += sR + "(" + std::to_string(cR) + "), ";
    // str += sQ + "(" + std::to_string(cQ) + "), ";
    // str += sL + "(" + std::to_string(cL) + ")]";
    str += sR + ", ";
    str += sQ + ", ";
    str += sL + "]";

    out << str;
    return out;
  }

  template <>
  std::ostream &operator<<<3>(std::ostream &out, const Singlet<3> &s) {
    auto [R, Q, L] = s.charges();
    auto [dimR, dimQ, dimL] = s.dimensions();

    int Rbar = Symbol3J::conjugate_representation<3>(R);
    int Lbar = Symbol3J::conjugate_representation<3>(L);

    std::string sR = std::to_string(dimR == 1 ? 0 : dimR);
    std::string sL = std::to_string(dimL == 1 ? 0 : dimL);

    sR = (R > Rbar) ? sR + "*" : sR + " ";
    sL = (L > Lbar) ? sL + "*" : sL + " ";

    std::string sQ;
    switch (Q) {
    case 0:
      sQ = s.qbar() ? "0*" : "0 ";
      break;
    case 1:
      sQ = "3 ";
      break;
    case 2:
      sQ = "3*";
      break;
    default:
      sQ = "error";
      break;
    }

    std::string str = "[ " + sR + ", " + sQ + ", " + sL + " ]";

    out << str;
    return out;
  }
} // namespace SUNCG::Space

template class SUNCG::Space::SUN_HilbertSpace<2>;
template class SUNCG::Space::SUN_HilbertSpace<3>;
