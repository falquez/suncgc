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

#ifndef _SUNCGC_HILBERTSPACE_H
#define _SUNCGC_HILBERTSPACE_H

#include <array>
#include <vector>

#include <SUNCGC/symbol3j.h>
#include <SUNCGC/vector.h>

namespace SUNCG::Space {

  template <int N>
  class Singlet {
    std::array<int, 3> charge;
    std::array<int, 3> dimension;
    std::array<double, 3> casimir2;
    bool _qbar;

  public:
    Singlet(const std::array<int, 3> &charge, const bool qbar = false) : charge{charge}, _qbar{qbar} {
      for (int i = 0; i < 3; i++) {
        dimension[i] = Symbol3J::representation_dimension<N>(charge[i]);
        casimir2[i] = Symbol3J::representation_casimir2<N>(charge[i]);
      }
    }

    bool operator<(const Singlet &vec) const {
      std::array<int, 3> c1{charge[0], (_qbar ? -1 : charge[1]), charge[2]};
      std::array<int, 3> c2{vec.charge[0], (vec._qbar ? -1 : vec.charge[1]), vec.charge[2]};

      return c1 < c2;
    }
    std::array<int, 3> charges() const { return charge; }
    int operator[](const int idx) const { return charge[idx]; }
    bool qbar() const { return _qbar; }
    std::array<double, 3> casimir() const { return casimir2; }
    std::array<int, 3> dimensions() const { return dimension; }
    // virtual std::string toString() const = 0;
  };

  template <int N>
  std::ostream &operator<<(std::ostream &out, const Singlet<N> &s);

  template <int N>
  class SUN_HilbertSpace : public HilbertSpace<Singlet<N>> {
    int maxR;

    std::map<std::array<int, 3>, unsigned int> index;
    std::vector<Singlet<N>> basis;

  public:
    SUN_HilbertSpace(const int maxR);

    const Singlet<N> &operator[](const unsigned int i) const override { return basis.at(i); }

    unsigned int operator()(const Singlet<N> &vec) const override;

    unsigned int dimension() const { return basis.size(); }
  };
} // namespace SUNCG::Space

#endif // _SUNCGC_HILBERTSPACE_H
