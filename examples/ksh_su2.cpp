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
#include <SUNCGC/vector.h>

#include <TNT/storage/storage.h>
#include <TNT/tensor/sparse/tensor.h>
#include <TNT/tensor/tensor.h>

namespace SUNCG::Operator {
  using namespace SUNCG::Space;

  template <typename F>
  class operatorE : public VectorOperator<Singlet<2>, F> {
  public:
    operatorE(const HilbertSpace<Singlet<2>> &space) : VectorOperator<Singlet<2>, F>{space} {};
    std::vector<Vector<F>> operator()(const Vector<F> &ket) const {
      std::vector<Vector<F>> result;
      auto C2 = this->space[ket[0]].casimir()[2];
      if (C2 > 1E-12)
        result.push_back(Vector<F>({ket[0]}, C2));

      return result;
    };
    std::string name() const { return std::string("operatorE"); };
  };

  template <typename F>
  class operatorF : public VectorOperator<Singlet<2>, F> {
  public:
    operatorF(const HilbertSpace<Singlet<2>> &space) : VectorOperator<Singlet<2>, F>{space} {};
    std::vector<Vector<F>> operator()(const Vector<F> &ket) const {
      std::vector<Vector<F>> result;
      auto C2 = this->space[ket[0]].casimir()[0];
      if (C2 > 1E-12)
        result.push_back(Vector<F>({ket[0]}, C2));

      return result;
    };
    std::string name() const { return std::string("operatorF"); };
  };

  template <typename F>
  class operatorN : public VectorOperator<Singlet<2>, F> {
  public:
    operatorN(const HilbertSpace<Singlet<2>> &space) : VectorOperator<Singlet<2>, F>{space} {};
    std::vector<Vector<F>> operator()(const Vector<F> &ket) const {
      std::vector<Vector<F>> result;
      auto Q = this->space[ket[0]].qbar() ? 2 : this->space[ket[0]].charges()[1];
      if (Q != 0)
        result.push_back(Vector<F>({ket[0]}, Q));
      return result;
    };
    std::string name() const { return std::string("operatorN"); };
  };

  template <typename F>
  class operatorW2 : public VectorOperator<Singlet<2>, F> {
  public:
    operatorW2(const HilbertSpace<Singlet<2>> &space) : VectorOperator<Singlet<2>, F>{space} {};
    std::vector<Vector<F>> operator()(const Vector<F> &ket) const {
      constexpr unsigned int N = 2;
      std::vector<Vector<F>> result;

      const unsigned int dimH = this->space.dimension();
      const int R = 1;

      // Read SUN singlets from tensor index
      const auto u1 = this->space[ket[0]];
      const auto u2 = this->space[ket[1]];

      // Read the charges
      const auto [F1, Q1, E1] = u1.charges();
      const auto [F2, Q2, E2] = u2.charges();

      const int E1bar = Symbol3J::conjugate_representation<N>(E1);
      const int E2bar = Symbol3J::conjugate_representation<N>(E2);

      if (E1bar == F2 && !u1.qbar() && (Q2 != 0 || u2.qbar())) {
        const int Q1p = (Q1 + 1) % N;
        const int Q2p = (Q2 + N - 1) % N;

        for (const auto [a, b, Fpbar] : Symbol3J::list_decompositions<N>(R, F2)) {
          const int Fp = Symbol3J::conjugate_representation<N>(Fpbar);

          const auto s1 = Singlet<N>({F1, Q1p, Fpbar}, Q1p == 0);
          const auto s2 = Singlet<N>({Fp, Q2p, E2});

          const unsigned int i1 = this->space(s1);
          const unsigned int i2 = this->space(s2);

          if (i1 < dimH || i2 < dimH) {
            const int dimFp = Symbol3J::representation_dimension<N>(Fp);
            const int dimF1 = Symbol3J::representation_dimension<N>(F1);
            const int dimF2 = Symbol3J::representation_dimension<N>(F2);
            const int dimE2 = Symbol3J::representation_dimension<N>(E2);
            const int dimQ1 = Symbol3J::representation_dimension<N>(Q1);
            const int dimQ2p = Symbol3J::representation_dimension<N>(Q2p);

            double phase_F1 = Symbol3J::weight_phase<N>(F1, dimF1 - 1);
            double phase_F2 = Symbol3J::weight_phase<N>(F2, dimF2 - 1);
            double phase_Fp = Symbol3J::weight_phase<N>(Fp, dimFp - 1);
            double phase_E2 = Symbol3J::weight_phase<N>(E2, dimE2 - 1);
            double phase_Q1 = Symbol3J::weight_phase<N>(Q1, dimQ1 - 1);
            double phase_Q2p = Symbol3J::weight_phase<N>(Q2p, dimQ2p - 1);

            double phase = phase_E2 - phase_F1 + phase_F2 - phase_Fp + phase_Q2p - phase_Q1;

            int iphase = std::abs(std::round(phase));
            double p = iphase % 2 ? -1.0 : 1.0;
            double c = std::sqrt(dimFp) * std::sqrt(dimF2) / (std::sqrt(dimF1) * std::sqrt(dimE2));

            result.push_back(Vector<F>({i1, i2}, p * c));
          } else {
            // std::cout << "out of bounds: " << u1 << u2 << "->" << s1 << s2 <<
            // std::endl;
            continue;
          }
        }
      }

      return result;
    }
    std::string name() const { return std::string("operatorW2"); }
  };

  template <typename F>
  class operatorP2 : public VectorOperator<Singlet<2>, F> {
  public:
    operatorP2(const HilbertSpace<Singlet<2>> &space) : VectorOperator<Singlet<2>, F>{space} {}
    std::vector<Vector<F>> operator()(const Vector<F> &ket) const {
      std::vector<Vector<F>> result;

      // Read SUN singlets from tensor index
      const auto u1 = this->space[ket[0]];
      const auto u2 = this->space[ket[1]];

      // Read the charges
      const auto [F1, Q1, E1] = u1.charges();
      const auto [F2, Q2, E2] = u2.charges();

      if (E1 == F2)
        result.push_back(Vector<F>({ket[0], ket[1]}, 1.0));

      return result;
    }
    std::string name() const { return std::string("operatorP2"); }
  };

} // namespace SUNCG::Operator

int main(int argc, char **argv) {
  using namespace TNT;

  int err = 0;

  if (argc != 2) {
    std::cout << argv[0] << " <maxR>";
    exit(0);
  }

  const int maxR = std::atoi(argv[1]);

  auto H = SUNCG::Space::SUN_HilbertSpace<2>(maxR);
  auto opN = SUNCG::Operator::operatorN<double>(H);
  auto opE = SUNCG::Operator::operatorE<double>(H);
  auto opF = SUNCG::Operator::operatorF<double>(H);
  auto opW = SUNCG::Operator::operatorW2<double>(H);
  auto opP = SUNCG::Operator::operatorP2<double>(H);

  auto dimH = H.dimension();

  std::cout << "dimH=" << dimH << std::endl;
  for (int i = 0; i < dimH; i++) {
    std::cout << H[i] << std::endl;
  }

  Tensor::Tensor<double> TN({dimH, dimH});
  Tensor::Tensor<double> TE({dimH, dimH});
  Tensor::Tensor<double> TF({dimH, dimH});
  Tensor::Tensor<double> TW({dimH, dimH, dimH, dimH});
  Tensor::Tensor<double> TP({dimH, dimH, dimH, dimH});

  for (unsigned int i1 = 0; i1 < dimH; i1++) {
    for (unsigned int i2 = 0; i2 < dimH; i2++) {
      SUNCG::Space::Vector<double> v1({i1}), v2({i2});
      if (auto v = v1 * opN(v2))
        TN[{i1, i2}] = *v;
      if (auto v = v1 * opE(v2))
        TE[{i1, i2}] = *v;
      if (auto v = v1 * opF(v2))
        TF[{i1, i2}] = *v;
    }
  }

  for (unsigned int i1 = 0; i1 < dimH; i1++) {
    for (unsigned int j1 = 0; j1 < dimH; j1++) {
      for (unsigned int i2 = 0; i2 < dimH; i2++) {
        for (unsigned int j2 = 0; j2 < dimH; j2++) {
          SUNCG::Space::Vector<double> v1({i1, j1}), v2({i2, j2});
          if (auto v = v1 * opW(v2))
            TW[{i1, j1, i2, j2}] = *v;
          if (auto v = v1 * opP(v2))
            TP[{i1, j1, i2, j2}] = *v;
        }
      }
    }
  }

  std::cout << TN << std::endl;
  std::cout << TE << std::endl;
  std::cout << TF << std::endl;

  auto W2 = TW + TW.conjugate().transpose({2, 3, 0, 1});
  // M = AxB
  // {1,_0_},{3,_2_} order is defined by T.matricize({_0_,_2_},{1,3})
  auto O1 = W2.matricize({1, 0}, {3, 2});

  std::cout << "W=" << std::endl;
  for (unsigned int i1 = 0; i1 < dimH; i1++) {
    for (unsigned int j1 = 0; j1 < dimH; j1++) {
      for (unsigned int i2 = 0; i2 < dimH; i2++) {
        for (unsigned int j2 = 0; j2 < dimH; j2++) {
          double v = W2[{i1, j1, i2, j2}];
          if (std::abs(v) > 1E-10) {
            std::cout << "(" << i1 + 1 << "," << j1 + 1 << "," << i2 + 1 << "," << j2 + 1 << ") = ";
            std::cout << H[i1] << H[j1] << "<-" << H[i2] << H[j2] << "=" << v << std::endl;
          }
        }
      }
    }
  }
  std::cout << "P=" << std::endl;
  for (unsigned int i1 = 0; i1 < dimH; i1++) {
    for (unsigned int j1 = 0; j1 < dimH; j1++) {
      for (unsigned int i2 = 0; i2 < dimH; i2++) {
        for (unsigned int j2 = 0; j2 < dimH; j2++) {
          double v = TP[{i1, j1, i2, j2}];
          if (std::abs(v) > 1E-10) {
            std::cout << "(" << i1 + 1 << "," << j1 + 1 << "," << i2 + 1 << "," << j2 + 1 << ") = ";
            std::cout << H[i1] << H[j1] << "<-" << H[i2] << H[j2] << "=" << v << std::endl;
          }
        }
      }
    }
  }

  std::cout << "W2=" << W2 << std::endl;
  std::cout << "O1=" << O1 << std::endl;

  for (unsigned int i = 0; i < dimH * dimH; i++) {
    for (unsigned int j = 0; j < dimH * dimH; j++) {
      if (std::abs(O1[{i, j}]) > 1E-10)
        std::cout << "{" << i + 1 << "," << j + 1 << "}->" << O1[{i, j}] << ",";
    }
  }
  std::cout << std::endl;

  auto uv = Tensor::kronecker_SVD(W2);

  auto lr = Tensor::kronecker_SVD(TP);

  std::cout << std::endl;
  for (int n = 0; n < uv.size(); n++) {
    std::cout << "U" << n << "={";
    for (unsigned int i = 0; i < dimH; i++) {
      for (unsigned int j = 0; j < dimH; j++) {
        if (std::abs(uv[n][0][{i, j}]) > 1E-10)
          std::cout << "{" << i + 1 << "," << j + 1 << "}->" << uv[n][0][{i, j}] << ",";
      }
    }
    std::cout << "}" << std::endl;
  }
  for (int n = 0; n < uv.size(); n++) {
    std::cout << "V" << n << "={";
    for (unsigned int i = 0; i < dimH; i++) {
      for (unsigned int j = 0; j < dimH; j++) {
        if (std::abs(uv[n][1][{i, j}]) > 1E-10)
          std::cout << "{" << i + 1 << "," << j + 1 << "}->" << uv[n][1][{i, j}] << ",";
      }
    }
    std::cout << "}" << std::endl;
  }

  std::cout << std::endl;
  for (int n = 0; n < lr.size(); n++) {
    std::cout << "L" << n << "={";
    for (unsigned int i = 0; i < dimH; i++) {
      for (unsigned int j = 0; j < dimH; j++) {
        if (std::abs(lr[n][0][{i, j}]) > 1E-10)
          std::cout << "{" << i + 1 << "," << j + 1 << "}->" << lr[n][0][{i, j}] << ",";
      }
    }
    std::cout << "}" << std::endl;
  }
  for (int n = 0; n < lr.size(); n++) {
    std::cout << "R" << n << "={";
    for (unsigned int i = 0; i < dimH; i++) {
      for (unsigned int j = 0; j < dimH; j++) {
        if (std::abs(lr[n][1][{i, j}]) > 1E-10)
          std::cout << "{" << i + 1 << "," << j + 1 << "}->" << lr[n][1][{i, j}] << ",";
      }
    }
    std::cout << "}" << std::endl;
  }

  std::string filename = "hamiltonian2_" + std::to_string(maxR) + ".hdf5";
  Storage::Storage storage(filename, Storage::FileMode::CreateOverwrite);

  storage.create_group("/HilbertSpace");
  storage.create("/HilbertSpace", Storage::Data::Metadata<unsigned int>{"dimension", dimH});
  {
    storage.create_group("/HilbertSpace/Operator/E2");
    Storage::Data::Dense<double> dense{TE.dimension(), TE.size()};
    TE.writeTo(dense.data.get());
    storage.create("/HilbertSpace/Operator/E2/0", dense);
  }
  {
    storage.create_group("/HilbertSpace/Operator/F2");
    Storage::Data::Dense<double> dense{TF.dimension(), TF.size()};
    TF.writeTo(dense.data.get());
    storage.create("/HilbertSpace/Operator/F2/0", dense);
  }
  {
    storage.create_group("/HilbertSpace/Operator/N");
    Storage::Data::Dense<double> dense{TN.dimension(), TN.size()};
    TF.writeTo(dense.data.get());
    storage.create("/HilbertSpace/Operator/N/0", dense);
  }
  storage.create_group("/HilbertSpace/Operator/U");
  for (unsigned int n = 0; n < uv.size(); n++) {
    Storage::Data::Dense<double> dense{uv[n][0].dimension(), uv[n][0].size()};
    uv[n][0].writeTo(dense.data.get());
    storage.create("/HilbertSpace/Operator/U/" + std::to_string(n), dense);
  }
  storage.create_group("/HilbertSpace/Operator/V");
  for (unsigned int n = 0; n < uv.size(); n++) {
    Storage::Data::Dense<double> dense{uv[n][1].dimension(), uv[n][1].size()};
    uv[n][1].writeTo(dense.data.get());
    storage.create("/HilbertSpace/Operator/V/" + std::to_string(n), dense);
  }
}
