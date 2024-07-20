// A Ptyhia8  test program to produce secondary X(3872)s from B decays
/*
  author : Serhat Istin
          istin@cern.ch
*/
// #define DEBUG_GEN

#include <algorithm>
#include <fmt/ranges.h>
#include <fmt/std.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <variant>

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <Math/Vector4D.h>

using namespace Pythia8;
using FourMom = ROOT::Math::PtEtaPhiEVector;

class HistogramRegistry {
  using Histogram = std::variant<TH1F, TH2F>;
  std::unordered_map<std::string, Histogram> m_histo_list;

public:
  template <typename K, typename... T> void Book(K key, T &&...args) {
    if constexpr (sizeof...(args) == 4) {
      m_histo_list[key] = TH1F(key, std::forward<T>(args)...);
    } else if constexpr (sizeof...(args) == 7) {
      m_histo_list[key] = TH2F(key, std::forward<T>(args)...);
    } else
      static_assert(false, "histogram book error ... ");
  }

  template <typename... T> void Fill(const std::string &key, T... val) {
    if (m_histo_list.find(key) == m_histo_list.end()) {
      throw std::runtime_error(fmt::format("{} not booked ? ", key));
    }
    if constexpr (sizeof...(val) == 1) {
      std::get<TH1F>(m_histo_list[key]).Fill(val...);
    } else if constexpr (sizeof...(val) == 2) {
      std::get<TH2F>(m_histo_list[key]).Fill(val...);
    } else {
      static_assert(false, "histogram fill error... : ");
    }
  }

  void Write(const char *fname) {
    auto *f = TFile::Open(fname, "RECREATE");
    for (const auto &[k, v] : m_histo_list) {
      std::visit([](auto &el) { el.Write(); }, v);
    }
    f->Close();
  }
};

enum pdgId : int {
  Bplus = 521,
  Bminus = -521,
  Kplus = 321,
  Kminus = -321,
  X3872 = 9120443,
  Pi0 = 111,
  Pi_plus = 211,
  Pi_minus = -211,
  Photon = 22,
  Psi = 443,
  Psi2S = 100443,
  D0 = 421,
  D0_bar = -421,
  Dplus = 411,
  Dminus = -411
};

// pick the last particle from the record as in Section-4 here :
// https://pythia.org/download/pdf/worksheet8183.pdf

auto find_particle = [](const auto *particles, int id) {
  auto rit =
      std::find_if(particles->rbegin(), particles->rend(),
                   [&id](const auto &p) { return std::abs(p.id()) == id; });
  return std::tuple{rit, rit == particles->rend() ? false : true};
};

auto print_daughters = [](const auto &record, const auto &ptcl,
                          bool all = false) {
  const auto &daughters =
      all ? ptcl.daughterListRecursive() : ptcl.daughterList();
  std::vector<std::string> names(daughters.size());
  std::transform(daughters.begin(), daughters.end(), names.begin(),
                 [&record](int i) { return (*record)[i].name(); });
  fmt::print("{} -> {}\n", ptcl.name(), names);
};

int main() {

  HistogramRegistry hists;

  hists.Book("h_photons_all_E", "E_{#gamma}", 100, 0, 8);
  hists.Book("h_all_mult", " particle multiplicity", 100, 0, 1000);
  hists.Book("h_photon_mult", "photon multiplicity", 100, 0, 1000);
  hists.Book("h_charged_mult", "charged particle multiplicity", 50, 0, 1000);

  hists.Book("h_B_pt", "B meson", 100, 0, 100);
  hists.Book("h_B_eta", "B meson", 100, -3, 3);
  hists.Book("h_B_phi", "B meson", 100, -4, 4);
  hists.Book("h_B_E", "B meson", 100, 0, 150);

  hists.Book("h_K_pt", "Kaon", 100, 0, 100);
  hists.Book("h_K_eta", "Kaon", 100, -3, 3);
  hists.Book("h_K_phi", "Kaon", 100, -4, 4);
  hists.Book("h_K_E", "Kaon", 100, 0, 150);

  hists.Book("h_X3872_pt", "X(3872)", 100, 0, 100);
  hists.Book("h_X3872_eta", "X(3872)", 100, -3, 3);
  hists.Book("h_X3872_phi", "X(3872)", 100, -4, 4);
  hists.Book("h_X3872_E", "X(3872)", 100, 0, 150);

  Pythia pythia;
  // Add X(3872) meson to Pythia's database or something ...
  pythia.particleData.addParticle(9120443, "X_3872", "X_3872_bar", 3, 0, 0,
                                  3.87169, 0.00122, 0, 0, 0);

  pythia.readFile("x3872_settings.cmnd");

  pythia.init();

  size_t nGenerated{};
  constexpr size_t nTarget{10000};

  while (nGenerated++ < nTarget) {
    if (!pythia.next()) {
      continue;
    }

    const auto *record = pythia.event.particles();
    // shortly if a B meson is found ...
    if (const auto &[B, found] = find_particle(record, pdgId::Bplus); found) {
      auto x3872 = record->at(B->daughter1());
      auto kaon = record->at(B->daughter2());
      if (std::abs(x3872.id()) != pdgId::X3872) {
        std::swap(x3872, kaon);
      }

      hists.Fill("h_B_pt", B->pT());
      hists.Fill("h_B_eta", B->eta());
      hists.Fill("h_B_phi", B->phi());
      hists.Fill("h_B_E", B->e());

      hists.Fill("h_K_pt", kaon.pT());
      hists.Fill("h_K_eta", kaon.eta());
      hists.Fill("h_K_phi", kaon.phi());
      hists.Fill("h_K_E", kaon.e());

      hists.Fill("h_X3872_pt", x3872.pT());
      hists.Fill("h_X3872_eta", x3872.eta());
      hists.Fill("h_X3872_phi", x3872.phi());
      hists.Fill("h_X3872_E", x3872.e());


#ifdef DEBUG_GEN
      fmt::print("{} -> {} + {} \n ", B->name(), x3872.name(), kaon.name());
      print_daughters(record, kaon);
      print_daughters(record, x3872, true);
#endif
    }

    // filter all final particles, fill some histograms  ...:
    size_t nPhotons{};
    std::for_each(record->begin(), record->end(), [&](const auto &p) {
      if (!p.isFinal())
        return;
      if (p.id() == pdgId::Photon) {
        ++nPhotons;
        hists.Fill("h_photons_all_E", p.e());
      }
    });

    hists.Fill("h_all_mult", pythia.event.nFinal());
    hists.Fill("h_charged_mult", pythia.event.nFinal(true));
    hists.Fill("h_photon_mult", nPhotons);

  } // event loop

  hists.Write("x_3872.root");
  pythia.stat();
  return 0;
}
