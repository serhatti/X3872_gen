// A Ptyhia8  test program to produce secondary X(3872)s from B decays
/*
  author : Serhat Istin
          istin@cern.ch
*/

#include <algorithm>
#include <fmt/ranges.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <variant>

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace Pythia8;

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
  std::vector<std::string> names;
  const auto &daughters =
      all ? ptcl.daughterListRecursive() : ptcl.daughterList();
  for (auto i : daughters) {
    names.push_back(std::move(record->at(i).name()));
  }
  fmt::print("{} -> {} \n", ptcl.name(), names);
};

int main() {
  HistogramRegistry hists;

  hists.Book("h_e_photons_all", "E_{#gamma}", 100, 0, 8);
  hists.Book("h_all_mult", " particle multiplicity", 50, 0, 1000);
  hists.Book("h_photon_mult", "photon multiplicity", 50, 0, 1000);
  hists.Book("h_charged_mult", "charged particle multiplicity", 50, 0, 1000);

  Pythia pythia;

  // p-p collisions at sqrt(s) in [GeV]
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");

  // setup proc. for B meson production
  pythia.readString("HardQCD:gg2bbbar = on");
  pythia.readString("HardQCD:qqbar2bbbar = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = on");
  pythia.readString("ProcessLevel:resonanceDecays = on");

  // Add X(3872) meson to Pythia's database or something ...
  pythia.particleData.addParticle(9120443, "X_3872", "X_3872_bar", 3, 0, 0,
                                  3.87169, 0.00122, 0, 0, 0);
  //..............X(3872) decay modes .....................
  // X(3872) --> gamma + gamma     channel
  pythia.readString("9120443:addChannel = 1 1 0 22 22");
  // add X(3872) --> gamma + Psi(2S)      channel
  pythia.readString("9120443:addChannel = 1 1 0 22 100443");
  //  X(3872) --> gamma + J/Psi(1S)     channel
  pythia.readString("9120443:addChannel = 1 1 0 22 443");
  // X(3872) --> D0 + D0bar + pi0      channel
  pythia.readString("9120443:addChannel = 1 1 0 421 -421 111");
  // X(3872) -->  D0 + D0_bar     channel
  pythia.readString("9120443:addChannel = 1 1 0 421 -421");
  // X(3872) -->  D+ + D-      channel
  pythia.readString("9120443:addChannel = 1 1 0 411 -411");

  pythia.readString("9120443:mayDecay = on");
  /* first turn of all decay modes for X(3872)
      then turn on as you like
  */
  pythia.readString("9120443:onMode = off");
  pythia.readString("9120443:onIfAll = 22 22");

  /* .............Settings for B mesons ................... */
  // turn-off b quark decay
  pythia.readString("4:mayDecay = off");
  // Force B decays to : B+ --> X(3872) + K+  + H.C
  pythia.readString("521:oneChannel = 1 1 0 9120443 321");

  //  See below with caution ! :
  /* https://pythia.org/latest-manual/ParticleDecays.html :
   * ii) The main switch for allowing this particle kind to decay must be on;
   * tested by the mayDecay() method of Event (and ParticleData). By default
   * this is defined as true for all particles with tau0 below 1000 mm, and
   * false for ones above, see the Particle Data Scheme. This means that mu^+-,
   * pi^+-, K^+-, K^0_L and n/nbar always remain stable unless decays are
   * explicity switched on, e.g. 211:mayDecay = true.
   */

  pythia.init();

  size_t nGenerated{};
  constexpr size_t nTarget{1000};

  while (nGenerated++ < nTarget) {
    if (!pythia.next()) {
      continue;
    }

    const auto *record = pythia.event.particles();
    // shortly if a B meson is found ...
    if (auto [B, found] = find_particle(record, pdgId::Bplus); found) {
      auto x3872 = record->at(B->daughter1());
      auto kaon = record->at(B->daughter2());

      if (std::abs(x3872.id()) != pdgId::X3872) {
        std::swap(x3872, kaon);
      }

      fmt::print("{} -> {} + {} \n ", B->name(), x3872.name(), kaon.name());
      print_daughters(record, kaon);
      print_daughters(record, x3872);
    }

    // all final particles fill some histograms etc ...:
    size_t nPhotons{};
    std::for_each(record->begin(), record->end(), [&](const auto &p) {
      if (!p.isFinal())
        return;
      if (p.id() == pdgId::Photon) {
        ++nPhotons;
        hists.Fill("h_e_photons_all", p.e());
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
