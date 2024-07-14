// A Ptyhia8  test program to produce secondary X(3872)s from B decays
/*
  author : Serhat Istin
          istin@cern.ch
*/

#include <algorithm>
#include <string>
#include <unordered_map>
#include <variant>
#include <format>


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
      m_histo_list[key] = TH2F(key, args...);
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

enum pdgIds : int {
  Bplus = 521,
  Bminus = -521,
  Kplus = 321,
  Kminus = -321,
  X3872 = 9120443,
  Pi0 = 111,
  Pi_plus = 211,
  Pi_minus = -211,
  Photon = 22
};

// pick the last particle from the record as in Section-4 here :
// https://pythia.org/download/pdf/worksheet8183.pdf
auto find_last_particle = [](const auto *particles, int id) {
  return std::find_if(particles->rbegin(), particles->rend(),
                      [&id](const auto &p) { return p.id() == id; });
};

auto find_last_particle_or_anti = [](const auto *particles, int id) {
  return std::find_if(particles->rbegin(), particles->rend(),
                      [&id](const auto &p) { return std::abs(p.id()) == id; });
};

auto print_daughters = [](const auto &record, const auto &ptcl) {
  std::cout << ptcl.name() << " ->";
  for (auto i : ptcl.daughterListRecursive()) {
    std::cout << " " << record->at(i).name() << " + ";
  }
  std::cout << " \n";
};

int main() {
  size_t nGenerated{};
  constexpr size_t nTarget{1000};
  HistogramRegistry hists;

  hists.Book("h_e_photons_all", "E_{#gamma}", 100, 0, 8);
  hists.Book("h_all_mult", " particle multiplicity", 50, 0, 1000);
  hists.Book("h_photon_mult", "photon multiplicity", 50, 0, 1000);
  hists.Book("h_charged_mult", "charged particle multiplicity", 50, 0, 1000);

  Pythia pythia;

  // Add X(3872) meson to Pythia's database or something ...
  pythia.particleData.addParticle(9120443, "X_3872", "X_3872", 3, 0, 0, 3.87169,
                                  0.00122);
  pythia.readString("9120443:mayDecay = on");

  // p-p collisions at sqrt(s) in [GeV]
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:MPI = off");

  // setup proc. for B meson production
  pythia.readString("HardQCD:gg2bbbar = on");
  pythia.readString("HardQCD:qqbar2bbbar = on");

  // turn-off b quark and charged pion decay
  pythia.readString("ProcessLevel:resonanceDecays = on");
  pythia.readString("4:mayDecay = off");
  pythia.readString("211:mayDecay = off");
  pythia.readString("443:mayDecay = off");

  //  See below with caution ! :

  /* https://pythia.org/latest-manual/ParticleDecays.html :
   * ii) The main switch for allowing this particle kind to decay must be on;
   * tested by the mayDecay() method of Event (and ParticleData). By default
   * this is defined as true for all particles with tau0 below 1000 mm, and
   * false for ones above, see the Particle Data Scheme. This means that mu^+-,
   * pi^+-, K^+-, K^0_L and n/nbar always remain stable unless decays are
   * explicity switched on, e.g. 211:mayDecay = true.
   */

  // only : B+ --> X(3872) + K+  + H.C
  pythia.readString("521:oneChannel = 1 1 0 9120443 321");
  //  X(3872) --> J/Psi(1S) + gamma
  pythia.readString("9120443:oneChannel = 1 1 0 443 22");
  // X(3872) --> Psi(2S) + gamma
  pythia.readString("9120443:addChannel = 1 1 0 100443 22");

  pythia.init();

  while (nGenerated++ < nTarget) {
    if (!pythia.next()) {
      continue;
    }

    const auto *particle_record = pythia.event.particles();
    // shortly if a B meson is found ...
    if (auto B_meson =
            find_last_particle_or_anti(particle_record, pdgIds::Bplus);
        B_meson != particle_record->rend()) {
      auto x3872 = particle_record->at(B_meson->daughter1());
      auto kaon = particle_record->at(B_meson->daughter2());

      if (std::abs(x3872.id()) != pdgIds::X3872) {
        std::swap(x3872, kaon);
      }

      std::cout << B_meson->name() << " --> " << x3872.name() << " + "
                << kaon.name() << "\n";

      print_daughters(particle_record, kaon);
      print_daughters(particle_record, x3872);
    }

    // all final particles fill some histograms etc ...:
    size_t nPhotons{};
    std::for_each(particle_record->begin(), particle_record->end(),
                  [&](const auto &p) {
                    if (!p.isFinal())
                      return;
                    if (p.id() == pdgIds::Photon) {
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
