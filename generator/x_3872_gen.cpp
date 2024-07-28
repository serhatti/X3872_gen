// A Ptyhia8  test program to produce secondary X(3872)s from B decays
/*
  author : Serhat Istin
          istin@cern.ch
*/
// #define DEBUG_GEN

#include <algorithm>
#include <fmt/ranges.h>
#include <format>
#include <ranges>
#include <tuple>

#include "Pythia8/Pythia.h"
#include "include/HistogramRegistry.h"

enum pdgId : int {
  Bplus = 521,
  Kplus = 321,
  X3872 = 9120443,
  Pi0 = 111,
  Pi_plus = 211,
  Photon = 22,
  Psi = 443,
  Psi2S = 100443,
  D0 = 421,
  Dplus = 411,
};

template <typename F, typename P = Pythia8::Particle>
concept ParticlePredicate = requires(F& predicate, P ptcl) {
  { predicate(ptcl) } -> std::same_as<bool>;
};

template <ParticlePredicate P>
auto filter_particles(const auto *record, P pred)

{
  std::vector<Pythia8::Particle> particles(record->size());
  auto it =
      std::copy_if(record->cbegin(), record->cend(), particles.begin(), pred);
  particles.resize(std::distance(particles.begin(), it));
  return particles;
};

template <ParticlePredicate P>
auto filter_indices(const auto *record, P pred)

{
  std::vector<size_t> indices(4);
  for (const auto &particle : *record) {
    if (pred(particle)) {
      indices.push_back(particle.index());
    }
  }
  return indices;
};

auto print_daughters = [](const auto *record, const auto &ptcl,
                          bool all = false) {
  const auto &daughters =
      all ? ptcl.daughterListRecursive() : ptcl.daughterList();
  std::vector<decltype(ptcl.name())> names(daughters.size());
  std::transform(daughters.begin(), daughters.end(), names.begin(),
                 [&record](int i) { return (*record)[i].name(); });
  fmt::println("{} -> {}", ptcl.name(), names);
};

auto is_b_meson = [](const auto &p) {
  return std::abs(p.id()) == pdgId::Bplus and
         (p.status() == -83 or p.status() == -84);
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
  hists.Book("h_B_N", "Bmeson", 10, 0, 10);

  hists.Book("h_K_pt", "Kaon", 100, 0, 100);
  hists.Book("h_K_eta", "Kaon", 100, -3, 3);
  hists.Book("h_K_phi", "Kaon", 100, -4, 4);
  hists.Book("h_K_E", "Kaon", 100, 0, 150);

  hists.Book("h_X3872_pt", "X(3872)", 100, 0, 100);
  hists.Book("h_X3872_eta", "X(3872)", 100, -3, 3);
  hists.Book("h_X3872_phi", "X(3872)", 100, -4, 4);
  hists.Book("h_X3872_E", "X(3872)", 100, 0, 150);

  hists.Book("h_X3872_d1_pt", "X(3872) daughter1", 100, 0, 100);
  hists.Book("h_X3872_d1_eta", "X(3872) daughter1", 100, -3, 3);
  hists.Book("h_X3872_d1_phi", "X(3872) daughter1", 100, -4, 4);
  hists.Book("h_X3872_d1_E", "X(3872) daughter1", 100, 0, 150);

  hists.Book("h_X3872_d2_pt", "X(3872) daughter2", 100, 0, 100);
  hists.Book("h_X3872_d2_eta", "X(3872) daughter2", 100, -3, 3);
  hists.Book("h_X3872_d2_phi", "X(3872) daughter2", 100, -4, 4);
  hists.Book("h_X3872_d2_E", "X(3872) daughter2", 100, 0, 150);

  Pythia8::Pythia pythia;
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
    const auto &b_mesons = filter_particles(record, is_b_meson);

    for (const auto &B : b_mesons) {
      auto kaon = record->at(B.daughter1());
      auto x3872 = record->at(B.daughter2());
      if (std::abs(x3872.id()) != pdgId::X3872) {
        std::swap(kaon, x3872);
      }
      hists.Fill("h_B_pt", B.pT());
      hists.Fill("h_B_eta", B.eta());
      hists.Fill("h_B_phi", B.phi());
      hists.Fill("h_B_E", B.e());
      hists.Fill("h_B_N", b_mesons.size());

      hists.Fill("h_X3872_pt", x3872.pT());
      hists.Fill("h_X3872_eta", x3872.eta());
      hists.Fill("h_X3872_phi", x3872.phi());
      hists.Fill("h_X3872_E", x3872.e());

      hists.Fill("h_K_pt", kaon.pT());
      hists.Fill("h_K_eta", kaon.eta());
      hists.Fill("h_K_phi", kaon.phi());
      hists.Fill("h_K_E", kaon.e());

      const auto &x3872_d1 = record->at(x3872.daughter1());
      const auto &x3872_d2 = record->at(x3872.daughter2());
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
