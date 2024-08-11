// A Ptyhia8  test program to produce secondary X(3872)s from B decays
/*
  author : Serhat Istin
          istin@cern.ch
*/
// #define DEBUG_GEN

#include <algorithm>
#include <fmt/ranges.h>
#include <ranges>

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

auto all = [](const auto &p) { return true; };

auto is_final = [](const auto &p) { return p.isFinal(); };

// pick a b meson produced via hadronization , B is also decayed and not present
// anymore
auto is_b_meson = [](const auto &p) {
  return std::abs(p.id()) == pdgId::Bplus;
};

auto is_photon = [](const auto &p) { return p.id() == pdgId::Photon; };

auto is_x3872 = [](const auto &p) { return std::abs(p.id()) == pdgId::X3872; };

auto is_kaon = [](const auto &p) { return std::abs(p.id()) == pdgId::Kplus; };

auto is_decayed = [](const auto &p) { return p.status() < 0; };

template <typename F, typename P = Pythia8::Particle>
concept ParticlePredicate = requires(F pred, P ptcl) {
  { pred(ptcl) } -> std::same_as<bool>;
};

template <typename Cont = std::vector<Pythia8::Particle>, ParticlePredicate P>
auto copy_daughters(const Cont &record, const auto &ptcl, P &&pred = all) {
  Cont daughters;
  for (auto i : ptcl.daughterListRecursive()) {
    if (!pred(record[i]))
      continue;
    daughters.push_back(ptcl);
  }
  return daughters;
}

auto print_daughters = [](const auto &record, const auto &ptcl,
                          bool all = false) {
  const auto &daughters =
      all ? ptcl.daughterListRecursive() : ptcl.daughterList();
  std::vector<decltype(ptcl.name())> names(daughters.size());
  std::transform(daughters.begin(), daughters.end(), names.begin(),
                 [&record](int i) { return record[i].name(); });
  fmt::println("{} -> {}", ptcl.name(), names);
};

int main() {

  HistogramRegistry hists;

  hists.Book("h_photons_all_E", "E_{#gamma}", 100, 0, 8);
  hists.Book("h_all_mult", " particle multiplicity", 500, 0, 500);
  hists.Book("h_photon_mult", "photon multiplicity", 500, 0, 500);

  hists.Book("h_B_pt", "B meson", 100, 0, 100);
  hists.Book("h_B_eta", "B meson", 100, -3, 3);
  hists.Book("h_B_phi", "B meson", 100, -4, 4);
  hists.Book("h_B_E", "B meson", 100, 0, 150);
  hists.Book("h_B_N", "Bmeson", 10, 0, 10);

  hists.Book("h_X3872_pt", "X(3872)", 100, 0, 100);
  hists.Book("h_X3872_eta", "X(3872)", 100, -3, 3);
  hists.Book("h_X3872_phi", "X(3872)", 100, -4, 4);
  hists.Book("h_X3872_E", "X(3872)", 100, 0, 150);
  hists.Book("h_X3872_N", "X(3872)", 10, 0, 10);

  hists.Book("h_X3872_d1_pt", "X(3872) daughter1", 100, 0, 60);
  hists.Book("h_X3872_d1_eta", "X(3872) daughter1", 100, -3, 3);
  hists.Book("h_X3872_d1_phi", "X(3872) daughter1", 100, -4, 4);
  hists.Book("h_X3872_d1_E", "X(3872) daughter1", 100, 0, 60);

  hists.Book("h_X3872_d2_pt", "X(3872) daughter2", 100, 0, 60);
  hists.Book("h_X3872_d2_eta", "X(3872) daughter2", 100, -3, 3);
  hists.Book("h_X3872_d2_phi", "X(3872) daughter2", 100, -4, 4);
  hists.Book("h_X3872_d2_E", "X(3872) daughter2", 100, 0, 60);

  Pythia8::Pythia pythia;
  // Add X(3872) meson to Pythia's database or something ...
  pythia.particleData.addParticle(pdgId::X3872, "X_3872", "X_3872_bar", 3, 0, 0,
                                  3.87169, 0.00122, 0, 0, 0);

  pythia.readFile("x3872_settings.cmnd");

  pythia.init();

  size_t nGenerated{};
  constexpr size_t nTarget{10000};

  while (nGenerated++ < nTarget) {
    if (!pythia.next()) {
      continue;
    }

    const auto &record = *(pythia.event.particles());
    using std::views::filter;
    auto decayed_particles = record | filter(is_decayed);
    auto b_mesons = decayed_particles | filter(is_b_meson);
    auto x3872s = decayed_particles | filter(is_x3872);
    auto final_particles = record | filter(is_final);
    auto all_photons = final_particles | filter(is_photon);

    for (const auto &B : b_mesons) {
      hists.Fill("h_B_pt", B.pT());
      hists.Fill("h_B_eta", B.eta());
      hists.Fill("h_B_phi", B.phi());
      hists.Fill("h_B_E", B.e());
    }

    for (const auto &x3872 : x3872s) {
      hists.Fill("h_X3872_pt", x3872.pT());
      hists.Fill("h_X3872_eta", x3872.eta());
      hists.Fill("h_X3872_phi", x3872.phi());
      hists.Fill("h_X3872_E", x3872.e());
      const auto &d1 = record[x3872.daughter1()];
      const auto &d2 = record[x3872.daughter2()];
      hists.Fill("h_X3872_d1_pt", d1.pT());
      hists.Fill("h_X3872_d1_eta", d1.eta());
      hists.Fill("h_X3872_d1_phi", d1.phi());
      hists.Fill("h_X3872_d1_E", d1.e());
      hists.Fill("h_X3872_d2_pt", d2.pT());
      hists.Fill("h_X3872_d2_eta", d2.eta());
      hists.Fill("h_X3872_d2_phi", d2.phi());
      hists.Fill("h_X3872_d2_E", d2.e());
#ifdef DEBUG_GEN
      fmt::println("d1={} d2={}", d1.name(), d2.name());
#endif
    }

    hists.Fill("h_B_N", std::ranges::distance(b_mesons));
    hists.Fill("h_X3872_N", std::ranges::distance(x3872s));
    hists.Fill("h_all_mult", std::ranges::distance(final_particles));
    hists.Fill("h_photon_mult", std::ranges::distance(all_photons));

  } // event loop

  hists.Write("x_3872.root");
  pythia.stat();
  return 0;
}
