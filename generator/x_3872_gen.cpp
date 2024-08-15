// A Ptyhia8  test program to produce secondary X(3872)s from B decays
/*
  author : Serhat Istin
          istin@cern.ch
*/
// #define DEBUG_GEN

#include <algorithm>
#include <cmath>
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
  Muon = 13,
  Electron = 11,
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

auto is_muon = [](const auto &p) { return p.id() == pdgId::Muon; };

auto is_x3872 = [](const auto &p) { return std::abs(p.id()) == pdgId::X3872; };

auto is_kaon = [](const auto &p) { return std::abs(p.id()) == pdgId::Kplus; };

auto is_psi = [](const auto &p) { return std::abs(p.id()) == pdgId::Psi; };

auto is_psi2s = [](const auto &p) { return std::abs(p.id()) == pdgId::Psi2S; };

auto is_decayed = [](const auto &p) { return p.status() < 0; };

auto d0 = [](const auto &p) {
  return std::sqrt(p.xDec() * p.xDec() + p.yDec() * p.yDec());
};

auto dist = [](const auto &p) {
  return std::sqrt(p.xDec() * p.xDec() + p.yDec() * p.yDec() +
                   p.zDec() * p.zDec());
};

template <typename F, typename P = Pythia8::Particle>
concept ParticlePredicate = requires(F pred, P ptcl) {
  { pred(ptcl) } -> std::same_as<bool>;
};

template <typename Cont = std::vector<Pythia8::Particle>, ParticlePredicate P>
auto copy_daughters(const Cont &record, const auto &ptcl, P &&pred = all) {
  Cont daughters;
  //for (auto i : ptcl.daughterListRecursive()) 
  for (auto i : ptcl.daughterList()) {
    if (!pred(record[i]))
      continue;
    daughters.push_back(record[i]);
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

  //todo : separated booking code 
  hists.Book("h_photons_all_E", "E_{#gamma}", 100, 0, 8);
  hists.Book("h_all_mult", " particle multiplicity", 500, 0, 500);
  hists.Book("h_photon_mult", "photon multiplicity", 500, 0, 500);

  hists.Book("h_B_pt", "B meson", 100, 0, 100);
  hists.Book("h_B_eta", "B meson", 100, -3, 3);
  hists.Book("h_B_phi", "B meson", 100, -4, 4);
  hists.Book("h_B_E", "B meson", 100, 0, 150);
  hists.Book("h_B_N", "Bmeson", 10, 0, 10);
  hists.Book("h_B_zProd", "B meson", 100, -100, 100);
  hists.Book("h_B_zDec", "B meson", 100, -100, 100);
  hists.Book("h_B_d0", "B meson", 100, 0, 100);

  hists.Book("h_X3872_pt", "X(3872)", 100, 0, 100);
  hists.Book("h_X3872_eta", "X(3872)", 100, -3, 3);
  hists.Book("h_X3872_phi", "X(3872)", 100, -4, 4);
  hists.Book("h_X3872_E", "X(3872)", 100, 0, 150);
  hists.Book("h_X3872_N", "X(3872)", 10, 0, 10);
  hists.Book("h_X3872_zProd", "X(3872)", 100, -100, 100);
  hists.Book("h_X3872_zDec", "X(3872)", 100, -100, 100);
  hists.Book("h_X3872_d0", "X(3872)", 100, 0, 100);
  hists.Book("h_X3872_dist", "X(3872)", 100, 0, 100);

  hists.Book("h_x3872_photon_pt", "X3872 #rightarrow #gamma", 100, 0, 60);
  hists.Book("h_x3872_photon_eta", "X3872 #rightarrow #gamma", 100, -3, 3);
  hists.Book("h_x3872_photon_phi", "X3872 #rightarrow #gamma", 100, -4, 4);
  hists.Book("h_x3872_photon_e", "X3872 #rightarrow #gamma", 100, 0, 60);

  hists.Book("h_jpsi_pt", "X3872 #rightarrow J/#psi", 100, 0, 60);
  hists.Book("h_jpsi_eta", "X3872 #rightarrow J/#psi", 100, -3, 3);
  hists.Book("h_jpsi_phi", "X3872 #rightarrow J/#psi", 100, -4, 4);
  hists.Book("h_jpsi_e", "X3872 #rightarrow J/#psi", 100, 0, 60);

  hists.Book("h_jpsi_muon_pt", "J/#psi #rightarrow #mu", 100, 0, 60);
  hists.Book("h_jpsi_muon_eta", "J/#psi #rightarrow #mu", 100, -3, 3);
  hists.Book("h_jpsi_muon_phi", "J/#psi #rightarrow #mu", 100, -4, 4);
  hists.Book("h_jpsi_muon_e", "J/#psi #rightarrow #mu", 100, 0, 60);

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
    auto b_meson_list = decayed_particles | filter(is_b_meson);
    auto x3872_list = decayed_particles | filter(is_x3872);
    auto psi2s_list = decayed_particles | filter(is_psi2s);
    auto final_particles = record | filter(is_final);
    auto all_photons = final_particles | filter(is_photon);

    for (const auto &B : b_meson_list) {
      hists.Fill("h_B_pt", B.pT());
      hists.Fill("h_B_eta", B.eta());
      hists.Fill("h_B_phi", B.phi());
      hists.Fill("h_B_E", B.e());
      hists.Fill("h_B_zProd", B.zProd());
      hists.Fill("h_B_zDec", B.zDec());
      hists.Fill("h_B_d0", d0(B));
    }

    for (const auto &x3872 : x3872_list) {
      hists.Fill("h_X3872_pt", x3872.pT());
      hists.Fill("h_X3872_eta", x3872.eta());
      hists.Fill("h_X3872_phi", x3872.phi());
      hists.Fill("h_X3872_E", x3872.e());
      hists.Fill("h_X3872_zProd", x3872.zProd());
      hists.Fill("h_X3872_zDec", x3872.zDec());
      hists.Fill("h_X3872_d0", d0(x3872));
      hists.Fill("h_X3872_dist", dist(x3872));
      for(const auto& ph : copy_daughters(record,x3872,is_photon)){
        hists.Fill("h_x3872_photon_pt", ph.pT());
        hists.Fill("h_x3872_photon_eta", ph.eta());
        hists.Fill("h_x3872_photon_phi", ph.phi());
        hists.Fill("h_x3872_photon_e", ph.e());
      }
    }

    for (const auto &jpsi : decayed_particles | filter(is_psi) ) {
      hists.Fill("h_jpsi_pt", jpsi.pT());
      hists.Fill("h_jpsi_eta", jpsi.eta());
      hists.Fill("h_jpsi_phi", jpsi.phi());
      hists.Fill("h_jpsi_e", jpsi.e());
      for (const auto &mu : copy_daughters(record, jpsi, is_muon)) {
        hists.Fill("h_jpsi_muon_pt", mu.pT());
        hists.Fill("h_jpsi_muon_eta", mu.eta());
        hists.Fill("h_jpsi_muon_phi", mu.phi());
        hists.Fill("h_jpsi_muon_e", mu.e());
      }

    }

    hists.Fill("h_B_N", std::ranges::distance(b_meson_list));
    hists.Fill("h_X3872_N", std::ranges::distance(x3872_list));
    hists.Fill("h_all_mult", std::ranges::distance(final_particles));
    hists.Fill("h_photon_mult", std::ranges::distance(all_photons));

  } // event loop

  hists.Write("x_3872.root");
  pythia.stat();
  return 0;
}
