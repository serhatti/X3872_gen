/*
  author : Serhat Istin
          istin@cern.ch
*/

#ifndef HISTOGRAMREGISTRY_H
#define HISTOGRAMREGISTRY_H

#include <string>
#include <type_traits>
#include <unordered_map>
#include <variant>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class HistogramRegistry {
  using Histogram = std::variant<TH1F, TH2F>;
  std::unordered_map<std::string, Histogram> m_histo_list;
  std::vector<std::string> m_keys;

public:
  template <typename K, typename... T> void Book(const K &key, T &&...args) {
    if constexpr (std::is_constructible_v<TH1F, K, T...>) {
      m_histo_list[key] = TH1F(key, std::forward<T>(args)...);
    } else if constexpr (std::is_constructible_v<TH2F, K, T...>) {
      m_histo_list[key] = TH2F(key, std::forward<T>(args)...);
    } else
      static_assert(false, "histogram book error ... ");
    m_keys.push_back(key);
  }

  template <typename... T> void Fill(const std::string &key, T... val) {
    if (m_histo_list.find(key) == m_histo_list.end()) {
      throw std::runtime_error(std::format("{} not boked ?", key));
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
    std::sort(m_keys.begin(), m_keys.end());
    auto *f = TFile::Open(fname, "RECREATE");
    for (const auto &k : m_keys) {
      std::visit([](auto &el) { el.Write(); }, m_histo_list[k]);
    }
    f->Close();
  }
};

#endif