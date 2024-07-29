/*
  author : Serhat Istin
          istin@cern.ch
*/


#ifndef HISTOGRAMREGISTRY_H
#define HISTOGRAMREGISTRY_H


#include <unordered_map>
#include <variant>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"


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
    if(m_histo_list.find(key) == m_histo_list.end()){
        throw std::runtime_error(std::format("{} not boked ?",key));
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


#endif