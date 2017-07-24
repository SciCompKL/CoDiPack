/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include <iomanip>
#include <sstream>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>

#include "../exceptions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

    enum class EntryType {
      Double,
      Int
    };

    typedef std::tuple<std::string, EntryType, size_t> Entry;

    struct ValueSection {
      std::string name;
      std::vector<Entry> data;
    };

    class TapeValues {
      private:
        std::vector<ValueSection> sections;

        std::vector<double> doubleData;
        std::vector<size_t> intData;

      public:
        void addSection(const std::string& name) {
          sections.resize(sections.size() + 1);

          sections.back().name = name;
        }

        void addData(const std::string& name, size_t value) {
          size_t pos = intData.size();
          intData.push_back(value);

          addDataInternal(std::make_tuple(name, EntryType::Int, pos));
        }

        void addData(const std::string& name, double value) {
          size_t pos = doubleData.size();
          doubleData.push_back(value);

          addDataInternal(std::make_tuple(name, EntryType::Double, pos));
        }

        template<typename Stream = std::ostream>
        void formatDefault(Stream& out = std::cout) const {

          const std::string hLine = "-------------------------------------\n";

          size_t maxNameSize = 0;
          size_t maxDataSize = 10;
          for(const ValueSection& section : sections) {
            for(const Entry& data : section.data) {
              maxNameSize = std::max(maxNameSize, std::get<0>(data).size());
              maxDataSize = std::max(maxDataSize, dataSize(data));
            }
          }

          out << hLine;
          for(const ValueSection& section : sections) {
            out << std::left << section.name << "\n";
            out << hLine;
            for(const Entry& data : section.data) {
              out << "  " << std::left << std::setw(maxNameSize) << std::get<0>(data) << " : ";
              formatValue(out, data, true, maxDataSize);
              out << "\n";
            }
            if(!section.data.empty()) {
              out << hLine;
            }
          }
        }

        template<typename Stream = std::ostream>
        void formatHeader(Stream& out = std::cout) const {

          bool first = true;
          for(const ValueSection& section : sections) {
            for(const Entry& data : section.data) {

              if(first) {
                first = false;
              } else {
                out << "; ";
              }
              out << section.name << "-" << std::get<0>(data);
            }
          }

          out << "\n";
        }

        template<typename Stream = std::ostream>
        void formatRow(Stream& out = std::cout) const {

          size_t maxDataSize = 10;
          for(const ValueSection& section : sections) {
            for(const Entry& data : section.data) {
              maxDataSize = std::max(maxDataSize, dataSize(data));
            }
          }

          bool first = true;
          for(const ValueSection& section : sections) {
            for(const Entry& data : section.data) {

              if(first) {
                first = false;
              } else {
                out << "; ";
              }
              formatValue(out, data, false, maxDataSize);
            }
          }

          out << "\n";
        }

      private:

        template<typename Stream>
        void formatValue(Stream& out, const Entry& data, bool outputType, int fieldSize) const {
          switch (std::get<1>(data)) {
          case EntryType::Int:
            out << std::right << std::setw(fieldSize) << intData[std::get<2>(data)];
            break;
          case EntryType::Double:
            out << std::right << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(fieldSize) << doubleData[std::get<2>(data)];
            if(outputType) {
              out << " MB";
            }
            break;
          default:
            CODI_EXCEPTION("Unimplemented switch case.");
            break;
          }
        }

        size_t dataSize(const Entry& data) const {
          std::stringstream ss;
          formatValue(ss, data, false, 0);

          return ss.str().size();
        }

        void addDataInternal(Entry entry){
          if(sections.empty()) {
            addSection("General");
          }

          sections.back().data.push_back(entry);
        }
    };
}
