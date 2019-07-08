/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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


  /**
   * @brief Entry kind for the tape value.
   *
   * Used for formatting.
   */
  enum class EntryType {
    Double,
    Int
  };

  /**
   * @brief Named entry in the tape values.
   *
   * The members are the name, the type and the location in the data arrays.
   */
  typedef std::tuple<std::string, EntryType, size_t> Entry;

  /**
   * @brief Data for one section in the tape values.
   */
  struct ValueSection {
    std::string name; /**< The name of the section */
    std::vector<Entry> data; /**< The data for the section */
  };

  /**
   * @brief Storage class for tape specific performance values.
   *
   * The class is used to gather information about the recorded tape and to output it.
   * Each tape can provide these values and every module of the tape will add its own values. The usual values are
   * the number of stored entries, the number of chunks and the allocated memory.
   *
   * In general class will also gather the total memory used and allocated.
   *
   * The function #formatDefault can be use to provide a pretty print of the values. #formatHeader and #formatRow output
   * the values in a csv table.
   */
  class TapeValues {
    private:
      std::vector<ValueSection> sections;  /**< All sections for this tape */

      std::vector<double> doubleData; /**< Stored double data */
      std::vector<size_t> intData; /**< Stored integer data */

      size_t usedMemoryIndex; /**< Index for the total used memory */
      size_t allocatedMemoryIndex; /**< Index for the total allocated memory */

    public:

      /**
       * @brief Initializes the total memory values.
       *
       * @param[in] tapeName  The name of the tape.
       */
      TapeValues(const std::string& tapeName) :
        sections(),
        doubleData(),
        intData(),
        usedMemoryIndex(0),
        allocatedMemoryIndex(1) {
        doubleData.push_back(0.0); // usedMemory
        doubleData.push_back(0.0); // allocated Memory

        addSection(tapeName);
        addDataInternal(std::make_tuple("Total memory used", EntryType::Double, usedMemoryIndex));
        addDataInternal(std::make_tuple("Total memory allocated", EntryType::Double, allocatedMemoryIndex));
      }

      /**
       * @brief Return the size of the used memory of the tape.
       * @return Memory in MB.
       */
      double getUsedMemorySize() {
        return doubleData[0];
      }

      /**
       * @brief Return the size of the allocated memory of the tape.
       * @return Memory in MB.
       */
      double getAllocatedMemorySize() {
        return doubleData[1];
      }

      /**
       * @brief Start a new section.
       *
       * All following addData calls will add there data to this section.
       *
       * @param[in] name  The name of the section.
       */
      void addSection(const std::string& name) {
        sections.resize(sections.size() + 1);

        sections.back().name = name;
      }

      /**
       * @brief Add a integer item to the currently active section.
       *
       * @param[in]  name  The name of the data item.
       * @param[in] value  The value of the data item.
       */
      void addData(const std::string& name, size_t value) {
        size_t pos = intData.size();
        intData.push_back(value);

        addDataInternal(std::make_tuple(name, EntryType::Int, pos));
      }

      /**
       * @brief Add a double item to the currently active section.
       *
       * The method will also update the global memory counts.
       *
       * @param[in]         name  The name of the data item.
       * @param[in]        value  The value of the data item.
       * @param[in]      usedMem  true if the value counts as used memory.
       * @param[in] allocatedMem  true if the value counts as allocated memory.
       */
      void addData(const std::string& name, double value, bool usedMem = false, bool allocatedMem = false) {
        size_t pos = doubleData.size();
        doubleData.push_back(value);

        addDataInternal(std::make_tuple(name, EntryType::Double, pos));

        if(usedMem) {
          doubleData[usedMemoryIndex] += value;
        }

        if(allocatedMem) {
          doubleData[allocatedMemoryIndex] += value;
        }

      }

      /**
       * @brief Add the default data of a data stream to the tape values.
       *
       * Adds total number of entries, number of chunks, used memory and allocated memory
       * as values.
       *
       * @param[in] stream  The data stream from which the values are taken.
       *
       * @tparm Stream  Either a ChunkVector or a SingleChunkVector
       */
      template<typename Stream>
      void addStreamData(const Stream& stream) {

        size_t numberOfChunks = stream.getNumChunks();
        size_t dataEntries    = stream.getDataSize();
        size_t entrySize      = Stream::ChunkType::EntrySize;

        double  memoryUsed  = (double)dataEntries*(double)entrySize* BYTE_TO_MB;
        double  memoryAlloc = (double)numberOfChunks*(double)stream.getChunkSize()*(double)entrySize* BYTE_TO_MB;

        addData("Total number", dataEntries);
        addData("Number of chunks", numberOfChunks);
        addData("Memory used", memoryUsed, true, false);
        addData("Memory allocated", memoryAlloc, false, true);
      }

      /**
       * @brief Output the default format.
       *
       * Each section header is enclosed in horizontal lines.
       * Each item is printed in by its name and value.
       * The end of a section is marked with a horizontal line.
       *
       *
       * \code{.txt}
       * -------------------------------------
       * CoDi Tape Statistics (ChunkTape)
       * -------------------------------------
       * Adjoint vector
       * -------------------------------------
       *   Number of Adjoints:      14517
       *   Memory allocated:         0.11 MB
       * -------------------------------------
       * \endcode
       *
       * @param[in,out] out  The stream for the output.
       *
       * @tparam Stream  Needs to implement stream operations.
       */
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

      /**
       * @brief Output a formatted header of the available data.
       *
       * The output format is a semicolon separated csv table.
       *
       * The header are generated like: &lt;section name&gt;-&lt;value name&gt;
       *
       * @param[in,out] out  The stream for the output.
       *
       * @tparam Stream  Needs to implement stream operations.
       */
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

      /**
       * @brief Output a formatted data row of the available data.
       *
       * The output format is a semicolon separated csv table.
       *
       * @param[in,out] out  The stream for the output.
       *
       * @tparam Stream  Needs to implement stream operations.
       */
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

      /**
       * @brief Helper function that combines the tape data with an MPI_Allreduce on MPI_COMM_WORLD.
       */
      void addData() {
#ifdef MPI_VERSION
        MPI_Allreduce(MPI_IN_PLACE, doubleData.data(), doubleData.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, intData.data(), intData.size(), MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
      }

    private:

      /**
       * @brief Format a value according to its type.
       *
       * @param[in,out]    out  The stream for the output.
       * @param[in]       data  The data to be formatted.
       * @param[in] outputType  If the type of the data should be written.
       * @param[in]  fieldSize  The minimum size of the formatted data.
       *
       * @tparam Stream  Needs to implement stream operations.
       */
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

      /**
       * @brief Get an estimate how the long the string for the formatted data will be.
       *
       * @param[in] data  The data for the estimation.
       * @return The formatted size of the value.
       */
      size_t dataSize(const Entry& data) const {
        std::stringstream ss;
        formatValue(ss, data, false, 0);

        return ss.str().size();
      }

      /**
       * @brief Add the data to the internal structure.
       *
       * Ensures that no data added to an empty section.
       *
       * @param[in] entry  The entry to be added.
       */
      void addDataInternal(Entry entry){
        if(sections.empty()) {
          addSection("General");
        }

        sections.back().data.push_back(entry);
      }
  };
}
