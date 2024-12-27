/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../../config.h"
#include "../../misc/exceptions.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Tape information that can be printed in a pretty print format or a table format.
   *
   * This structure is generated for tapes with the ReverseTapeInterface::getTapeValues() function. The tape provides
   * the information for all internal data structures and the user can then output this information for further
   * analysis. To make the output as useful as possible, tapes should provide information about all DataInterface
   * members, arrays, and IndexManagerInterface members.
   *
   * - Functions for adding data:
   *   - addDoubleEntry(): Add a double entry. If this a memory entry, it can be added automatically to the global
   *                       counters. Memory is computed in MB.
   *   - addLongEntry(): Add a long entry.
   *   - addUnsignedLongEntry(): Add unsigned long entry.
   *   - addSection(): Add a new section under which all following entries are added.
   *
   * - Format data:
   *   - formatDefault(): Default human readable format. One row per entry.
   *   - formatHeader(): Output the header for the table output.
   *   - formatRow(): Output the data in this object in one row. One column per entry.
   *
   * - Misc:
   *   - combineData(TapeVales const& other): Perform element-wise addition with other tape values.
   *   - combineData(): Deprecated. Kept for backwards compatibility. Perform an MPI_Allreduce on MPI_COMM_WORLD.
   *   - combineDataMPI(): Perform an MPI_Allreduce on a given communicator.
   *   - getAllocatedMemorySize(): Get the allocated memory size.
   *   - getUsedMemorySize(): Get the used memory size.
   */
  struct TapeValues {
    public:
      enum class LocalReductionOperation {
        Sum,
        Max
      };

    private:
      enum class EntryType {
        Double,
        Long,
        UnsignedLong
      };

      struct Entry {
        public:
          std::string name;
          EntryType type;
          LocalReductionOperation operation;
          size_t pos;

          Entry() : name(), type(), operation(), pos() {}

          Entry(std::string const& name, EntryType const& type, LocalReductionOperation const& operation,
                size_t const& pos)
              : name(name), type(type), operation(operation), pos(pos) {}
      };

      struct Section {
        public:
          std::string name;
          std::vector<Entry> data;

          Section() : name(), data() {}

          Section(std::string const& name) : name(name), data() {}
      };

      std::vector<Section> sections;

      std::vector<double> doubleData;
      std::vector<long> longData;
      std::vector<unsigned long> unsignedLongData;

      size_t usedMemoryIndex;
      size_t allocatedMemoryIndex;

    public:

      /// Constructor
      TapeValues(std::string const& tapeName)
          : sections(), doubleData(), longData(), unsignedLongData(), usedMemoryIndex(0), allocatedMemoryIndex(1) {
        addSection(tapeName);
        addEntryInternal("Total memory used", EntryType::Double, LocalReductionOperation::Sum, doubleData, 0.0);
        addEntryInternal("Total memory allocated", EntryType::Double, LocalReductionOperation::Sum, doubleData, 0.0);
      }

      /*******************************************************************************/
      /// @name Add data
      /// @{

      /// Add double entry. If it is a memory entry, it should be in bytes.
      void addDoubleEntry(std::string const& name, double const& value,
                          LocalReductionOperation operation = LocalReductionOperation::Sum, bool usedMem = false,
                          bool allocatedMem = false) {
        addEntryInternal(name, EntryType::Double, operation, doubleData, value);

        if (usedMem) {
          doubleData[usedMemoryIndex] += value;
        }

        if (allocatedMem) {
          doubleData[allocatedMemoryIndex] += value;
        }
      }

      /// Add long entry.
      void addLongEntry(std::string const& name, long const& value,
                        LocalReductionOperation operation = LocalReductionOperation::Sum) {
        addEntryInternal(name, EntryType::Long, operation, longData, value);
      }

      /// Add section. All further entries are added under this section.
      void addSection(std::string const& name) {
        sections.push_back(Section(name));
      }

      /// Add unsigned long entry.
      void addUnsignedLongEntry(std::string const& name, unsigned long const& value,
                                LocalReductionOperation operation = LocalReductionOperation::Sum) {
        addEntryInternal(name, EntryType::UnsignedLong, operation, unsignedLongData, value);
      }

      /// @}
      /*******************************************************************************/
      /// @name Format data
      /// @{

      /// Output in a human readable format. One row per entry.
      template<typename Stream = std::ostream>
      void formatDefault(Stream& out = std::cout) const {
        std::string const hLine = "-------------------------------------\n";

        size_t maxNameSize = getMaximumNameLength();
        size_t maxValueSize = std::max((size_t)10, getMaximumValueLength());

        out << hLine;
        for (Section const& section : sections) {
          out << std::left << section.name << "\n";
          out << hLine;

          for (Entry const& entry : section.data) {
            out << "  " << std::left << std::setw(maxNameSize) << entry.name << " : "
                << formatEntry(entry, maxValueSize) << "\n";
          }

          if (!section.data.empty()) {
            out << hLine;
          }
        }
      }

      /// Output the header for a table output.
      template<typename Stream = std::ostream>
      void formatHeader(Stream& out = std::cout) const {
        bool first = true;
        for (Section const& section : sections) {
          for (Entry const& entry : section.data) {
            if (first) {
              first = false;
            } else {
              out << "; ";
            }
            out << section.name << "-" << entry.name;
          }
        }

        out << "\n";
      }

      /// Output this data in one row. One entry per column.
      template<typename Stream = std::ostream>
      void formatRow(Stream& out = std::cout) const {
        size_t maxValueSize = std::max((size_t)10, getMaximumValueLength());

        bool first = true;
        for (Section const& section : sections) {
          for (Entry const& entry : section.data) {
            if (first) {
              first = false;
            } else {
              out << "; ";
            }
            out << formatEntry(entry, maxValueSize);
          }
        }

        out << "\n";
      }

      /// @}
      /*******************************************************************************/
      /// @name Misc.
      /// @{

      /// Perform entry-wise additions.
      void combineData(TapeValues const& other) {
        // Size check for the number of sections.
        codiAssert(this->sections.size() == other.sections.size());

        for (size_t section = 0; section < this->sections.size(); ++section) {
          auto& thisSection = this->sections[section];
          auto const& otherSection = other.sections[section];

          // Basic check to ensure that we combine identically structured tape values.
          codiAssert(thisSection.name == otherSection.name);

          // Size check for the number of entries.
          codiAssert(thisSection.data.size() == otherSection.data.size());

          for (size_t entry = 0; entry < thisSection.data.size(); ++entry) {
            auto& thisEntry = thisSection.data[entry];
            auto const& otherEntry = otherSection.data[entry];

            // Basic checks to ensure that we combine identically structured tape values.
            codiAssert(thisEntry.name == otherEntry.name);
            codiAssert(thisEntry.type == otherEntry.type);
            codiAssert(thisEntry.operation == otherEntry.operation);

            switch (thisEntry.type) {
              case EntryType::Double: {
                performLocalReduction(this->doubleData[thisEntry.pos], other.doubleData[otherEntry.pos],
                                      thisEntry.operation);
                break;
              }
              case EntryType::Long: {
                performLocalReduction(this->longData[thisEntry.pos], other.longData[otherEntry.pos],
                                      thisEntry.operation);
                break;
              }
              case EntryType::UnsignedLong: {
                performLocalReduction(this->unsignedLongData[thisEntry.pos], other.unsignedLongData[otherEntry.pos],
                                      thisEntry.operation);
                break;
              }
            }
          }
        }
      }

      /** @brief Perform an MPI_Allreduce with MPI_COMM_WORLD.
       *
       * This method is deprecated and only kept for backwards compatibility. combineDataMPI should be used instead.
       */
      void combineData() {
#ifdef MPI_VERSION
        combineDataMPI(MPI_COMM_WORLD);
#endif
      }

      /// Perform an MPI_Allreduce with the given communicator.
#ifdef MPI_VERSION
      void combineDataMPI(MPI_Comm communicator) {
        MPI_Allreduce(MPI_IN_PLACE, doubleData.data(), doubleData.size(), MPI_DOUBLE, MPI_SUM, communicator);
        MPI_Allreduce(MPI_IN_PLACE, longData.data(), longData.size(), MPI_LONG, MPI_SUM, communicator);
        MPI_Allreduce(MPI_IN_PLACE, unsignedLongData.data(), unsignedLongData.size(), MPI_UNSIGNED_LONG, MPI_SUM,
                      communicator);
      }
#else
      template<typename Comm>
      void combineDataMPI(Comm communicator) {
        CODI_UNUSED(communicator);
      }
#endif

      /// Get the allocated memory in bytes.
      double getAllocatedMemorySize() {
        return doubleData[allocatedMemoryIndex];
      }

      /// Get the used memory in bytes.
      double getUsedMemorySize() {
        return doubleData[usedMemoryIndex];
      }

      /// @}

    private:

      template<typename T>
      void performLocalReduction(T& lhs, T const& rhs, LocalReductionOperation operation) {
        switch (operation) {
          case LocalReductionOperation::Sum: {
            lhs += rhs;
            break;
          }
          case LocalReductionOperation::Max: {
            lhs = std::max(lhs, rhs);
            break;
          }
        }
      }

      template<typename T>
      void addEntryInternal(std::string const& name, EntryType const& type, LocalReductionOperation const& operation,
                            std::vector<T>& vector, T const& value) {
        size_t entryPos = vector.size();
        vector.push_back(value);

        if (sections.empty()) {
          addSection("General");
        }

        sections.back().data.push_back(Entry(name, type, operation, entryPos));
      }

      std::string formatEntry(Entry const& entry, int maximumFieldSize) const {
        return formatEntryFull(entry, true, maximumFieldSize);
      }

      std::string formatEntryFull(Entry const& entry, bool outputType, int maximumFieldSize) const {
        std::stringstream ss;

        switch (entry.type) {
          case EntryType::Double: {
            double formattedData = doubleData[entry.pos];
            std::string typeString = "";

            if (outputType) {
              formatSizeHumanReadable(formattedData, typeString);
            }
            ss << std::right << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(maximumFieldSize)
               << formattedData << typeString;
          } break;
          case EntryType::Long:
            ss << std::right << std::setw(maximumFieldSize) << longData[entry.pos];
            break;
          case EntryType::UnsignedLong:
            ss << std::right << std::setw(maximumFieldSize) << unsignedLongData[entry.pos];
            break;
          default:
            CODI_EXCEPTION("Unimplemented switch case.");
            break;
        }

        return ss.str();
      }

      size_t formatEntryLength(Entry const& entry) const {
        return formatEntryFull(entry, false, 0).size();
      }

      void formatSizeHumanReadable(double& size, std::string& type) const {
        char const* const typeList[] = {"B", "KB", "MB", "GB", "TB"};
        size_t typeListSize = sizeof(typeList) / sizeof(typeList[0]);

        size_t pos = 0;
        while (pos < typeListSize && size > 1024.0) {
          size /= 1024.0;
          pos += 1;
        }

        type = " ";
        type += typeList[pos];
      }

      size_t getMaximumNameLength() const {
        size_t maxLength = 0;
        for (Section const& section : sections) {
          for (Entry const& data : section.data) {
            maxLength = std::max(maxLength, data.name.size());
          }
        }

        return maxLength;
      }

      size_t getMaximumValueLength() const {
        size_t maxLength = 0;
        for (Section const& section : sections) {
          for (Entry const& data : section.data) {
            maxLength = std::max(maxLength, formatEntryLength(data));
          }
        }

        return maxLength;
      }
  };
}
