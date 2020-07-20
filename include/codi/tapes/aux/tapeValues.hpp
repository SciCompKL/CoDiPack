#pragma once

#include <iomanip>
#include <sstream>
#include <string>
#include <sstream>
#include <vector>


#include "../../config.h"
#include "../../aux/exceptions.hpp"

/** \copydoc codi::Namespace */
namespace codi {



  struct TapeValues {
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
          size_t pos;

          Entry() : name(), type(), pos() {}

          Entry(std::string const& name, EntryType const& type, size_t const& pos) :
            name(name),
            type(type),
            pos(pos) {}
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

      static double constexpr BYTE_TO_MB = 1.0 / 1024.0 / 1024.0;

      TapeValues(std::string const& tapeName) :
        sections(),
        doubleData(),
        longData(),
        unsignedLongData(),
        usedMemoryIndex(0),
        allocatedMemoryIndex(1) {

        addSection(tapeName);
        addEntryInternal("Total memory used", EntryType::Double, doubleData, 0.0);
        addEntryInternal("Total memory allocated", EntryType::Double, doubleData, 0.0);
      }

      void addDoubleEntry(std::string const& name, double const& value, bool usedMem = false, bool allocatedMem = false) {

        addEntryInternal(name, EntryType::Double, doubleData, value);

        if(usedMem) {
          doubleData[usedMemoryIndex] += value;
        }

        if(allocatedMem) {
          doubleData[allocatedMemoryIndex] += value;
        }
      }

      void addLongEntry(std::string const& name, long const& value) {
        addEntryInternal(name, EntryType::Long, longData, value);
      }

      void addSection(std::string const& name) {
        sections.push_back(Section(name));
      }

      void addUnsignedLongEntry(std::string const& name, unsigned long const& value) {
        addEntryInternal(name, EntryType::UnsignedLong, unsignedLongData, value);
      }

      void combineData() {
#ifdef MPI_VERSION
        MPI_Allreduce(MPI_IN_PLACE, doubleData.data(), doubleData.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, longData.data(), longData.size(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, unsignedLongData.data(), unsignedLongData.size(), MPI_UNSINGED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
      }

      template<typename Stream = std::ostream>
      void formatDefault(Stream& out = std::cout) const {

        std::string const hLine = "-------------------------------------\n";

        size_t maxNameSize = getMaximumNameLength();
        size_t maxValueSize = std::max((size_t)10, getMaximumValueLength());

        out << hLine;
        for(Section const& section : sections) {
          out << std::left << section.name << "\n";
          out << hLine;

          for(Entry const& entry : section.data) {
            out << "  "
                << std::left << std::setw(maxNameSize) << entry.name
                << " : "
                << formatEntry(entry, maxValueSize) << "\n";
          }

          if(!section.data.empty()) {
            out << hLine;
          }
        }
      }

      template<typename Stream = std::ostream>
      void formatHeader(Stream& out = std::cout) const {

        bool first = true;
        for(Section const& section : sections) {
          for(Entry const& entry : section.data) {

            if(first) {
              first = false;
            } else {
              out << "; ";
            }
            out << section.name << "-" << entry.name;
          }
        }

        out << "\n";
      }

      template<typename Stream = std::ostream>
      void formatRow(Stream& out = std::cout) const {

        size_t maxValueSize = std::max((size_t)10, getMaximumValueLength());

        bool first = true;
        for(Section const& section : sections) {
          for(Entry const& entry : section.data) {

            if(first) {
              first = false;
            } else {
              out << "; ";
            }
            out << formatEntry(entry, maxValueSize);
          }
        }

        out << "\n";
      }

      double getUsedMemorySize() {
        return doubleData[0];
      }

      double getAllocatedMemorySize() {
        return doubleData[1];
      }


    private:

      template<typename T>
      void addEntryInternal(std::string const& name, EntryType const& type, std::vector<T>& vector, T const& value) {
        size_t entryPos = vector.size();
        vector.push_back(value);

        if(sections.empty()) {
          addSection("General");
        }

        sections.back().data.push_back(Entry(name, type, entryPos));
      }


      std::string formatEntry(Entry const& entry, int maximumFieldSize) const {
        return formatEntryFull(entry, true, maximumFieldSize);
      }

      std::string formatEntryFull(Entry const& entry, bool outputType, int maximumFieldSize) const {
        std::stringstream ss;

        switch (entry.type) {
          case EntryType::Double:
            ss << std::right << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(maximumFieldSize) << doubleData[entry.pos];
            if(outputType) {
              ss << " MB";
            }
            break;
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

      size_t getMaximumNameLength() const {
        size_t maxLength = 0;
        for(Section const& section : sections) {
          for(Entry const& data : section.data) {
            maxLength = std::max(maxLength, data.name.size());
          }
        }

        return maxLength;
      }

      size_t getMaximumValueLength() const {
        size_t maxLength = 0;
        for(Section const& section : sections) {
          for(Entry const& data : section.data) {
            maxLength = std::max(maxLength, formatEntryLength(data));
          }
        }

        return maxLength;
      }
  };
}
