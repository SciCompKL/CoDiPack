#pragma once

#include <set>

#include "../../config.h"
#include "../aux/tapeConfiguration.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct DataManagementTapeInterface {
    public:

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void swap(DataManagementTapeInterface& other);
      void resetHard();
      void deleteAdjointVector();

      void writeToFile(const std::string& filename) const ;
      void readFromFile(const std::string& filename);
      void deleteData();

      std::set<ConfigurationOption> const& getAvailableOptions() const;
      size_t getOption(ConfigurationOption option) const;
      bool hasOption(ConfigurationOption option) const;
      void setOption(ConfigurationOption option, size_t value);
  };
}
