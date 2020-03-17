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

      void writeToFile(const std::string& filename);
      void readFromFile(const std::string& filename);
      void deleteData();

      std::set<ConfigurationOption> const& getAvailableOptions();
      size_t getOption(ConfigurationOption option);
      bool hasOption(ConfigurationOption option);
      void setOption(ConfigurationOption option, size_t value);
  };
}
