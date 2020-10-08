#pragma once

#include <set>

#include "../../config.h"
#include "../aux/tapeParameters.hpp"

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

      void writeToFile(std::string const& filename) const ;
      void readFromFile(std::string const& filename);
      void deleteData();

      std::set<TapeParameters> const& getAvailableParameters() const;
      size_t getParameter(TapeParameters parameter) const;
      bool hasParameter(TapeParameters parameter) const;
      void setParameter(TapeParameters parameter, size_t value);
  };
}
