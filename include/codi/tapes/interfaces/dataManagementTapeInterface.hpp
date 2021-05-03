#pragma once

#include <set>

#include "../../config.h"
#include "../aux/tapeParameters.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Provides file IO, information about internal tape vectors and allows to clear tape data.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * \section fileIO File IO functions
   * This interface offers advanced data management capabilities for the tape. The file IO routines provide
   * the capability to write the internal tape data to the disk. The goal of this is moving the tape temporarily from
   * RAM to disk. After writing the tape with writeToFile(), a call to deleteData() ensures that all internal data that
   * was written to disk is freed so that the RAM footprint is minimized. Usually, neither management data nor external
   * function data are exported. This means that the same tape that called writeToFile() has to call readFromFile() and
   * that offloaded tapes are not meaningful across multiple executions of the application.
   *
   * \section parameters Parameters functions
   * The parameter functions provide access to the sizes of the internal tape implementations. For most of the
   * parameters, they also allow the resizing of the underlying data. There are a few parameters that are read only and a
   * CODI_EXCEPTION is thrown if the user tries to set them. See the documentation of the enumerators in #TapeParameters
   * for details about each parameters.
   *
   * getParameter() and setParameter() will throw a CODI_EXCEPTION if the parameter is not defined for the tape. Which
   * parameters are defined can either be checked with hasParameter() or looked up in the list returned by
   * getAvailableParameters().
   *
   * \section misc Misc. functions
   * Some other functions for tape data management. Please see the function documentation.
   */
  struct DataManagementTapeInterface {
    public:

      /*******************************************************************************/
      /// @name Interface: File IO

      void writeToFile(std::string const& filename) const;  ///< See \ref fileIO.
      void readFromFile(std::string const& filename);       ///< See \ref fileIO.
      void deleteData();                                    ///< See \ref fileIO.

      /*******************************************************************************/
      /// @name Interface: Parameters

      std::set<TapeParameters> const& getAvailableParameters() const;  ///< See \ref parameters.
      size_t getParameter(TapeParameters parameter) const;             ///< See \ref parameters.
      bool hasParameter(TapeParameters parameter) const;               ///< See \ref parameters.
      void setParameter(TapeParameters parameter, size_t value);       ///< See \ref parameters.

      /*******************************************************************************/
      /// @name Interface: Misc

      void swap(DataManagementTapeInterface& other);  ///< Swap all data with an other tape.
      void resetHard();  ///< Delete everything and return to the state after construction, as far as possible.
      void deleteAdjointVector();  ///< Delete the adjoint vector.
  };
}
