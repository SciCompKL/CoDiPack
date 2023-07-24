/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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

#include <set>

#include "../../config.h"
#include "../../traits/realTraits.hpp"
#include "../misc/tapeParameters.hpp"
#include "../misc/vectorAccessInterface.hpp"

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
   * parameters, they also allow the resizing of the underlying data. There are a few parameters that are read only and
   * a CODI_EXCEPTION is thrown if the user tries to set them. See the documentation of the enumerators in
   * #TapeParameters for details about each parameters.
   *
   * getParameter() and setParameter() will throw a CODI_EXCEPTION if the parameter is not defined for the tape. Which
   * parameters are defined can either be checked with hasParameter() or looked up in the list returned by
   * getAvailableParameters().
   *
   * \section vectorAccess Adjoint vector access
   * The function createVectorAccess() provides access to the internal vectors of the tape, usually the adjoint vector
   * and if available the primal value vector. If a generalized adjoint vector should be used, then the function
   * createCustomAdjointVectorAccess() is used. This is the same functionality that is also used in external functions
   * called from an evaluation with a custom adjoint vector.
   *
   * Instances of both methods have to be deleted with the deleteVectorAccess().
   *
   * Implementations may return different types that implement the same interface. Capturing these with auto may
   * improve the performance by eliminating virtual function calls.
   *
   * \section adjointMgmt Adjoint vector management
   * Tapes manage their internal adjoint vector automatically. This covers all routines offered by the tape itself. This
   * interface exposes parts of this adjoint vector management for external algorithms that build on top of a tape. See
   * also codi::AdjointsManagement.
   *
   * The functions resizeAdjointVector() and deleteAdjointVector() allow for memory optimizations. resizeAdjointVector()
   * can be used to guarantee a sufficient adjoint vector size for subsequent access without bounds checking.
   * deleteAdjointVector() frees the memory consumed by the adjoints.
   *
   * beginUseAdjointVector() and endUseAdjointVector() allow for guarding the adjoint vector against resizing, in a way
   * that is consistent with the internal adjoint vector safeguarding. See codi::InternalAdjointsInterface for a
   * description of the "in use" mechanism. In particular, the adjoint vector is "in use" whenever there is read or
   * write access to adjoint variables. As long as the adjoint vector is "in use", we cannot reallocate it. This is
   * important in multithreaded applications where multiple tapes compete for using and resizing the same adjoint
   * vector. Multiple threads can use the adjoint vector simultaneously. Attempts to use and resize the adjoint vector
   * from different threads will be resolved by means of this safeguarding mechanism. An attempt to resize the adjoint
   * vector from a thread while it has also declared usage results in a deadlock. The user of this interface is
   * responsible for avoiding this, that is, after a thread calls beginUseAdjointVector(), this thread must not call
   * tape methods that involve resizing and must not call resizeAdjointVector() until after a call to
   * endUseAdjointVector().
   *
   * \section misc Misc. functions
   * Some other functions for tape data management. Please see the function documentation.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Identifier>
  struct DataManagementTapeInterface {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See DataManagementTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See DataManagementTapeInterface.

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
      /// @name Interface: Adjoint vector access

      VectorAccessInterface<Real, Identifier>* createVectorAccess();  ///< See \ref vectorAccess.

      /// See \ref vectorAccess.
      template<typename Adjoint>
      VectorAccessInterface<Real, Identifier>* createVectorAccessCustomAdjoints(Adjoint* data);

      void deleteVectorAccess(VectorAccessInterface<Real, Identifier>* access);  ///< See \ref vectorAccess.

      /*******************************************************************************/
      /// @name Interface: Adjoint vector management

      void resizeAdjointVector();    ///< Explicitly trigger resizing of the adjoint vector. See \ref adjointMgmt.
      void deleteAdjointVector();    ///< Delete the adjoint vector. See \ref adjointMgmt.
      void beginUseAdjointVector();  ///< Declare that the adjoint vector is being used. See \ref adjointMgmt.
      void endUseAdjointVector();    ///< Declare that the adjoint vector is no longer used. See \ref adjointMgmt.

      /*******************************************************************************/
      /// @name Interface: Misc

      void swap(DataManagementTapeInterface& other);  ///< Swap all data with an other tape.

      /**
       * @brief Delete everything and return to the state after construction, as far as possible.
       *
       * Unlike other reset methods, this methods involves resizing the adjoint vector, this is not optional. Therefore,
       * no codi::AdjointsManagement parameter is offered.
       */
      void resetHard();
  };
}
