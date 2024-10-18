/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <iostream>
#include <memory>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../data/position.hpp"
#include "../io/tapeReaderWriterInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Write a tape to a file using a writer from the \ref codi::TapeWriterInterface. When reading a stored tape,
   * create a primal or a Jacobian statement on the tape.
   *
   * See \ref codi::TapeWriterInterface and \ref codi::TapeReaderInterface for details.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape usually, chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier, typename T_Position>
  struct ReadWriteTapeInterface {
    public:
      using Real = CODI_DD(T_Real, double);                 ///< See ReadWriteTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);         ///< See ReadWriteTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);        ///< See ReadWriteTapeInterface.
      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See ReadWriteTapeInterface.

      using Type = CODI_ANY;     ///< Provided by the template arguments in the writeTape methods.
      using EvalHandle = void*;  ///< See PrimalValueTapeTypes.

      using WriterInterfaceUniquePtr = std::unique_ptr<TapeWriterInterface<Type>>;  ///< Unique pointer to writer.
      using WriterInterface = TapeWriterInterface<Type>;                            ///< Writer interface.

      /*******************************************************************************/
      /// @name Tape writing

      /// For full-tape writers using a smart pointer.
      template<typename Type>
      void writeTape(WriterInterfaceUniquePtr writer);
      /// For partial-tape writers using a smart pointer.
      template<typename Type>
      void writeTape(WriterInterfaceUniquePtr writer, Position const& start, Position const& end);

      /// For full-tape writers using a manually generated writers.
      template<typename Type>
      void writeTape(WriterInterface& writeToFileRef);
      /// For partial-tape writers using a manually generated writers.
      template<typename Type>
      void writeTape(WriterInterface& writeToFileRef, Position const& start, Position const& end);

      /// For full or partial tapes using a pointer to the writer.
      template<typename Type>
      void writeTape(WriterInterface* writer, Position const& start, Position const& end);
      /*******************************************************************************/
      /// @name Tape Reading

      /// Initialize the statement from a file. It is especially important that the lhsIndex is valid when using this
      /// method. This overload is used for Jacobian tapes.
      void createStatementManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size,
                                 Real const* jacobians, Identifier const* rhsIdentifiers);

      ///  Initialize the statement and the rhs vectors from a file. It is especially important that the Identifiers are
      ///  valid when using this method. This overload is used for Primal tapes.
      void createStatementManual(Identifier const& lhsIndex, Real const& lhsValue,
                                 Config::ArgumentSize const& nActiveValues, Identifier const* const rhsIdentifiers,
                                 Config::ArgumentSize const& nPassiveValues, Real const* const rhsPrimals,
                                 Config::ArgumentSize const& nConstants, Real const* const rhsConstant,
                                 EvalHandle const& evalHandle);
  };
}
