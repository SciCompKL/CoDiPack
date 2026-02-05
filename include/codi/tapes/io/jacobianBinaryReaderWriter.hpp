/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "jacobianBaseReaderWriter.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief  Writes a Jacobian tape in a binary format.
   *
   * The writer also generates an additional file <tt>"filename"IO.dat</tt> for the inputs and outputs of the tape.
   *
   * The format of the produced <tt>"filename".dat</tt> file is as follows:
   * \code
   * lhsIdentifier(Identifier) numberOfArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) rhsJacobian(Real) ...
   * [repeats numberOfArguments times]) \endcode
   *
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct JacobianBinaryTapeWriter : public CommonBaseTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Tape = typename Type::Tape;                           ///< See TapeWriterInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeWriterInterface.
      using Real = typename Type::Real;                           ///< See TapeWriterInterface.

      FILE* fileHandleBin = nullptr;  ///< File handle.

      /// Constructor.
      JacobianBinaryTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                               std::vector<Identifier> const& out)
          : CommonBaseTapeWriter<T_Type>(name, in, out) {};

      /// \copydoc codi::TapeWriterInterface::start()
      void start(Tape& tape) {
        this->printIoBinary(tape);

        this->openFile(fileHandleBin, this->fileName, "wb");
      }

      /**
       * \copydoc codi::TapeWriterInterface::writeStatement(Identifier const&, size_t&, Real const* const,
       *                                                    Identifier const* const, Config::ArgumentSize const&)
       */
      void writeStatement(Identifier const& curLhsIdentifier, size_t& curJacobianPos, Real const* const rhsJacobians,
                          Identifier const* const rhsIdentifiers, Config::ArgumentSize const& nJacobians) {
        fwrite(&curLhsIdentifier, sizeof(Identifier), 1, fileHandleBin);
        fwrite(&nJacobians, sizeof(Config::ArgumentSize), 1, fileHandleBin);
        if (nJacobians == Config::StatementInputTag) CODI_Unlikely {
          // Do nothing.
        } else CODI_Likely {
          fwrite(&rhsIdentifiers[curJacobianPos], sizeof(Identifier), nJacobians, fileHandleBin);
          fwrite(&rhsJacobians[curJacobianPos], sizeof(Real), nJacobians, fileHandleBin);
        }
      }

      /// \copydoc codi::TapeWriterInterface::finish()
      void finish() {
        fclose(fileHandleBin);
      }
  };

  /**
   * @brief  Reads and restores a binary file for a Jacobian tape
   *
   * The reader uses the <tt>"filename".dat</tt> file to restore the statements and the <tt>"filename"IO.dat</tt> to
   * restore the inputs and outputs of the tape.
   *
   * The format of the read <tt>"filename".dat</tt> file is as follows:
   * \code
   * lhsIdentifier(Identifier) numberOfArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) rhsJacobian(Real) ...
   * [repeats numberOfArguments times]) \endcode
   *
   * See codi::TapeReaderInterface for a general description on how to use tape readers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  struct JacobianBinaryTapeReader : public JacobianBaseTapeReader<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeReaderInterface.
      using Tape = typename Type::Tape;                           ///< See TapeReaderInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeReaderInterface.
      using Real = typename Type::Real;                           ///< See TapeReaderInterface.

      JacobianBinaryTapeReader() : JacobianBaseTapeReader<T_Type>() {}  ///< Constructor.

      /// \copydoc codi::TapeReaderInterface::readFile()
      void readFile(std::string const& name) {
        FILE* fileHandleReadBin = nullptr;

        bool isFirstIdentifier = true;
        Identifier lowestIndex = 0;

        Identifier lhsIdentifier;
        Config::ArgumentSize nArgs;

        std::vector<Identifier> rhsIdentifiers(Config::MaxArgumentSize, 0);
        std::vector<Real> rhsJacobians(Config::MaxArgumentSize, 0);

        this->fileName = name;
        this->restoreIoBinary();
        this->tape.getIndexManager().updateLargestCreatedIndex(this->largestIndex);

        this->openFile(fileHandleReadBin, this->fileName, "rb");

        while (fread(&lhsIdentifier, sizeof(Identifier), 1, fileHandleReadBin) == 1) {
          fread(&nArgs, sizeof(Config::ArgumentSize), 1, fileHandleReadBin);
          if (nArgs == Config::StatementLowLevelFunctionTag) CODI_Unlikely {
            // TODO.
          } else if (nArgs == Config::StatementInputTag) CODI_Unlikely {
            // Do nothing.
          } else CODI_Likely {
            fread(rhsIdentifiers.data(), sizeof(Identifier), nArgs, fileHandleReadBin);
            fread(rhsJacobians.data(), sizeof(Real), nArgs, fileHandleReadBin);
          }
          this->registerStatement(lhsIdentifier, nArgs, rhsIdentifiers, rhsJacobians, lowestIndex, isFirstIdentifier);
        }

        /* Update the user provided IO with a potential offset in the linear case. The registerStatement() returns
           the offset.*/
        this->updateUserIO(lowestIndex);

        fclose(fileHandleReadBin);
      }
  };
}
