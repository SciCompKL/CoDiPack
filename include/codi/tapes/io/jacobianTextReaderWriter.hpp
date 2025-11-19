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

#include "../../config.h"
#include "jacobianBaseReaderWriter.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief  Writes a Jacobian tape in a text format.
   *
   * The format of the produced <tt>"filename".txt</tt> file is as follows:
   * \code
   * lhsIdentifier(Identifier) numberOfArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) rhsJacobian(Real) ...
   * [repeats numberOfArguments times]) \endcode
   *
   * The writer also generates an additional file <tt>"filename"IO.txt</tt> for the inputs and outputs of the tape.
   *
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct JacobianTextTapeWriter : public CommonBaseTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Tape = typename Type::Tape;                           ///< See TapeWriterInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeWriterInterface.
      using Real = typename Type::Real;                           ///< See TapeWriterInterface.

      FILE* fileHandleTxt = nullptr;  ///< File handle.

      bool printIO;                      ///< Flag to enable the writing of the IO file.
      bool printColumnNames;             ///< Flag to enable column names.
      bool printInputStatements = true;  ///< Flag to enable the output of input statements.

      /// Constructor.
      JacobianTextTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                             std::vector<Identifier> const& out, bool ifIO, bool ifColumnNames)
          : CommonBaseTapeWriter<T_Type>(name, in, out), printIO(ifIO), printColumnNames(ifColumnNames) {};

      /// Set if input statements should be printed.
      void setInputStatementOutput(bool value) {
        printInputStatements = value;
      }

      /// \copydoc codi::TapeWriterInterface::start()
      void start(Tape& tape) {
        if (printIO) {
          this->printIoText(tape);
        }

        this->openFile(fileHandleTxt, this->fileName, "w");
        if (printColumnNames) {
          fprintf(fileHandleTxt, "|  LHS Index  |  # of Args  |  RHS Indices  | RHS Jacobian Values |");
        }
      }

      /**
       * \copydoc codi::TapeWriterInterface::writeStatement(Identifier const&, size_t&, Real const* const,
       *                                                    Identifier const* const, Config::ArgumentSize const&)
       */
      void writeStatement(Identifier const& curLhsIdentifier, size_t& curJacobianPos, Real const* const rhsJacobians,
                          Identifier const* const rhsIdentifiers, Config::ArgumentSize const& nJacobians) {
        if (nJacobians == Config::StatementInputTag) CODI_Unlikely {
          if (printInputStatements) {
            fprintf(fileHandleTxt, "\n%d  %hhu  []", curLhsIdentifier, static_cast<Config::ArgumentSize>(nJacobians));
          }
        } else CODI_Likely {
          fprintf(fileHandleTxt, "\n%d  %hhu  [", curLhsIdentifier, static_cast<Config::ArgumentSize>(nJacobians));

          for (size_t argCount = 0; argCount < nJacobians; argCount++) {
            fprintf(fileHandleTxt, " %d ", rhsIdentifiers[curJacobianPos + argCount]);
          }

          fprintf(fileHandleTxt, "]  [");
          for (size_t argCount = 0; argCount < nJacobians; argCount++) {
            fprintf(fileHandleTxt, " %0.12e ", rhsJacobians[curJacobianPos + argCount]);
          }

          fprintf(fileHandleTxt, "]");
        }
      }

      /// \copydoc codi::TapeWriterInterface::finish()
      void finish() {
        fclose(fileHandleTxt);
      }
  };

  /**
   * @brief  Reads and restores a text file for a Jacobian tape
   *
   * The reader uses the <tt>"filename".txt</tt> file to restore the statements and the <tt>"filename"IO.txt</tt> to
   * restore the inputs and outputs of the tape.
   *
   * The format of the read <tt>"filename".txt</tt> file is as follows:
   * \code
   * lhsIdentifier(Identifier) numberOfArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) rhsJacobian(Real) ...
   * [repeats numberOfArguments times]) \endcode
   *
   * See codi::TapeReaderInterface for a general description on how to use tape readers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  struct JacobianTextTapeReader : public JacobianBaseTapeReader<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeReaderInterface.
      using Tape = typename Type::Tape;                           ///< See TapeReaderInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeReaderInterface.
      using Real = typename Type::Real;                           ///< See TapeReaderInterface.

      JacobianTextTapeReader() : JacobianBaseTapeReader<T_Type>() {}  ///< Constructor.

      /// \copydoc codi::TapeReaderInterface::readFile()
      void readFile(std::string const& name) {
        FILE* fileHandleReadTxt = nullptr;

        bool isFirstIdentifier = true;
        Identifier lowestIndex = 0;

        Identifier lhsIdentifier;
        Config::ArgumentSize nArgs;

        std::vector<Identifier> rhsIdentifiers(Config::MaxArgumentSize, 0);
        std::vector<Real> rhsJacobians(Config::MaxArgumentSize, 0);
        this->fileName = name;

        this->restoreIoText();
        this->tape.getIndexManager().updateLargestCreatedIndex(this->largestIndex);

        this->openFile(fileHandleReadTxt, this->fileName, "r+");

        // Read and discard column titles.
        fscanf(fileHandleReadTxt, "|  LHS Index  |  # of Args  |  RHS Indices  | RHS Jacobian Values |");

        // Read the statements.
        while (fscanf(fileHandleReadTxt, "\n%d  %hhu  [", &lhsIdentifier, &nArgs) == 2) {
          if (nArgs == Config::StatementLowLevelFunctionTag) CODI_Unlikely {
            // TODO.
          } else if (nArgs == Config::StatementInputTag) CODI_Unlikely {
            // Do nothing.
          } else CODI_Likely {
            for (size_t argumentCount = 0; argumentCount < nArgs; argumentCount++) {
              fscanf(fileHandleReadTxt, " %d ", &rhsIdentifiers[argumentCount]);
            }
            fscanf(fileHandleReadTxt, "]  [");
            for (size_t argumentCount = 0; argumentCount < nArgs; argumentCount++) {
              fscanf(fileHandleReadTxt, " %lf ", &rhsJacobians[argumentCount]);
            }
          }
          fscanf(fileHandleReadTxt, "]");

          this->registerStatement(lhsIdentifier, nArgs, rhsIdentifiers, rhsJacobians, lowestIndex, isFirstIdentifier);
        }

        /* Update the user provided IO with a potential offset in the linear case. The registerStatement() returns
           the offset.*/
        this->updateUserIO(lowestIndex);

        fclose(fileHandleReadTxt);
      }
  };
}
