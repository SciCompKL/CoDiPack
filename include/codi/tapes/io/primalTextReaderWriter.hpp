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
#include "primalBaseReaderWriter.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief  Writes a primal value tape in a text format.
   *
   * The format of the produced <tt>"filename".txt</tt> file is as follows:
   *  \code
   * lhsIdentifier(Identifier) primalValue(Real) nPassiveValues(Config::ArgumentSize)
   * numberOfActiveArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) ... [repeats numberOfArguments times])
   * (passiveValues(Real) ... [repeats nPassiveValues times]) numberOfConstantArguments(Config::ArgumentSize)
   * (constantValues(Real) ... [repeats numberOfConstantArguments times]) evalHandleKey(EvalHandle)
   * \endcode
   *
   * The writer also generates an additional file <tt>"filename"IO.txt</tt> for the inputs and outputs, and a
   * <tt>"filename"Primals.txt</tt> for a sparse vector containing the primal values.
   *
   * A header file is produced that generates the evaluation handles in the new context <tt>"filename".hpp</tt>.
   *
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct PrimalTextTapeWriter : public PrimalBaseTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Tape = typename Type::Tape;                           ///< See TapeWriterInterface.

      using Identifier = typename Type::Identifier;  ///< See TapeWriterInterface.
      using Real = typename Type::Real;              ///< See TapeWriterInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See TapeWriterInterface.

      FILE* fileHandleTxt = nullptr;  ///< File handle.

      bool printIO;           ///< Flag to enable the writing of the IO file.
      bool printColumnNames;  ///< Flag to enable column names.

      /// Constructor.
      PrimalTextTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                           std::vector<Identifier> const& out, bool ifIO, bool ifColumnNames)
          : PrimalBaseTapeWriter<T_Type>(name, in, out), printIO(ifIO), printColumnNames(ifColumnNames) {}

      /// \copydoc codi::TapeWriterInterface::start()
      void start(Tape& tape) {
        if (Tape::TapeTypes::IsStaticIndexHandler) {
          printPrimals(tape);
        }

        if (printIO) {
          this->printIoText(tape);
        }

        this->openFile(fileHandleTxt, this->fileName, "w");

        if (printColumnNames) {
          fprintf(fileHandleTxt,
                  "| # of LHS Args | LHS Indices | Primal Values | # of Passive Args | # of Active Args | RHS Indices "
                  "| RHS Primal Values | # of Constants | Constants | Statement Key |");
        }
      }

      /**
       * \copydoc codi::TapeWriterInterface::writeStatement(WriteInfo const&, Identifier const*, Real const*,
       *                                                    Config::ArgumentSize const&, Identifier const* const,
       *                                                    Real const* const, Real const* const, EvalHandle)
       */
      void writeStatement(WriteInfo const& info, Identifier const* lhsIdentifiers, Real const* primalValues,
                          Config::ArgumentSize const& nPassiveValues, Identifier const* const rhsIdentifiers,
                          Real const* const passiveValues, Real const* const constantValues,
                          EvalHandle stmtEvalHandle) {
        // nOutputArguments
        fprintf(fileHandleTxt, "\n%hhu [", static_cast<Config::ArgumentSize>(info.numberOfOutputArguments));

        // Lhs Identifiers.
        for (size_t lhsCount = 0; lhsCount < info.numberOfOutputArguments; lhsCount++) {
          fprintf(fileHandleTxt, " %d ", lhsIdentifiers[lhsCount]);
        }
        fprintf(fileHandleTxt, "]  [");

        // Lhs primal values.
        for (size_t lhsCount = 0; lhsCount < info.numberOfOutputArguments; lhsCount++) {
          fprintf(fileHandleTxt, " %0.12e ", primalValues[lhsCount]);
        }
        fprintf(fileHandleTxt, "] ");

        // nPassiveValues
        fprintf(fileHandleTxt, "%hhu ", static_cast<Config::ArgumentSize>(nPassiveValues));

        if (nPassiveValues == Config::StatementInputTag) CODI_Unlikely {
          // Do nothing.
        } else CODI_Likely {
          // nActiveValues.
          fprintf(fileHandleTxt, " %hhu  [", static_cast<Config::ArgumentSize>(info.numberOfActiveArguments));

          // Rhs Identifiers.
          for (size_t rhsCount = 0; rhsCount < info.numberOfActiveArguments; rhsCount++) {
            fprintf(fileHandleTxt, " %d ", rhsIdentifiers[rhsCount]);
          }
          fprintf(fileHandleTxt, "]  [");

          // Rhs Passive Primal Values.
          for (size_t passiveCount = 0; passiveCount < nPassiveValues; passiveCount++) {
            fprintf(fileHandleTxt, " %0.12e ", passiveValues[passiveCount]);
          }

          // nConstants.
          fprintf(fileHandleTxt, "]  %hhu  [", static_cast<Config::ArgumentSize>(info.numberOfConstantArguments));

          // Rhs Constants.
          for (size_t constantCount = 0; constantCount < info.numberOfConstantArguments; constantCount++) {
            fprintf(fileHandleTxt, " %0.12e ", constantValues[constantCount]);
          }

          fprintf(fileHandleTxt, "]");
        }

        size_t evalHandleIndex = this->getEvalHandleIndex(stmtEvalHandle, info.stmtExpression);

        // Stmt handle.
        fprintf(fileHandleTxt, "  [ %zu ]", evalHandleIndex);
      }

      /// Print the primal vector in sparse vector representation.
      void printPrimals(Tape& tape) {
        FILE* filePrimalTxt = nullptr;
        size_t nPrimals = tape.getParameter(TapeParameters::PrimalSize);

        std::string fileNamePrimals = this->modifyFileName("Primals.txt");

        this->openFile(filePrimalTxt, fileNamePrimals, "w");

        // Print out the primals in a sparse vector representation.
        fprintf(filePrimalTxt, "%zu\n", nPrimals);
        for (size_t indexCounter = 0; indexCounter < nPrimals; indexCounter++) {
          Real const curPrimal = tape.getPrimal(indexCounter);
          if (curPrimal != 0) {
            fprintf(filePrimalTxt, "%zu ", indexCounter);
            fprintf(filePrimalTxt, "%0.12e\n", curPrimal);
          }
        }

        fclose(filePrimalTxt);
      }

      /// \copydoc codi::TapeWriterInterface::finish()
      void finish() {
        this->generateHandleCreatorFile();
        fclose(fileHandleTxt);
      }
  };

  /**
   * @brief  Reads and restores a text file for a primal value tape
   *
   * The reader uses the <tt>"filename".txt</tt> file to restore the statements and the <tt>"filename"IO.txt</tt> to
   * restore the inputs and outputs of the tape.
   *
   * The format of the read <tt>"filename".txt</tt> file is as follows:
   * \code
   * lhsIdentifier(Identifier) primalValue(Real) nPassiveValues(Config::ArgumentSize)
   * numberOfActiveArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) ... [repeats numberOfArguments times])
   * (passiveValues(Real) ... [repeats nPassiveValues times]) numberOfConstantArguments(Config::ArgumentSize)
   * (constantValues(Real) ... [repeats numberOfConstantArguments times]) evalHandleKey(EvalHandle)
   * \endcode
   *
   * The primal values are read from the <tt>"filename"Primals.txt</tt> file.
   *
   * The EvalHandles are created in the <tt>"filename".hpp</tt> file, which should be included for the read process.
   *
   * See codi::TapeReaderInterface for a general description on how to use primal value tape readers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  struct PrimalTextTapeReader : public CommonBaseTapeReader<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Tape = typename Type::Tape;                           ///< See TapeReaderInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeReaderInterface.
      using Real = typename Type::Real;                           ///< See TapeReaderInterface.

      std::vector<typename Tape::EvalHandle> evalHandles;  ///< Contains the unique evalHandles.

      PrimalTextTapeReader(std::vector<typename Tape::EvalHandle> const& handles)
          : CommonBaseTapeReader<T_Type>(), evalHandles(handles) {}  ///< Constructor.

      /// \copydoc codi::TapeReaderInterface::readFile()
      void readFile(std::string const& name) {
        FILE* fileHandleReadTxt = nullptr;

        Config::ArgumentSize nOutputValues;
        std::vector<Identifier> lhsIdentifiers(Config::MaxArgumentSize, 0);
        std::vector<Real> primalValues(Config::MaxArgumentSize, 0.0);
        Config::ArgumentSize nPassiveValues;
        Config::ArgumentSize nActiveValues;
        std::vector<Identifier> rhsIdentifiers(Config::MaxArgumentSize, 0);
        std::vector<Real> rhsPrimalValues(Config::MaxArgumentSize, 0);
        Config::ArgumentSize nConstants;
        std::vector<Real> constants(Config::MaxArgumentSize, 0);

        EvalHandleKey evalHandleKey;
        this->fileName = name;

        // Restore IO and Primals
        this->restoreIoText();
        this->tape.getIndexManager().updateLargestCreatedIndex(this->largestIndex);

        if (Tape::TapeTypes::IsStaticIndexHandler) {
          restorePrimals();
        }
        this->openFile(fileHandleReadTxt, this->fileName, "r+");

        // Read and discard column titles
        fscanf(fileHandleReadTxt,
               "| # of LHS Args | LHS Indices | Primal Values | # of Passive Args | # of Active Args | RHS Indices "
               "| RHS Primal Values | # of Constants | Constants | Statement Key |");
        // Read Statements
        while (fscanf(fileHandleReadTxt, "\n%hhu [", &nOutputValues) == 1) {
          // Lhs identifiers
          for (size_t argumentCount = 0; argumentCount < nOutputValues; argumentCount++) {
            fscanf(fileHandleReadTxt, " %d ", &lhsIdentifiers[argumentCount]);
          }
          fscanf(fileHandleReadTxt, "]  [");

          // Lhs primal values
          for (size_t argumentCount = 0; argumentCount < nOutputValues; argumentCount++) {
            fscanf(fileHandleReadTxt, " %lf ", &primalValues[argumentCount]);
          }

          // nPassive values
          fscanf(fileHandleReadTxt, "] %hhu ", &nPassiveValues);

          if (nPassiveValues == Config::StatementLowLevelFunctionTag) CODI_Unlikely {
            // TODO.
          } else if (nPassiveValues == Config::StatementInputTag) CODI_Unlikely {
            // Do nothing.
          } else CODI_Likely {
            fscanf(fileHandleReadTxt, " %hhu  [", &nActiveValues);

            // Rhs identifiers
            for (size_t argumentCount = 0; argumentCount < nActiveValues; argumentCount++) {
              fscanf(fileHandleReadTxt, " %d ", &rhsIdentifiers[argumentCount]);
            }
            fscanf(fileHandleReadTxt, "]  [");
            // Rhs primals
            for (size_t argumentCount = 0; argumentCount < nPassiveValues; argumentCount++) {
              fscanf(fileHandleReadTxt, " %lf ", &rhsPrimalValues[argumentCount]);
            }
            fscanf(fileHandleReadTxt, "]  %hhu  [", &nConstants);
            // Constants
            for (size_t argumentCount = 0; argumentCount < nConstants; argumentCount++) {
              fscanf(fileHandleReadTxt, " %lf ", &constants[argumentCount]);
            }
            fscanf(fileHandleReadTxt, "]");
          }
          // EvalHandleKey
          fscanf(fileHandleReadTxt, "  [ %zu ]", &evalHandleKey);
          this->tape.createStatementManual(nOutputValues, lhsIdentifiers.data(), primalValues.data(), nActiveValues,
                                           rhsIdentifiers.data(), nPassiveValues, rhsPrimalValues.data(), nConstants,
                                           constants.data(), evalHandles[evalHandleKey]);
        }

        fclose(fileHandleReadTxt);
      }

      /// Read and record the primal values from the sparse primal vector stored in <tt>"filename"Primals.txt</tt>
      void restorePrimals() {
        FILE* filePrimalTxt = nullptr;

        size_t nPrimals;
        Identifier curIdentifier;
        Real curPrimal;

        this->openFile(filePrimalTxt, this->modifyFileName("Primals.txt"), "r+");

        // Restore the primal vector from a sparse vector.
        fscanf(filePrimalTxt, "%zu\n", &nPrimals);

        // Resize the primal vector.
        this->tape.setParameter(TapeParameters::PrimalSize, nPrimals);

        while (fscanf(filePrimalTxt, "%d ", &curIdentifier) == 1) {
          fscanf(filePrimalTxt, "%lf\n", &curPrimal);
          this->tape.setPrimal(curIdentifier, curPrimal);
        }
      }
  };
}
