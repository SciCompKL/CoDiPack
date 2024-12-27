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

#include "../../config.h"
#include "primalBaseReaderWriter.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief  Writes a primal value tape in a binary format.
   *
   * The format of the produced <tt>"filename".dat</tt> file is as follows:
   *  \code
   * lhsIdentifier(Identifier) primalValue(Real) nPassiveValues(Config::ArgumentSize)
   * numberOfActiveArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) ... [repeats numberOfArguments times])
   * (passiveValues(Real) ... [repeats nPassiveValues times]) numberOfConstantArguments(Config::ArgumentSize)
   * (constantValues(Real) ... [repeats numberOfConstantArguments times]) evalHandleKey(EvalHandle)
   * \endcode
   *
   * The writer also generates an additional file <tt>"filename"IO.dat</tt> for the inputs and outputs, and a
   * <tt>"filename"Primals.dat</tt> for a sparse vector containing the primal values.
   *
   * A header file is produced that generates the evaluation handles in the new context <tt>"filename".hpp</tt>.
   *
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct PrimalBinaryTapeWriter : public PrimalBaseTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Tape = typename Type::Tape;                           ///< See TapeWriterInterface.

      using Identifier = typename Type::Identifier;  ///< See TapeWriterInterface.
      using Real = typename Type::Real;              ///< See TapeWriterInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See TapeWriterInterface.

      FILE* fileHandleBin = nullptr;  ///< File handle.

      /// Constructor.
      PrimalBinaryTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                             std::vector<Identifier> const& out)
          : PrimalBaseTapeWriter<T_Type>(name, in, out) {};

      /// \copydoc codi::TapeWriterInterface::start()
      void start(Tape& tape) {
        if (Tape::TapeTypes::IsStaticIndexHandler) {
          printPrimals(tape);
        }

        this->printIoBinary(tape);
        this->openFile(fileHandleBin, this->fileName, "wb");
      }

      /**
       * \copydoc codi::TapeWriterInterface::writeStatement(WriteInfo const&, Identifier const&, Real const&,
       *                                                    Config::ArgumentSize const&, size_t const&,
       *                                                    Identifier const* const, size_t const&,
       *                                                    Real const* const, size_t&, Real const* const,
       *                                                    EvalHandle)
       */
      void writeStatement(WriteInfo const& info, Identifier const& curLhsIdentifier, Real const& primalValue,
                          Config::ArgumentSize const& nPassiveValues, size_t const& curRhsIdentifiersPos,
                          Identifier const* const rhsIdentifiers, size_t const& curPassiveValuePos,
                          Real const* const passiveValues, size_t& curConstantPos, Real const* const constantValues,
                          EvalHandle stmtEvalHandle) {
        fwrite(&curLhsIdentifier, sizeof(Identifier), 1, fileHandleBin);
        fwrite(&primalValue, sizeof(Real), 1, fileHandleBin);
        fwrite(&nPassiveValues, sizeof(Config::ArgumentSize), 1, fileHandleBin);

        if (nPassiveValues == Config::StatementInputTag) CODI_Unlikely {
          // Do nothing.
        } else CODI_Likely {
          fwrite(&info.numberOfActiveArguments, sizeof(Config::ArgumentSize), 1, fileHandleBin);
          fwrite(&rhsIdentifiers[curRhsIdentifiersPos], sizeof(Identifier), info.numberOfActiveArguments,
                 fileHandleBin);
          fwrite(&passiveValues[curPassiveValuePos], sizeof(Real), nPassiveValues, fileHandleBin);
          fwrite(&info.numberOfConstantArguments, sizeof(Config::ArgumentSize), 1, fileHandleBin);
          fwrite(&constantValues[curConstantPos], sizeof(Real), info.numberOfConstantArguments, fileHandleBin);
        }
        // Register the statement expression (in the form of a string) for the evalHandleKey.

        size_t evalHandleIndex = this->getEvalHandleIndex(stmtEvalHandle, info.stmtExpression);
        fwrite(&evalHandleIndex, sizeof(size_t), 1, fileHandleBin);
      }

      /// Print the primal vector in sparse vector representation.
      void printPrimals(Tape& tape) {
        FILE* filePrimalBin = nullptr;
        size_t nPrimals = tape.getParameter(TapeParameters::PrimalSize);

        std::string fileNamePrimals = this->modifyFileName("Primals.dat");

        this->openFile(filePrimalBin, fileNamePrimals, "wb");

        // Print out the primals in a sparse vector representation.
        fwrite(&nPrimals, sizeof(size_t), 1, filePrimalBin);
        for (size_t indexCounter = 0; indexCounter < nPrimals; indexCounter++) {
          Real const curPrimal = tape.getPrimal(indexCounter);
          if (curPrimal != 0) {
            fwrite(&indexCounter, sizeof(Identifier), 1, filePrimalBin);
            fwrite(&curPrimal, sizeof(Real), 1, filePrimalBin);
          }
        }

        fclose(filePrimalBin);
      }

      /// \copydoc codi::TapeWriterInterface::finish()
      void finish() {
        // Generate the handleCreator .hpp file.
        this->generateHandleCreatorFile();

        fclose(fileHandleBin);
      }
  };

  /**
   * @brief  Reads and restores a binary file for a Primal tape
   *
   * The reader uses the <tt>"filename".dat</tt> file to restore the statements and the <tt>"filename"IO.dat</tt> to
   * restore the inputs and outputs of the tape.
   *
   * The format of the read <tt>"filename".dat</tt> file is as follows:
   * \code
   * lhsIdentifier(Identifier) primalValue(Real) nPassiveValues(Config::ArgumentSize)
   * numberOfActiveArguments(Config::ArgumentSize) (rhsIdentifiers(Identifier) ... [repeats numberOfArguments times])
   * (passiveValues(Real) ... [repeats nPassiveValues times]) numberOfConstantArguments(Config::ArgumentSize)
   * (constantValues(Real) ... [repeats numberOfConstantArguments times]) evalHandleKey(EvalHandle)
   * \endcode
   *
   * The primal values are read from the <tt>"filename"Primals.dat</tt> file.
   *
   * The EvalHandles are created in the <tt>"filename".hpp</tt> file, which should be included for the read process.
   *
   * See codi::TapeReaderInterface for a general description on how to use primal value tape readers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  struct PrimalBinaryTapeReader : public CommonBaseTapeReader<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeReaderInterface.
      using Tape = typename Type::Tape;                           ///< See TapeReaderInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeReaderInterface.
      using Real = typename Type::Real;                           ///< See TapeReaderInterface.

      std::vector<typename Tape::EvalHandle> evalHandles;  ///< Contains the unique evalHandles.
      /// Constructor.
      PrimalBinaryTapeReader(std::vector<typename Tape::EvalHandle> const& handles)
          : CommonBaseTapeReader<T_Type>(), evalHandles(handles) {}

      /// \copydoc codi::TapeReaderInterface::readFile()
      void readFile(std::string const& name) {
        FILE* fileHandleReadBin = nullptr;

        Identifier lhsIdentifier;
        Real primalValue;
        Config::ArgumentSize nPassiveValues;
        Config::ArgumentSize nActiveValues;
        std::vector<Identifier> rhsIdentifiers(Config::MaxArgumentSize, 0);
        std::vector<Real> rhsPrimalValues(Config::MaxArgumentSize, 0);
        Config::ArgumentSize nConstants;
        std::vector<Real> constants(Config::MaxArgumentSize, 0);

        EvalHandleKey evalHandleKey;
        this->fileName = name;

        // Restore IO and the primal vector, if the tape is a reuseType.
        this->restoreIoBinary();
        this->tape.getIndexManager().updateLargestCreatedIndex(this->largestIndex);

        if (Tape::TapeTypes::IsStaticIndexHandler) {
          restorePrimals();
        }

        this->openFile(fileHandleReadBin, this->fileName, "rb");

        while (fread(&lhsIdentifier, sizeof(Identifier), 1, fileHandleReadBin) == 1) {
          fread(&primalValue, sizeof(Real), 1, fileHandleReadBin);
          fread(&nPassiveValues, sizeof(Config::ArgumentSize), 1, fileHandleReadBin);

          if (nPassiveValues == Config::StatementLowLevelFunctionTag) CODI_Unlikely {
            // TODO.
          } else if (nPassiveValues == Config::StatementInputTag) CODI_Unlikely {
            // Do nothing.
          } else CODI_Likely {
            fread(&nActiveValues, sizeof(Config::ArgumentSize), 1, fileHandleReadBin);
            fread(rhsIdentifiers.data(), sizeof(Identifier), nActiveValues, fileHandleReadBin);
            fread(rhsPrimalValues.data(), sizeof(Real), nPassiveValues, fileHandleReadBin);
            fread(&nConstants, sizeof(Config::ArgumentSize), 1, fileHandleReadBin);
            fread(constants.data(), sizeof(Real), nConstants, fileHandleReadBin);
          }

          fread(&evalHandleKey, sizeof(EvalHandleKey), 1, fileHandleReadBin);
          this->tape.createStatementManual(lhsIdentifier, primalValue, nActiveValues, rhsIdentifiers.data(),
                                           nPassiveValues, rhsPrimalValues.data(), nConstants, constants.data(),
                                           evalHandles[evalHandleKey]);
        }

        fclose(fileHandleReadBin);
      }

      /// Read and record the primal values from the sparse primal vector stored in <tt>"filename"Primals.dat</tt>
      void restorePrimals() {
        FILE* filePrimalBin = nullptr;

        size_t nPrimals;
        Identifier curIdentifier;
        Real curPrimal;

        this->openFile(filePrimalBin, this->modifyFileName("Primals.dat"), "rb");

        // Restore the primal vector from a sparse vector.
        fread(&nPrimals, sizeof(size_t), 1, filePrimalBin);
        // Resize the primal vector.
        this->tape.setParameter(TapeParameters::PrimalSize, nPrimals);

        while (fread(&curIdentifier, sizeof(Identifier), 1, filePrimalBin) == 1) {
          fread(&curPrimal, sizeof(Real), 1, filePrimalBin);
          this->tape.setPrimal(curIdentifier, curPrimal);
        }
      }
  };
}
