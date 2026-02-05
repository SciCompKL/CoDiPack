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
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../data/position.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  using EvalHandleKey = size_t;  ///< Key for the evalHandle lookup.

  /// Used to select the type of writer that should be generated.
  enum class FileType {
    Text,
    Binary,
    Graph,
    Math,
    Invalid
  };

  /**
   * @brief Used by the math statements to record the type of each identifier. This information is then used in the math
   * representation.
   */
  enum class IdentifierType {
    Input,
    Output,
    Temp
  };

  /**
   * @brief This class is used during the writing process of a primal value tape. The WriteInfo is returned by
   * StatementEvaluator call with StatementCall::WriteInformation.
   */
  struct WriteInfo {
      size_t numberOfOutputArguments;    ///< Number of output arguments.
      size_t numberOfActiveArguments;    ///< Number of active arguments.
      size_t numberOfConstantArguments;  ///< Number of constant arguments.
      std::string stmtExpression;        ///< Used to generate a .hpp file for reading back a primal value tape.
      std::string mathRepresentation;  ///< The math representation of a statement used in codi::PrimalValueGraphWriter
                                       ///< and in the codi::MathRepWriter.
  };

  /**
   * @brief  The interface used by all the tape writers. The tape writers are used to generate text, binary, graphical
   * or the math statement files from the current tape. The text and binary files can be used to restore a tape using
   * the \ref codi::TapeReaderInterface in a new context. The writers are supported for both Jacobian and primal value
   * tapes.
   *
   * The user does not access the methods of this interface directly, instead the standard steps are as follows:
   *
   * 1) Record the inputs and outputs of the tape in std::vectors.
   *
   * 2) Use the \ref codi::createWriter() function with a fileName, the input and output vectors and the desired
   * Codi::FileType to create an automatic writer.
   *
   * 3) Call the \ref codi::ReadWriteTapeInterface::writeTape() method with the generated writer.
   *
   * 4) The \ref codi::TapeWriterInterface::start() is called once, after which the relevant \ref
   * codi::TapeWriterInterface::writeStatement() method is called for each of the statements on the tape. After all the
   * statements have been added, the \ref codi::TapeWriterInterface::finish method() is called and the file has been
   * generated.
   *
   * \code{.cpp}
   *   std::vector inputs, outputs;
   *
   *   // **Follow the normal procedure to add statements to the tape, to register inputs and outputs, and to evaluate
   *   // the tape.
   *
   *   // Generate a text writer for the current tape and call the write tape method to create the .txt file
   *   tape.writeTape(codi::createWriter<Real>("example.txt", inputs, outputs, codi::FileType::Text));
   *
   * \endcode
   *
   * The createWriter function does not have to be used. Writers can instead be created and managed manually. The manual
   * writers are then passed to the writeTape method. The rest of the procedure remains the same.
   *
   * \code{.cpp}
   *   std::vector inputs, outputs;
   *
   *   // **Follow the normal procedure to add statements to the tape, to register inputs and outputs, and to evaluate
   *   // the tape.
   *
   *   // Manually create a binary writer for a Jacobian tape.
   *   codi::JacobianBinaryWriter<Real> writerHandle{"example.dat", inputs, outputs};
   *
   *   // Generate a binary file for the current tape.
   *   tape.writeTape(writerHandle);
   *
   * \endcode
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct TapeWriterInterface {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< The evaluation type.
      using Tape = typename Type::Tape;                           ///< The tape type that is to be written out.

      using Identifier = typename Type::Identifier;  ///< Identifier for the internal management, e.g. int.
      using Real = typename Type::Real;              ///< Primal computation type, e.g. double.

      using EvalHandle = typename Tape::EvalHandle;  ///< Evaluation handle used for primal value tapes.

      virtual ~TapeWriterInterface() {};  ///< Destructor.

      /**
       * @brief  Called once at the beginning of the tape write process. Should initialize all required data structures
       * and files.
       */
      virtual void start(Tape& tape) {
        CODI_UNUSED(tape);
      }

      /**
       * @brief  Called for each statement. The method writes the current statement to the file. This
       * overload is used for the Jacobian writers.
       */
      virtual void writeStatement(Identifier const& curLhsIdentifier, size_t& curJacobianPos,
                                  Real const* const rhsJacobians, Identifier const* const rhsIdentifiers,
                                  Config::ArgumentSize const& nJacobians) {
        CODI_UNUSED(curLhsIdentifier, curJacobianPos, rhsJacobians, rhsIdentifiers, nJacobians);
      }

      /**
       * @brief  Called for each statement. The method writes the current statement to the file. This
       * overload is used for the primal value writers and contains additional arguments, such as the WriteInfo.
       */
      virtual void writeStatement(WriteInfo const& info, Identifier const* lhsIdentifiers, Real const* lhsPrimalValues,
                                  Config::ArgumentSize const& nPassiveValues, Identifier const* const rhsIdentifiers,
                                  Real const* const passiveValues, Real const* const constantValues,
                                  EvalHandle stmtEvalHandle) {
        CODI_UNUSED(info, lhsIdentifiers, lhsPrimalValues, nPassiveValues, rhsIdentifiers, passiveValues,
                    constantValues, stmtEvalHandle);
      }

      /// Used for statements that contain a low level function.
      virtual void writeLowLevelFunction(LowLevelFunctionEntry<Tape, Real, Identifier> const* func,
                                         ByteDataView& data) {
        CODI_UNUSED(func, data);
      }

      /// After all the statements have been written, the finish method finalizes the writing process.
      virtual void finish() {}
  };

  /**
   * @brief  The interface is used by all the tape readers. The tape readers are used to restore a tape from either a
   * text file or a binary file that was created by the \ref codi::TapeWriterInterface. The readers produce an identical
   * tape in a new context, which can be evaluated and new statements can be added.
   *
   * The readers are supported for both Jacobian and primal value tapes, however reading a primal value tape requires
   * additional steps to restore the EvalHandles.
   *
   * For reading and restoring a Jacobian tape file, the following steps are followed.
   *
   * 1) Use the \ref codi::readTapeFile() with the name of the text or binary tape file that was previously created.
   * Ensure that the file with the "IO.txt" or "IO.dat" suffix is also available in the same directory.
   *
   * 2) The restored tape can then be accessed using the \ref codi::TapeReaderInterface::getTape() method.
   *
   * 3) Use the \ref codi::TapeReaderInterface::getInputs() or  \ref codi::TapeReaderInterface::getOutputs() to set
   * seed the gradients before evaluating the restored tape.
   *
   * 4) The same input and output vectors can be use to access the adjoint values after the evaluation.
   *
   *  \code{.cpp}
   * // Restore the tape from the text or binary file;
   * std::string fileName{"Example.txt"}
   * auto reader = codi::readTapeFile<Real>(fileName);
   *
   * // Use the reader to get the tape
   * Tape& tape = reader->getTape();
   *
   * // Seed the desired IO using the reader.
   * tape.gradient(reader->getOutputs()[0]) = 1;
   *
   * // Evaluate the tape
   * tape.evaluate();
   *
   * // Get the resulting gradients
   * std::cout << tape.gradient(textRead->getInputs()[0]) << std::endl;
   * \endcode
   *
   * The procedure for restoring a primal value tape is similar, however there is a header file <tt>"filename".hpp</tt>
   * that the codi::TapeWriterInterface produces, which needs to be added to reproduce the EvalHandles in the new
   * context.
   *
   * The example below shows the additional steps for a primal value tape
   *  \code{.cpp}
   *
   * // Include the header file that was produced by the primal value writer
   * #include "Example.hpp"
   *
   * // Call the createEvalHandles function that is contained in the header file
   * std::vector<typename Tape::EvalHandle> evalHandles = fileNameCreateEvalHandles<codi::Tape>();
   *
   * // Restore the tape from the text or binary file. The evalHandles are added to the function arguments.
   * std::string fileName{"Example.txt"}
   * auto reader = codi::readTapeFile<Real>(fileName, evalHandles);
   *
   * // Use the reader to get the tape
   * Tape& tape = reader->getTape();
   *
   * // Seed the desired IO using the reader.
   * tape.gradient(reader->getOutputs()[0]) = 1;
   *
   * // Evaluate the tape
   * tape.evaluate();
   *
   * // Get the resulting gradients
   * std::cout << tape.gradient(textRead->getInputs()[0]) << std::endl;
   * \endcode
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  struct TapeReaderInterface {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< The evaluation type.
      using Tape = typename Type::Tape;                           ///< The tape type that is to be restored.
      using Identifier = typename Type::Identifier;               ///< Identifier for the internal management, e.g. int.
      using Real = typename Type::Real;                           ///< Primal computation type, e.g. double.

      virtual ~TapeReaderInterface() {}  ///< Destructor

      /// This method uses the the fileName to reproduce a valid tape.
      virtual void readFile(std::string const& name) {
        CODI_UNUSED(name);
      }

      virtual Tape& getTape() = 0;  ///< Used to get a reference to the restored tape.

      virtual std::vector<Identifier>& getInputs() = 0;  ///< Used to get the restored inputs of the tape.

      virtual std::vector<Identifier>& getOutputs() = 0;  ///< Used to get the restored outputs of the
                                                          ///< tape.
  };
}
