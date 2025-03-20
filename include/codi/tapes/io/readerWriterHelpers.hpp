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

#include <map>

#include "../../config.h"
#include "graphWriters.hpp"
#include "jacobianBinaryReaderWriter.hpp"
#include "jacobianTextReaderWriter.hpp"
#include "mathRepWriter.hpp"
#include "primalBinaryReaderWriter.hpp"
#include "primalTextReaderWriter.hpp"
#include "tapeReaderWriterInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper for creating a unique pointer. (std::make_unique is not available in C++11.)
  template<typename T, typename... Args>
  std::unique_ptr<T> make_unique_helper(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  /**
   * @brief  The createWriter() function is used to generate an automatic writer using the FileType and the TapeTraits.
   *
   * The method uses default values for the writers with parameters, such as printJacobians for the graphWriter.
   *
   * Smart Pointers are used to avoid manual memory management.
   *
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  std::unique_ptr<TapeWriterInterface<T_Type>> createWriter(const std::string& fileName,
                                                            std::vector<typename T_Type::Identifier>& inputVariables,
                                                            std::vector<typename T_Type::Identifier>& outputVariables,
                                                            FileType selectedType) {
    using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< The evaluation type.
    using Tape = typename Type::Tape;                           ///< The tape type that is to be written out.

    std::unique_ptr<TapeWriterInterface<T_Type>> writer = nullptr;
    static constexpr bool IsPrimalValueTape = TapeTraits::IsPrimalValueTape<Tape>::value;

    switch (selectedType) {
      case FileType::Text:
        if (IsPrimalValueTape) {
          writer = make_unique_helper<codi::PrimalTextTapeWriter<Type>>(fileName, inputVariables, outputVariables, true,
                                                                        true);
        } else {
          writer = make_unique_helper<codi::JacobianTextTapeWriter<Type>>(fileName, inputVariables, outputVariables,
                                                                          true, true);
        }
        break;

      case FileType::Binary:
        if (IsPrimalValueTape) {
          writer = make_unique_helper<codi::PrimalBinaryTapeWriter<Type>>(fileName, inputVariables, outputVariables);
        } else {
          writer = make_unique_helper<codi::JacobianBinaryTapeWriter<Type>>(fileName, inputVariables, outputVariables);
        }
        break;

      case FileType::Graph:
        if (IsPrimalValueTape) {
          writer = make_unique_helper<codi::PrimalGraphTapeWriter<Type>>(fileName, inputVariables, outputVariables);
        } else {
          writer =
              make_unique_helper<codi::JacobianGraphTapeWriter<Type>>(fileName, inputVariables, outputVariables, true);
        }
        break;

      case FileType::Math:
        if (IsPrimalValueTape) {
          writer = make_unique_helper<codi::MathRepWriter<Type>>(fileName, inputVariables, outputVariables);
        } else {
          CODI_EXCEPTION("The MathRepWriter is not supported for Jacobian tapes.");
        }
        break;

      default:
        CODI_EXCEPTION("A valid file format was not selected.");
        break;
    }

    return writer;
  }

  /**
   * @brief  Uses the fileName to read and restore a Jacobian tape. The file extension is used to determine wether the
   * tape file is a binary (.dat) or text (.txt) file.
   *
   * This overload is specifically for a Jacobian tape, where the EvalHandles are not required.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  std::unique_ptr<TapeReaderInterface<T_Type>> readTapeFile(std::string const& fileName) {
    // This function uses the file extension provided in the fileName to determine the type of TapeReader.

    /// The reader that will contain the read IO and restored tape.
    std::unique_ptr<TapeReaderInterface<T_Type>> reader = nullptr;
    std::string fileExtension;  ///< The extension is used to determine the file type.

    size_t dotPosition = fileName.find_last_of(".");
    if (dotPosition != std::string::npos) {
      fileExtension = fileName.substr(dotPosition + 1);
    }

    if (fileExtension == "txt") {
      reader = make_unique_helper<JacobianTextTapeReader<T_Type>>();
    } else if (fileExtension == "dat") {
      reader = make_unique_helper<JacobianBinaryTapeReader<T_Type>>();
    } else {
      CODI_EXCEPTION("The given file type is not supported.");
    }
    reader->readFile(fileName);

    return reader;
  }

  /**
   * @brief  Uses the fileName to read and restore a primal value tape. The file extension is used to determine whether
   * the tape file is a binary (.dat) or text (.txt) file.
   *
   * This overload is specifically for a primal value tape, which requires the EvalHandles as an argument. They can be
   * created with a call to \c createEvalHandles from the generated header file. See TapeReaderInterface for an example.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  std::unique_ptr<TapeReaderInterface<T_Type>> readTapeFile(
      std::string const& fileName, std::vector<typename T_Type::Tape::EvalHandle> const& evalHandles) {
    // This function uses the file extension provided in the fileName to determine the type of TapeReader.

    /// The reader that will contain the read IO and restored tape.
    std::unique_ptr<TapeReaderInterface<T_Type>> reader = nullptr;
    std::string fileExtension;  ///< The extension is used to determine the file type.

    size_t dotPosition = fileName.find_last_of(".");
    if (dotPosition != std::string::npos) {
      fileExtension = fileName.substr(dotPosition + 1);
    }

    if (fileExtension == "txt") {
      reader = make_unique_helper<PrimalTextTapeReader<T_Type>>(evalHandles);
    } else if (fileExtension == "dat") {
      reader = make_unique_helper<PrimalBinaryTapeReader<T_Type>>(evalHandles);
    } else {
      CODI_EXCEPTION("The given file type is not supported.");
    }

    reader->readFile(fileName);

    return reader;
  }
}
