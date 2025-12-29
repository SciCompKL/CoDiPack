/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
#include "../../misc/demangleName.hpp"
#include "commonReaderWriterBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief  This base class is used to implement the automatic generation of the .hpp file that restores the evaluation
   * handles in the reading process.
   *
   * See codi::TapeWriterInterface for a description on how to use the primal tape writers with the generated header
   * file.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct PrimalBaseTapeWriter : public CommonBaseTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Tape = typename Type::Tape;                           ///< See TapeWriterInterface.

      using Identifier = typename Type::Identifier;  ///< See TapeWriterInterface.
      using Real = typename Type::Real;              ///< See TapeWriterInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See TapeWriterInterface.

      std::map<EvalHandle, size_t> existingEvalHandles;  ///< Contains the existing handles and the assigned index.
      size_t evalHandleCount;                            ///< Count the number of unique evalHandles.
      std::vector<std::string> evalHandleStatements;     ///< The unique evalHandleStatements.

      /// Constructor.
      PrimalBaseTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                           std::vector<Identifier> const& out)
          : CommonBaseTapeWriter<T_Type>(name, in, out),
            existingEvalHandles(),
            evalHandleCount(0),
            evalHandleStatements() {};

    protected:

      /// This method is used to generate an .hpp file which creates the necessary EvalHandles in the reading process.
      void generateHandleCreatorFile() {
        FILE* fileHandleCreator = nullptr;
        std::string handleCreatorFileName = this->modifyFileName(".hpp");
        this->openFile(fileHandleCreator, handleCreatorFileName, "w");

        std::string functionName = this->modifyFileName("CreateEvalHandles");
        std::string fileHeader =
            "#include <codi.hpp>\n\n"
            "template <typename Tape>\n"
            "std::vector<typename Tape::EvalHandle> " +
            functionName +
            "(){\n\n"
            "  std::vector<typename Tape::EvalHandle> evalHandles;\n"
            "  using Impl = ";

        // Add the Impl type to the fileHeader
        fileHeader += demangleName<Tape>() + ";\n\n";
        // Resize the evalHandles
        fileHeader += "  evalHandles.resize(" + std::to_string(evalHandleCount) + ");\n";

        // Stmt strings
        std::string frontOfStmt = "  evalHandles[";
        std::string backOfStmt = "] = Tape::StatementEvaluator::template createHandle<";

        // Create fileFooter
        std::string fileFooter = "return evalHandles;\n}";

        // Header print
        fprintf(fileHandleCreator, "%s", fileHeader.c_str());

        // Use CreateHandle for all the statements stored in the lookup table.
        for (size_t handleCounter = 0; handleCounter < evalHandleStatements.size(); handleCounter++) {
          fprintf(fileHandleCreator, "%s%zu%s%s>();\n", frontOfStmt.c_str(), handleCounter, backOfStmt.c_str(),
                  evalHandleStatements[handleCounter].c_str());
        }
        // Footer print
        fprintf(fileHandleCreator, "\n  %s", fileFooter.c_str());
      }

      /// Get the index for an evalHandle.
      size_t getEvalHandleIndex(EvalHandle const evalHandle, std::string const& evalStatement) {
        auto it = existingEvalHandles.find(evalHandle);
        if (it != existingEvalHandles.end()) {
          return it->second;
        } else {
          existingEvalHandles[evalHandle] = evalHandleCount;
          evalHandleStatements.push_back(evalStatement);
          evalHandleCount += 1;
          return evalHandleCount - 1;
        }
      }
  };

}
