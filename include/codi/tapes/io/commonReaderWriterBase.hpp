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

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "../../config.h"
#include "tapeReaderWriterInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief Used to implement methods common to both the tape readers and the tape writers.
   */
  struct CommonReaderWriterMethods {
      std::string fileName;  ///< The base file name provided by the user.

      /// Constructor.
      CommonReaderWriterMethods(std::string const& name) : fileName(name) {};
      /// Constructor.
      CommonReaderWriterMethods() : fileName("") {};

      /// Remove the file extension and replace it with a new suffix.
      std::string modifyFileName(std::string const& suffix) {
        std::string::size_type sepPosition = fileName.rfind('.');
        return fileName.substr(0, sepPosition) + suffix;
      }

      /// Open a file and check for success. Failure terminates and prints an error.
      ///
      /// \c mode is open mode from \c fopen.
      void openFile(FILE*& fileHandle, std::string const& name, std::string const& mode) {
        fileHandle = fopen(name.c_str(), mode.c_str());
        if (nullptr == fileHandle) {
          CODI_EXCEPTION("Could not open file %s", name.c_str());
        }
      }
  };

  /**
   * @brief  This class is a common base for all the writers and produces a IO file that contains the input and output
   * variables of the current tape. The IO file is written in a binary or text format.
   *
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct CommonBaseTapeWriter : public TapeWriterInterface<T_Type>, public CommonReaderWriterMethods {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeWriterInterface.
      using Tape = typename Type::Tape;                           ///< See TapeWriterInterface.

      std::vector<Identifier> inputVariables;   ///< The identifiers which have been registered as inputs.
      std::vector<Identifier> outputVariables;  ///< The identifiers which have been registered as outputs.

      /// Constructor.
      CommonBaseTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                           std::vector<Identifier> const& out)
          : CommonReaderWriterMethods(name), inputVariables(in), outputVariables(out) {};

    protected:

      /// Generate the IO file in a text format.
      void printIoText(Tape& tape) {
        FILE* ioFileHandleTxt = nullptr;
        this->openFile(ioFileHandleTxt, modifyFileName("IO.txt"), "w");

        fprintf(ioFileHandleTxt, "%zu Inputs = ", inputVariables.size());
        for (size_t inputCount = 0; inputCount < inputVariables.size(); inputCount++) {
          fprintf(ioFileHandleTxt, "%d   ", inputVariables[inputCount]);
        }
        fprintf(ioFileHandleTxt, "\n%zu Outputs = ", outputVariables.size());
        for (size_t outputCount = 0; outputCount < outputVariables.size(); outputCount++) {
          fprintf(ioFileHandleTxt, "%d   ", outputVariables[outputCount]);
        }
        fprintf(ioFileHandleTxt, "\nLargest Index = %d", tape.getIndexManager().getLargestCreatedIndex());
        fclose(ioFileHandleTxt);
      }

      /// Generate the IO file in a binary format.
      void printIoBinary(Tape& tape) {
        FILE* ioFileHandleBinary = nullptr;

        this->openFile(ioFileHandleBinary, modifyFileName("IO.dat"), "wb");

        size_t nInputs = inputVariables.size();
        size_t nOutputs = outputVariables.size();

        fwrite(&nInputs, sizeof(nInputs), 1, ioFileHandleBinary);
        fwrite(&inputVariables[0], sizeof(Identifier), nInputs, ioFileHandleBinary);
        fwrite(&nOutputs, sizeof(nOutputs), 1, ioFileHandleBinary);
        fwrite(&outputVariables[0], sizeof(Identifier), nOutputs, ioFileHandleBinary);
        Identifier largestIndex = tape.getIndexManager().getLargestCreatedIndex();
        fwrite(&largestIndex, sizeof(Identifier), 1, ioFileHandleBinary);

        fclose(ioFileHandleBinary);
      }
  };

  /**
   * @brief  Used to restore the IO from the <tt>"filename"IO.dat</tt> or <tt>"filename"IO.txt</tt> files. Also provides
   * the get methods from the codi::tapeReaderInterface.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  struct CommonBaseTapeReader : public TapeReaderInterface<T_Type>, public CommonReaderWriterMethods {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeReaderInterface.
      using Tape = typename Type::Tape;                           ///< See TapeReaderInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeReaderInterface.

      std::vector<Identifier> inputVariables;   ///< Stores the restored input variables from the tape file.
      std::vector<Identifier> outputVariables;  ///< Stores the restored output variables from the tape file.
      size_t nInputs;                           ///< Size of the input vector.
      size_t nOutputs;                          ///< Size of the output vector.

      Tape tape;  ///< The newly resorted tape.

      Identifier largestIndex;                                          ///< The largest index on the stored tape.
      CommonBaseTapeReader() : CommonReaderWriterMethods(), tape() {};  ///< Constructor.

    protected:

      /// Restore the IO for the text readers.
      void restoreIoText() {
        FILE* ioFileHandleTxt;
        this->openFile(ioFileHandleTxt, modifyFileName("IO.txt"), "r+");

        fscanf(ioFileHandleTxt, "%zu Inputs = ", &nInputs);
        inputVariables.resize(nInputs, 0);
        for (size_t inputCount = 0; inputCount < nInputs; inputCount++) {
          fscanf(ioFileHandleTxt, "%d   ", &inputVariables[inputCount]);
        }
        fscanf(ioFileHandleTxt, "\n%zu Outputs = ", &nOutputs);
        outputVariables.resize(nOutputs, 0);
        for (size_t outputCount = 0; outputCount < nOutputs; outputCount++) {
          fscanf(ioFileHandleTxt, "%d   ", &outputVariables[outputCount]);
        }
        fscanf(ioFileHandleTxt, "\nLargest Index = %d", &largestIndex);
        fclose(ioFileHandleTxt);
      }

      /// Restore the IO for the binary readers.
      void restoreIoBinary() {
        FILE* fileIoHandleBin;
        this->openFile(fileIoHandleBin, modifyFileName("IO.dat"), "rb");

        fread(&nInputs, sizeof(nInputs), 1, fileIoHandleBin);
        inputVariables.resize(nInputs, 0);
        fread(inputVariables.data(), sizeof(Identifier), nInputs, fileIoHandleBin);
        fread(&nOutputs, sizeof(nOutputs), 1, fileIoHandleBin);
        outputVariables.resize(nOutputs, 0);
        fread(outputVariables.data(), sizeof(Identifier), nOutputs, fileIoHandleBin);
        fread(&largestIndex, sizeof(Identifier), 1, fileIoHandleBin);
      }

      /// \copydoc codi::TapeReaderInterface::getTape()
      Tape& getTape() {
        return tape;
      }

      /// \copydoc codi::TapeReaderInterface::getInputs()
      std::vector<Identifier> const& getInputs() const& {
        return inputVariables;
      }

      /// \copydoc codi::TapeReaderInterface::getOutputs()
      std::vector<Identifier> const& getOutputs() const& {
        return outputVariables;
      }

      /// This method is used to remove any offset and to update the largest created index.
      void updateUserIO(Identifier const& linearOffset) {
        if (Tape::TapeTypes::IsLinearIndexHandler && !TapeTraits::IsPrimalValueTape<Tape>::value) {
          // For the linear case, user IO is updated with a potential offset.
          for (size_t inputCount = 0; inputCount < inputVariables.size(); inputCount++) {
            inputVariables[inputCount] -= linearOffset;
          }

          for (size_t outputCount = 0; outputCount < outputVariables.size(); outputCount++) {
            outputVariables[outputCount] -= linearOffset;
          }
        }
      }
  };

  /**
   * @brief This base class is used to modify the math representation of a statement.
   *
   * Additionally, the class is used by the graphical writers to create nodes and edges.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct CommonTextTapeWriter : public CommonBaseTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Base = CommonBaseTapeWriter<T_Type>;                  ///< See CommonBaseTapeWriter.

      using Tape = typename Type::Tape;              ///< See TapeWriterInterface.
      using Identifier = typename Type::Identifier;  ///< See TapeWriterInterface.
      using Real = typename Type::Real;              ///< See TapeWriterInterface.

      bool writeDotHeaderFooter;  ///< This flag sets toggles the heading and color index.

      FILE* fileHandleGraph = nullptr;  ///< The handle for the writer.

      std::vector<Identifier> identifierExtensions;  ///< Used to differentiate multiple instances of the same
                                                     ///< identifier.
      std::vector<IdentifierType> identifierType;    ///< Used to differentiate between an input, output and temp
                                                     ///< variable.

      CommonTextTapeWriter(bool writeDotHeaderFooter, std::string const& name, std::vector<Identifier> const& in,
                           std::vector<Identifier> const& out)
          : Base(name, in, out),
            writeDotHeaderFooter(writeDotHeaderFooter),
            fileHandleGraph(nullptr),
            identifierExtensions(0),
            identifierType(0) {}  ///< Constructor.

      /// \copydoc codi::TapeWriterInterface::start()
      void start(Tape& tape) {
        // Resize vectors to maximum index.
        identifierExtensions.resize(tape.getIndexManager().getLargestCreatedIndex() + 1, 0);
        identifierType.resize(tape.getIndexManager().getLargestCreatedIndex() + 1, IdentifierType::Temp);

        // Record the input and output identifiers in the identifierType vector. This is used later to colour code the
        // nodes and avoids searching through the input and output vectors for each statement.
        for (size_t inputCount = 0; inputCount < this->inputVariables.size(); inputCount++) {
          identifierType[this->inputVariables[inputCount]] = IdentifierType::Input;
        }

        for (size_t outputCount = 0; outputCount < this->outputVariables.size(); outputCount++) {
          identifierType[this->outputVariables[outputCount]] = IdentifierType::Output;
        }

        this->openFile(fileHandleGraph, this->fileName, "w");

        if (writeDotHeaderFooter) {
          // Print out the header and add a color index.
          fprintf(fileHandleGraph,
                  "digraph Net {\nInputs [label = \"Inputs\", color=\"blue\"];\nOutputs [label = \"Outputs\", "
                  "color=\"red\"];\nInter [label = \"Inter\"];\n");
        }
      }

      /// \copydoc codi::TapeWriterInterface::finish()
      void finish() {
        if (writeDotHeaderFooter) {
          fprintf(this->fileHandleGraph, "}");
        }
        fclose(this->fileHandleGraph);
      }

    protected:
      /// Add the identifier extension of the identifier to the node name.
      std::string formatNodeName(Identifier const& identifier, int extensionOffset = 0) {
        return "A" + std::to_string(identifier) + "_" +
               std::to_string(identifierExtensions[identifier] + extensionOffset);
      }

      /// Returns the color for a given identifier.
      std::string nodeColorProperties(Identifier const& identifier) {
        if (identifierType[identifier] == IdentifierType::Input) {
          return "blue";
        } else if (identifierType[identifier] == IdentifierType::Output) {
          return "red";
        } else {
          return "black";
        }
      }

      /// Creates a new node for a given Identifier and label.
      void createNode(Identifier const& identifier, std::string const& label) {
        std::string node = formatNodeName(identifier, 1) + " [label = \"" + label + "\", color=\"" +
                           nodeColorProperties(identifier) + "\"];\n";
        fprintf(this->fileHandleGraph, "%s", node.c_str());
      }

      /// Return a string with the current identifier type and the identifier value.
      std::string formatNodeLabel(Identifier const& identifier) {
        std::string result;
        if (this->identifierType[identifier] == IdentifierType::Input) {
          result = "X";
        } else if (this->identifierType[identifier] == IdentifierType::Output) {
          result = "Y";
        } else {
          result = "T";
        }

        result += std::to_string(identifier);

        return result;
      }

      /// Replaces all substrings of \c search with \c replace in \c str.
      void replaceAll(std::string& str, const std::string& search, const std::string& replace) {
        size_t start_pos = 0;
        while ((start_pos = str.find(search, start_pos)) != std::string::npos) {
          str.replace(start_pos, search.length(), replace);
          start_pos += replace.length();
        }
      }

      /**
       * @brief Replaces all general identifiers in the math representation with the input, output or temporary
       * annotation. E.g. <tt>x42 -> Y42</tt> for a variable tagged as output.
       */
      std::string modifyMathRep(std::string const& mathRep, Identifier const& lhsIdentifier,
                                Identifier const* const rhsIdentifiers, size_t const& nActiveValues) {
        std::string result = formatNodeLabel(lhsIdentifier) + " = " + mathRep;
        // Iterate through the RHS and replace x1..xn with the identifier type and the identifier value.
        for (size_t curArg = 0; curArg < nActiveValues; curArg++) {
          std::string searchString = "x" + std::to_string(rhsIdentifiers[curArg]);
          std::string replaceString = formatNodeLabel(rhsIdentifiers[curArg]);

          replaceAll(result, searchString, replaceString);
        }

        return result;
      }

      /// Ensure that all the nodes on the rhs have been placed in the .dot file before creating edges to them.
      void placeUnusedRhsNodes(Identifier const* const rhsIdentifiers, Config::ArgumentSize const& nArguments) {
        // Check if the identifierExtension is zero for any of the rhsIdentifiers. A zero extension indicates that the
        // node has not been placed. The type of the identifier is then checked to add the correct colour coding.
        for (size_t argCount = 0; argCount < nArguments; argCount++) {
          if (this->identifierExtensions[rhsIdentifiers[argCount]] == 0) {
            createNode(rhsIdentifiers[argCount], formatNodeLabel(rhsIdentifiers[argCount]));
            // Increment the extension of the newly placed identifier.
            this->identifierExtensions[rhsIdentifiers[argCount]] += 1;
          }
        }
      }

      /// Return a string that creates an edge between two nodes in the .dot language.
      void createEdge(Identifier const& from, Identifier const& to, std::string const& label = "") {
        std::string edge = formatNodeName(from) + " -> " + formatNodeName(to, 1);

        if (0 != label.size()) {
          edge += " [label=\"" + label + "\"];";
        }
        edge += "\n";
        fprintf(this->fileHandleGraph, "%s", edge.c_str());
      }
  };
}
