/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

/*
 * In order to include this file the user has to define the preprocessor macro ROOT_VECTOR.
 *
 * ROOT_VECTOR defines the name of the root vector for all the data of the tape.
 *
 * All these macros are undefined at the end of the file.
 *
 * It defines the methods writeToFile, readFromFile, deleteData, resetHard.
 */

#ifndef ROOT_VECTOR
  #error Please define the name of the root vector.
#endif

  private:

  // ----------------------------------------------------------------------
  // Private functions of the module
  // ----------------------------------------------------------------------

    /**
     * @brief Helper function, that writes a chunk to the io handle
     *
     * @param[in]      chunk  The chunk which is written by the io handle.
     * @param[in,out] handle  The io handle
     */
    static void writeFunction(const ChunkInterface* chunk, CoDiIoHandle& handle) {
      chunk->writeData(handle);
    }

    /**
     * @brief Helper function, that reads a chunk from the io handle
     *
     * @param[in,out]  chunk  The chunk which is read by the io handle.
     * @param[in,out] handle  The io handle
     */
    static void readFunction(ChunkInterface* chunk, CoDiIoHandle& handle) {
      chunk->readData(handle);
    }

    /**
     * @brief Helper function, that removes the data of all the chunks without deleting the position information.
     */
    static void deleteFunction(ChunkInterface* chunk) {
      chunk->deleteData();
    }

  public:

  // ----------------------------------------------------------------------
  // Public functions for the tape interface
  // ----------------------------------------------------------------------

    /**
     * @brief Write a binary blob of the whole tape data.
     *
     * Only the data of the chunks is written to the file.
     * The position information is for example not written.
     *
     * The data file can only be read by this tape and can not be used
     * by another tape.
     *
     * @param[in] filename  The name of the file.
     */
    void writeToFile(const std::string& filename) {
      CoDiIoHandle ioHandle(filename, true);

      // we ignore the external function vector here because the data there should not be written
      ROOT_VECTOR.forEachChunk(writeFunction, true, ioHandle);
    }

    /**
     * @brief Read a binary blob of the whole tape data.
     *
     * Only the data of the chunks is read from the file.
     * The position information is for example not read.
     *
     * See also writeToFile
     *
     * @param[in] filename  The name of the file.
     */
    void readFromFile(const std::string& filename) {
      CoDiIoHandle ioHandle(filename, false);

      ROOT_VECTOR.forEachChunk(readFunction, true, ioHandle);
    }

    /**
     * @brief Delete all the data of the chunks such that the data is released.
     *
     * This will leave the structure in an invalid state. Only after a call
     * to readFromFile the state of the structure is valid again.
     *
     */
    void deleteData() {
      ROOT_VECTOR.forEachChunk(deleteFunction, true);
    }

    /**
     * @brief Reset the position of the tape to the zero position and release all acquired data.
     *
     * For chunk tapes all chunks are released.
     * For unchecked tapes, the size is set to zero.
     */
    void resetHard() {

      reset();

      cleanTapeBase();
      ROOT_VECTOR.resetHard();
    }

#undef ROOT_VECTOR
