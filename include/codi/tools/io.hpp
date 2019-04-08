/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#pragma once

#include <iostream>
#include <stdio.h>
#include <errno.h>
#include <string>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Types of CoDiPack io errors.
   *
   * Mode: The wrong mode was used on the file. e.g. Write on a file opened for read.
   * Open: File could not be opened.
   * Write: Error during the write of some data. e.g. No space left.
   * Read: Error during the read of some data. e.g. eof reached.
   */
  enum struct IoError {
    Mode,
    Open,
    Write,
    Read
  };

  /**
   * @brief Common exception for all io errors.
   *
   * This exception needs to be caught if io operation are performed on any tape.
   */
  struct IoException {

      /** @brief The text of the error. */
      std::string text;

      /** @brief The id of the error. */
      IoError id;

      /**
       * @brief Create a new exception from the given text.
       *
       * The message of the error is text. If appendErrno is true, then the
       * error message from errno is converted to a string and appended.
       *
       * @param[in]          id  The id of the error.
       * @param[in]        text  The text of the error message.
       * @param[in] appendErrno  Flag if the error from errno should be appended to the text.
       */
      IoException(IoError id, const std::string& text, bool appendErrno) :
        text(text),
        id(id)
      {
        if(appendErrno) {
          this->text += " (Internal error: ";
          this->text += strerror(errno);
          this->text += ")";
        }
      }
  };

  /**
   * @brief Contains methods for writing and reading data to and from files.
   *
   * The handle provides a save way to open a file and write or read from that file.
   * The file is opened in binary mode.
   */
  class CoDiIoHandle {

      /** @brief The handle for the file. */
      FILE* fileHandle;

      /** @brief The write mode of the file. Used for error checking. */
      bool writeMode;

    public:

      /**
       * @brief Create a handle from the given file and with the specified mode
       *
       * The file is opened in binary mode.
       *
       * @param[in] file  The name of the file.
       * @param[in] write  If the file is opened for reading. Otherwise for writing.
       */
      CoDiIoHandle(const std::string& file, bool write) {
        writeMode = write;
        fileHandle = NULL;

        if(write) {
          fileHandle = fopen(file.c_str(), "wb");
        } else {
          fileHandle = fopen(file.c_str(), "rb");
        }

        if(NULL == fileHandle) {
          throw IoException(IoError::Open , "Could not open file: " + file, true);
        }
      }

      /**
       * @brief Close the file.
       */
      ~CoDiIoHandle() {
        if(NULL != fileHandle) {
          fclose(fileHandle);
        }
      }

      /**
       * @brief Write a blob of data to the file.
       *
       * @param[in]   data  The data that is writen to the file.
       * @param[in] length  The number of items of the array.
       *
       * @tparam Data  The type of the data items.
       */
      template<typename Data>
      void writeData(const Data* data, const size_t length) {
        if(writeMode) {
          size_t s = fwrite(data, sizeof(Data), length, fileHandle);

          if(s != length) {
            throw IoException(IoError::Read, "Wrong number of bytes written.", true);
          }
        } else {
          throw IoException(IoError::Mode, "Using write io handle in wrong mode.", false);
        }
      }

      /**
       * @brief Read a blob of data from the file.
       *
       * @param[out]  data  The data that is read from the file.
       * @param[in] length  The number of items of the array.
       *
       * @tparam Data  The type of the data items.
       */
      template<typename Data>
      void readData(Data* data, const size_t length) {
        if(!writeMode) {
          size_t s = fread(data, sizeof(Data), length, fileHandle);

          if(s != length) {
            throw IoException(IoError::Read, "Wrong number of bytes read.", false);
          }
        } else {
          throw IoException(IoError::Mode, "Using read io handle in wrong mode.", false);
        }
      }
  };
}
