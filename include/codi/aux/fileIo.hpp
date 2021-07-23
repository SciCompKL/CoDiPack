
#pragma once

#include <errno.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <string>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Possible IO errors.
  enum struct IoError
  {
    Mode,
    Open,
    Write,
    Read
  };

  /// IoException for CoDiPack.
  struct IoException {
    public:

      std::string text;  ///< Textual description.
      IoError id;        ///< Exception ID.

      /// Constructor
      IoException(IoError id, std::string const& text, bool appendErrno) : text(text), id(id) {
        if (appendErrno) {
          this->text += " (Internal error: ";
          this->text += strerror(errno);
          this->text += ")";
        }
      }
  };

  /**
   * @brief Helper structure for writing binary data.
   *
   * Exceptions are thrown if:
   *  - file could not be opened,
   *  - file is used in the wrong mode,
   *  - number of bytes read/written is wrong.
   */
  struct FileIo {
    private:

      FILE* fileHandle;  ///< File handle
      bool writeMode;    ///< true = write, false = read

    public:

      /// Constructor
      /// Will throw an IoException if file cannot be opened.
      FileIo(std::string const& file, bool write) {
        writeMode = write;
        fileHandle = nullptr;

        if (write) {
          fileHandle = fopen(file.c_str(), "wb");
        } else {
          fileHandle = fopen(file.c_str(), "rb");
        }

        if (nullptr == fileHandle) {
          throw IoException(IoError::Open, "Could not open file: " + file, true);
        }
      }

      /// Destructor
      ~FileIo() {
        if (nullptr != fileHandle) {
          fclose(fileHandle);
        }
      }

      /// Write data to a file.
      /// Will throw an IoException if not in write mode or if the number of bytes written is wrong.
      template<typename Data>
      void writeData(Data const* data, size_t const length) {
        if (writeMode) {
          size_t s = fwrite(data, sizeof(Data), length, fileHandle);

          if (s != length) {
            throw IoException(IoError::Read, "Wrong number of bytes written.", true);
          }
        } else {
          throw IoException(IoError::Mode, "Using write io handle in wrong mode.", false);
        }
      }

      /// Read data from a file.
      /// Will throw an IoException if not in read mode or if the number of bytes read is wrong.
      template<typename Data>
      void readData(Data* data, size_t const length) {
        if (!writeMode) {
          size_t s = fread(data, sizeof(Data), length, fileHandle);

          if (s != length) {
            throw IoException(IoError::Read, "Wrong number of bytes read.", false);
          }
        } else {
          throw IoException(IoError::Mode, "Using read io handle in wrong mode.", false);
        }
      }
  };
}
