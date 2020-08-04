
#pragma once

#include <iostream>
#include <stdio.h>
#include <errno.h>
#include <string>
#include <string.h>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  enum struct IoError {
    Mode,
    Open,
    Write,
    Read
  };

  struct IoException {
    public:

      std::string text;
      IoError id;

      IoException(IoError id, std::string const& text, bool appendErrno) :
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

  struct FileIo {
    private:

      FILE* fileHandle;
      bool writeMode;

    public:

      FileIo(std::string const& file, bool write) {
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

      ~FileIo() {
        if(NULL != fileHandle) {
          fclose(fileHandle);
        }
      }

      template<typename Data>
      void writeData(Data const* data, size_t const length) {
        if(writeMode) {
          size_t s = fwrite(data, sizeof(Data), length, fileHandle);

          if(s != length) {
            throw IoException(IoError::Read, "Wrong number of bytes written.", true);
          }
        } else {
          throw IoException(IoError::Mode, "Using write io handle in wrong mode.", false);
        }
      }

      template<typename Data>
      void readData(Data* data, size_t const length) {
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
