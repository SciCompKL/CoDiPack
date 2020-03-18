
#pragma once

#include <iostream>
#include <stdio.h>
#include <errno.h>
#include <string>
#include <string.h>

#include "../config.h"
#include "macros.h"

/** \copydoc codi::Namespace */
namespace codi {

  enum struct IoError {
    Mode,
    Open,
    Write,
    Read
  };

  struct IoException {

      std::string text;
      IoError id;

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

  struct FileIo {

      FILE* fileHandle;
      bool writeMode;

    public:

      FileIo(const std::string& file, bool write) {
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
