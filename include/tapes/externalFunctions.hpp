/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: TODO
 */

#pragma once

namespace codi {

  struct ExternalFunction {
    typedef void (*CallFunction)(void*);
    typedef void (*DeleteFunction)(void*);

  private:
    CallFunction func;
    DeleteFunction deleteCheckpoint;

    void* checkpoint;

  public:
    ExternalFunction(){}
    ExternalFunction(CallFunction func, void* checkpoint, DeleteFunction deleteCheckpoint) :
      func(func),
      deleteCheckpoint(deleteCheckpoint),
      checkpoint(checkpoint){}

    void deleteData() {
      if (deleteCheckpoint != NULL){
        deleteCheckpoint(checkpoint);
        checkpoint = NULL;
      }
    }

    void evaluate() {
      if(NULL != func) {
        func(checkpoint);
      }
    }
  };

  template<typename Data>
  class ExternalFunctionDataHelper {
  public:
    typedef void (*CallFunction)(Data*);
    typedef void (*DeleteFunction)(Data*);

  private:
    CallFunction func;
    DeleteFunction deleteData;

    Data* data;

  public:
    ExternalFunctionDataHelper(CallFunction func, Data* data, DeleteFunction deleteData) :
      func(func),
      deleteData(deleteData),
      data(data){}


    static void callFunction(void* data) {
      ExternalFunctionDataHelper<Data>* castData = cast(data);
      castData->func(castData->data);
    }

    static void deleteFunction(void* data) {
      ExternalFunctionDataHelper<Data>* castData = cast(data);
      castData->deleteData(castData->data);

      // delete self
      delete castData;
    }

  private:

    static ExternalFunctionDataHelper<Data>* cast(void* data) {
      return (ExternalFunctionDataHelper<Data>*)data;
    }

  };
}
