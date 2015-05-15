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
      ExternalFunction(){}
      ExternalFunction(void(*func_)(void*),
                   void* checkpoint_,void (*delete_checkpoint_)(void*))
          :func(func_), checkpoint(checkpoint_) ,
           delete_checkpoint(delete_checkpoint_){}

      void (*func)(void*);

      void* checkpoint;

      void (*delete_checkpoint)(void*);

      ~ExternalFunction(){
        if (delete_checkpoint != NULL){
          delete_checkpoint(checkpoint);
        }
      }
  };
}
