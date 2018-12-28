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

#include "chunk.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Storage for pointer to chunk data.
   *
   * The pointer handle stores a pointer to each data item of the
   * chunk data, that is provided as a template.
   *
   * @tparam ChunkData  The chunk type needs implement the interface ChunkInterface.
   */
  template<typename ChunkData>
  struct PointerHandle {

      /**
       * @brief Set the internal pointers to the the data of the chunk at the
       * specific position.
       *
       * @param[in]     dataPos  The position of the data in the chunk that is set to the pointers.
       * @param[in,out]   chunk  The chunk from which the data is gained.
       */
      void setPointers(const size_t& dataPos, ChunkData* chunk);

      /**
       * @brief Call the function object with the arguments and the pointers from this handle object.
       *
       * The call is:
       * func(<pointers> , <args>);
       *
       * @param[in]     func  The function object which is called with the data.
       * @param[in,out] args  The additional arguments for the function call.
       */
      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args);

  };

  /**
   * @brief Specialization for PointerHandle with a Chunk1 type.
   *
   * @tparam Data1  The first data type for the chunk.
   */
  template<typename Data1>
  struct PointerHandle<Chunk1<Data1> > {
      Data1* p1; /**< Pointer for the first data item. */

      /**
       * @brief Set the internal pointers to the the data of the chunk at the
       * specific position.
       *
       * @param[in]     dataPos  The position of the data in the chunk that is set to the pointers.
       * @param[in,out]   chunk  The chunk from which the data is gained.
       */
      void setPointers(const size_t& dataPos, Chunk1<Data1>* chunk) {
        chunk->dataPointer(dataPos, p1);
      }

      /**
       * @brief Call the function object with the arguments and the pointers from this handle object.
       *
       * The call is:
       * func(<pointers> , <args>);
       *
       * @param[in]     func  The function object which is called with the data.
       * @param[in,out] args  The additional arguments for the function call.
       */
      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, std::forward<Args>(args)...);
      }

      /**
       * @brief Call reverse evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the reverse evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., p1);
      }

      /**
       * @brief Call forward evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the forward evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., p1);
      }
  };

  /**
   * @brief Specialization for PointerHandle with a Chunk2 type.
   *
   * @tparam Data1  The first data type for the chunk.
   * @tparam Data2  The second data type for the chunk.
   */
  template<typename Data1, typename Data2>
  struct PointerHandle<Chunk2<Data1, Data2> > {
      Data1* p1; /**< Pointer for the first data item. */
      Data2* p2; /**< Pointer for the second data item. */

      /**
       * @brief Set the internal pointers to the the data of the chunk at the
       * specific position.
       *
       * @param[in]     dataPos  The position of the data in the chunk that is set to the pointers.
       * @param[in,out]   chunk  The chunk from which the data is gained.
       */
      void setPointers(const size_t& dataPos, Chunk2<Data1, Data2>* chunk) {
        chunk->dataPointer(dataPos, p1, p2);
      }

      /**
       * @brief Call the function object with the arguments and the pointers from this handle object.
       *
       * The call is:
       * func(<pointers> , <args>);
       *
       * @param[in]     func  The function object which is called with the data.
       * @param[in,out] args  The additional arguments for the function call.
       */
      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, std::forward<Args>(args)...);
      }

      /**
       * @brief Call reverse evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the reverse evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., p1, p2);
      }

      /**
       * @brief Call forward evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the forward evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., p1, p2);
      }
  };

  /**
   * @brief Specialization for PointerHandle with a Chunk3 type.
   *
   * @tparam Data1  The first data type for the chunk.
   * @tparam Data2  The second data type for the chunk.
   * @tparam Data3  The third data type for the chunk.
   */
  template<typename Data1, typename Data2, typename Data3>
  struct PointerHandle<Chunk3<Data1, Data2, Data3> > {
      Data1* p1; /**< Pointer for the first data item. */
      Data2* p2; /**< Pointer for the second data item. */
      Data3* p3; /**< Pointer for the third data item. */

      /**
       * @brief Set the internal pointers to the the data of the chunk at the
       * specific position.
       *
       * @param[in]     dataPos  The position of the data in the chunk that is set to the pointers.
       * @param[in,out]   chunk  The chunk from which the data is gained.
       */
      void setPointers(const size_t& dataPos, Chunk3<Data1, Data2, Data3>* chunk) {
        chunk->dataPointer(dataPos, p1, p2, p3);
      }

      /**
       * @brief Call the function object with the arguments and the pointers from this handle object.
       *
       * The call is:
       * func(<pointers> , <args>);
       *
       * @param[in]     func  The function object which is called with the data.
       * @param[in,out] args  The additional arguments for the function call.
       */
      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, p3, std::forward<Args>(args)...);
      }

      /**
       * @brief Call reverse evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the reverse evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., p1, p2, p3);
      }

      /**
       * @brief Call forward evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the forward evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., p1, p2, p3);
      }
  };

  /**
   * @brief Specialization for PointerHandle with a Chunk4 type.
   *
   * @tparam Data1  The first data type for the chunk.
   * @tparam Data2  The second data type for the chunk.
   * @tparam Data3  The third data type for the chunk.
   * @tparam Data4  The fourth data type for the chunk.
   */
  template<typename Data1, typename Data2, typename Data3, typename Data4>
  struct PointerHandle<Chunk4<Data1, Data2, Data3, Data4> > {
      Data1* p1; /**< Pointer for the first data item. */
      Data2* p2; /**< Pointer for the second data item. */
      Data3* p3; /**< Pointer for the third data item. */
      Data4* p4; /**< Pointer for the fourth data item. */

      /**
       * @brief Set the internal pointers to the the data of the chunk at the
       * specific position.
       *
       * @param[in]     dataPos  The position of the data in the chunk that is set to the pointers.
       * @param[in,out]   chunk  The chunk from which the data is gained.
       */
      void setPointers(const size_t& dataPos, Chunk4<Data1, Data2, Data3, Data4>* chunk) {
        chunk->dataPointer(dataPos, p1, p2, p3, p4);
      }

      /**
       * @brief Call the function object with the arguments and the pointers from this handle object.
       *
       * The call is:
       * func(<pointers> , <args>);
       *
       * @param[in]     func  The function object which is called with the data.
       * @param[in,out] args  The additional arguments for the function call.
       */
      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, p3, p4, std::forward<Args>(args)...);
      }

      /**
       * @brief Call reverse evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the reverse evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., p1, p2, p3, p4);
      }

      /**
       * @brief Call forward evaluation on the nested vector with the pointers from this handle.
       *
       * @param[in,out] nested  The nested vector on which the forward evaluation is called.
       * @param[in,out]   args  The additional arguments for the reverse evaluation.
       */
      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., p1, p2, p3, p4);
      }
  };
}
