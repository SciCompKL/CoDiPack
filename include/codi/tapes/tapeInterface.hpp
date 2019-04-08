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

#include "../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Interface common to all tapes.
   *
   * The basic interface each tape has to implement. It defines functions which are used
   * by the active type to signal the tape when an active type is created or destroyed. For
   * each operation the active type calls also the store operation to inform the tape that
   * an expression is assigned to an active type.
   *
   * @tparam             Real  Floating point type of the gradients.
   * @tparam GradientDataType  The data the tape uses to identify each active variable
   *                             and where the tape can store information about the
   *                             gradient.
   * @tparam GradientValueType The value type that is used for the gradient calculation.
   */
  template <typename Real, typename GradientDataType, typename GradientValueType>
  class TapeInterface {
  public:

    /**
     * @brief The data for the gradient information of the tape.
     *
     * Each tape can define data for the gradient which each active type will
     * define in its own class. The tape can use this data to identify each active
     * type and compute the gradient information.
     */
    typedef GradientDataType GradientData;

    /**
     * @brief The actual values for the gradient information of the tape.
     *
     * The actual floating point or vector valued data for the gradient information.
     */
    typedef GradientValueType GradientValue;

    /**
     * This functions are called from the expression templates. They tell the
     * tape about the operations which is evaluated.
     */

    /**
     * @brief Set the value of the expression and evaluate the gradient.
     *
     * This function is called for every operation which resembles
     *
     *   lhs = rhs;
     *
     * Store sets the primal value of the operation which is stored and the tape can
     * do the evaluation of the gradient and perform other appropriate steps.
     *
     * @param[out]          lhsValue    The primal value of the lhs. This value is set to the value
     *                                  of the right hand side.
     * @param[out]   lhsGradientData    The gradient data of the lhs. The tape will update the gradient data
     *                                  according to the expression on the right hand side.
     * @param[in]                rhs    The right hand side expression of the assignment.
     *
     * @tparam Rhs The expression on the rhs of the statement.
     */
    template<typename Rhs>
    void store(Real& lhsValue, GradientData& lhsGradientData, const Rhs& rhs);

    /**
     * @brief Add a jacobi of 1.0 to the tape.
     *
     * The optimized version of push jacobi which signals the tape that the jacobi
     * corresponding to the tape is 1.0.
     *
     * @param[in,out]      data  A handle for data the tape can use for the evaluation.
     * @param[in]         value  The value of the active type which pushes the jacobi.
     * @param[in]  gradientData  The gradient data of the active type which pushes the jacobi.
     *
     * @tparam Data  The type of the data for the tape.
     */
    template<typename Data>
    void pushJacobi(Data& data, const Real& value, const GradientData& gradientData);

    /**
     * @brief Add a jacobi to the tape.
     *
     * The general version of push jacobi which signals the tape that the jacobi is used
     * in the evaluation and needs to be evaluated or stored.
     *
     * @param[in,out]      data  A handle for data the tape can use for the evaluation.
     * @param[in]        jacobi  The value of the jacobi.
     * @param[in]         value  The value of the active type which pushes the jacobi.
     * @param[in]  gradientData  The gradient data of the active type which pushes the jacobi.
     *
     * @tparam Data  The type of the data for the tape.
     */
    template<typename Data>
    void pushJacobi(Data& data, const Real& jacobi, const Real& value, const GradientData& gradientData);

    /**
     * @brief Called in the construction of a active type.
     *
     * The tape can initialize its gradient data for the active type.
     *
     * @param[in,out]        value  The value of the active type.
     * @param[in,out] gradientData  The gradient data which needs to be initialized.
     */
    virtual void initGradientData(Real& value, GradientData& gradientData) = 0;

    /**
     * @brief Destroy the gradient data of a active type.
     *
     * The tape can destroy the gradient data of the active type. This method is
     * called in the destructor of an active type.
     *
     * @param[in,out]        value  The value of the active type.
     * @param[in,out] gradientData  The gradient data which needs to be destroyed.
     */
    virtual void destroyGradientData(Real& value, GradientData& gradientData) = 0;

    /**
     * @brief Checks if all entries of the gradient are zero.
     *
     * @param[in] gradientData  The corresponding gradient data for the gradient.
     * @return true if all entries are zero.
     */
    virtual bool isGradientTotalZero(const GradientData& gradientData) = 0;

    /*
     * Access functions for the gradient information.
     */

    /**
     * @brief Set the gradient of the gradient data.
     *
     * The tape can set the gradient which corresponds to the gradient data.
     *
     * @param[in,out]  value  The gradient data of the active type.
     * @param[in]   gradient  The new gradient value for the active type.
     */
    virtual void setGradient(GradientData& value, const GradientValue& gradient) = 0;

    /**
     * @brief Get the gradient of the gradient data.
     *
     * The tape  returns the gradient which corresponds to the gradient data.
     *
     * @param[in] value  The gradient data of the active type.
     *
     * @return The gradient which belongs to the active type.
     */
    virtual GradientValue getGradient(const GradientData& value) const = 0;

    /**
     * @brief Get the gradient of the gradient data as a reference.
     *
     * The tape  returns the gradient which corresponds to the gradient data as a reference.
     *
     * @param[in] value  The gradient data of the active type.
     *
     * @return The gradient which belongs to the active type as a reference.
     */
    virtual GradientValue& gradient(GradientData& value) = 0;

    /**
     * @brief Get the gradient of the gradient data as a constant reference.
     *
     * The tape  returns the gradient which corresponds to the gradient data as a constant reference.
     *
     * @param[in] value  The gradient data of the active type.
     *
     * @return The gradient which belongs to the active type as a constant reference.
     */
    virtual const GradientValue& gradient(const GradientData& value) const = 0;

    /**
     * @brief Check if the gradient data has a nontrivial value.
     *
     * @param[in] value  The gradient data of the active type.
     *
     * @return Returns true if it has a nontrivial value. Otherwise returns false.
     */
    virtual bool isActive(const GradientData& value) const = 0;

  };
}
