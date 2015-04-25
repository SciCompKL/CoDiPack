#pragma once

#include "../configure.h"

namespace codi {

  /**
   * @brief Interface common to all tapes.
   *
   * The basic interface each tape has to implement. It defines functions which are used
   * by the active type to signal the tape when an active type is created or destroyed. For
   * each operation the active type calls also the store operation to inform the tape that
   * an expression is assigned to an active type.
   *
   * @template GradientDataType  The data the tape uses to identify each active variable
   *                             and where the tape can store information about the
   *                             gradient.
   */
  template <typename GradientDataType>
  class TapeInterface {
  public:

    /**
     * @brief The data for the gradient information of the tape.
     *
     * Each tape can define a data for the gradient which each active type will
     * define in its own class. The tape can use this data to identify each active
     * type and compute the gradient information.
     */
    typedef GradientDataType GradientData;

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
     */
    template<typename Rhs>
    void store(Real& lhsValue, GradientData& lhsGradientData, const Rhs& rhs);

    /**
     * @brief Add a jacobi of 1.0 to the tape.
     *
     * The optimized version of push jacobi which signals the tape that the jacobi
     * corresponding to the tape is 1.0.
     *
     * @param[in/out]   gradient  A handle to the gradient of the operation. Mostly used for the forward mode.
     * @param[in]          value  The value of the active type which pushes the jacobi.
     * @param[in]   gradientData  The gradient data of the active type which pushes the jacobi.
     */
    virtual void pushJacobi(Real& gradient, const Real& value, const GradientData& gradientData) = 0;

    /**
     * @brief Add a jacobi to the tape.
     *
     * The general version of push jacobi which signals the tape that the jacobi is used
     * in the evaluation and needs to be evaluated or stored.
     *
     * @param[in/out]   gradient  A handle to the gradient of the operation. Mostly used for the forward mode.
     * @param[in]          value  The value of the active type which pushes the jacobi.
     * @param[in]   gradientData  The gradient data of the active type which pushes the jacobi.
     */
    virtual void pushJacobi(Real& gradient, const Real& jacobi, const Real& value, const GradientData& gradientData) = 0;

    /**
     * @brief Called in the construction of a active type.
     *
     * The tape can initialize its gradient data for the active type.
     *
     * @param[in/out]        value  The value of the active type.
     * @param[in/out] gradientData  The gradient data which needs to be initialized.
     */
    virtual void initGradientData(Real& value, GradientData& gradientData) = 0;

    /**
     * @brief Destroy the gradient data of a active type.
     *
     * The tape can destory the gradient data of the active type. This method is
     * called in the destructor of an active type.
     *
     * @param[in/out]        value  The value of the active type.
     * @param[in/out] gradientData  The gradient data which needs to be destroyed.
     */
    virtual void destroyGradientData(Real& value, GradientData& gradientData) = 0;


    /**
     * Access functions for the gradient information.
     */

    /**
     * @brief Set the gradient of the gradient data.
     *
     * The tape can set the gradient which corresponds to the gradient data.
     *
     * @param[in/out]  value  The gradient data of the active type.
     * @param[in]   gradient  The new gradient value for the active type.
     */
    virtual void setGradient(GradientData& value, const Real& gradient) = 0;

    /**
     * @brief Get the gradient of the gradient data.
     *
     * The tape  returns the gradient which corresponds to the gradient data.
     *
     * @param[in] value  The gradient data of the active type.
     *
     * @return The gradient which belongs to the active type.
     */
    virtual Real getGradient(const GradientData& value) const = 0;

  };
}