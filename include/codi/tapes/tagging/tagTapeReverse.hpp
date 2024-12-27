/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../interfaces/fullTapeInterface.hpp"
#include "../misc/adjointVectorAccess.hpp"
#include "tagTapeBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Tape for tagging variables and find errors in the AD workflow.
   *
   * Mimics a CoDiPack reverse tape.
   *
   * See TagTapeBase for detailed information and functionality.
   *
   * @tparam T_Real  The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_tag   The type of the tag, usually int.
   */
  template<typename T_Real, typename T_Tag>
  struct TagTapeReverse : public FullTapeInterface<T_Real, T_Real, TagData<T_Tag>, EmptyPosition>,
                          public TagTapeBase<T_Real, T_Tag, T_Real, TagTapeReverse<T_Real, T_Tag>> {
      using Real = CODI_DD(T_Real, double);  ///< See TagTapeReverse.
      using Tag = CODI_DD(T_Tag, int);       ///< See TagTapeReverse.

      /// Required definition for event system.
      struct TapeTypes {
          /// Required definition for event system.
          struct IndexManager {
              /// Required definition for event system.
              using Index = int;
          };
      };

      using Gradient = Real;            ///< See TapeTypesInterface.
      using Identifier = TagData<Tag>;  ///< See TapeTypesInterface.
      using Position = EmptyPosition;   ///< See TapeTypesInterface.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      using Base = TagTapeBase<T_Real, T_Tag, T_Real, TagTapeReverse>;  ///< Base class abbreviation.

    private:
      bool active;  ///< Tape activity.

      Real tempPrimal;        ///< Temporary for primal values.
      Gradient tempGradient;  ///< Temporary for gradient values.

      std::set<TapeParameters> parameters;  ///< Temporary for tape parameters.

    public:

      /// Constructor.
      TagTapeReverse() : Base(), active(), tempPrimal(), tempGradient(), parameters() {}

      /*******************************************************************************/
      /// @name CustomAdjointVectorEvaluationTapeInterface interface implementation
      /// @{

      /// Do nothing.
      template<typename Adjoint>
      void evaluate(Position const& start, Position const& end, Adjoint* data) {
        CODI_UNUSED(start, end, data);
      }

      /// Do nothing.
      template<typename Adjoint>
      void evaluateForward(Position const& start, Position const& end, Adjoint* data) {
        CODI_UNUSED(start, end, data);
      }

      /// @}
      /*******************************************************************************/
      /// @name DataManagementTapeInterface interface implementation
      /// @{

      /// Do nothing.
      void writeToFile(std::string const& filename) const {
        CODI_UNUSED(filename);
      }

      /// Do nothing.
      void readFromFile(std::string const& filename) {
        CODI_UNUSED(filename);
      }

      /// Do nothing.
      void deleteData() {}

      /// Empty set.
      std::set<TapeParameters> const& getAvailableParameters() const {
        return parameters;
      }

      /// Do nothing.
      size_t getParameter(TapeParameters parameter) const {
        CODI_UNUSED(parameter);

        return 0;
      }

      /// Do nothing.
      bool hasParameter(TapeParameters parameter) const {
        CODI_UNUSED(parameter);

        return false;
      }

      /// Do nothing.
      void setParameter(TapeParameters parameter, size_t value) {
        CODI_UNUSED(parameter, value);
      }

      /// Do nothing.
      VectorAccessInterface<Real, Identifier>* createVectorAccess() {
        return nullptr;
      }

      /// Do nothing.
      template<typename Adjoint>
      VectorAccessInterface<Real, Identifier>* createVectorAccessCustomAdjoints(Adjoint* data) {
        CODI_UNUSED(data);
        return nullptr;
      }

      /// Do nothing.
      void deleteVectorAccess(VectorAccessInterface<Real, Identifier>* access) {
        delete access;
      }

      /// Swap members.
      void swap(TagTapeReverse& other) {
        std::swap(active, other.active);
        std::swap(parameters, other.parameters);
        Base::swap(other);
      }
      void resetHard() {}            ///< Do nothing.

      void deleteAdjointVector()   {} ///< Do nothing.
      void resizeAdjointVector()   {} ///< Do nothing.
      void beginUseAdjointVector() {} ///< Do nothing.
      void endUseAdjointVector()   {} ///< Do nothing.

      /// @}
      /*******************************************************************************/
      /// @name ExternalFunctionTapeInterface interface implementation
      /// @{

      /// Verifies tag properties.
      template<typename Lhs>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, TagTapeReverse, Lhs>& value) {
        registerInput(value);

        return Real();
      }

      /// Do nothing.
      void pushExternalFunction(ExternalFunction<TagTapeReverse> const& extFunc) {
        CODI_UNUSED(extFunc);
      }

      /// @}
      /*******************************************************************************/
      /// @name ForwardEvaluationTapeInterface interface implementation
      /// @{

      /// Do nothing.
      void evaluateForward(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);
      }

      /// Do nothing.
      void evaluateForward() {}

      /// @}
      /*******************************************************************************/
      /// @name GradientAccessTapeInterface interface implementation
      /// @{

      /// Verify tag.
      void setGradient(Identifier const& identifier, Gradient const& gradient,
                       AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(gradient, adjointsManagement);

        Base::verifyTagAndProperties(identifier.tag, 0.0, identifier.properties);
      }

      /// Verify tag.
      Gradient const& getGradient(Identifier const& identifier,
                                  AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const {
        CODI_UNUSED(adjointsManagement);

        Base::verifyTagAndProperties(identifier.tag, 0.0, identifier.properties);

        return tempGradient;
      }

      /// Verify tag.
      Gradient& gradient(Identifier const& identifier,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);

        Base::verifyTagAndProperties(identifier.tag, 0.0, identifier.properties);

        return tempGradient;
      }

      /// Verify tag.
      Gradient const& gradient(Identifier const& identifier,
                               AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const {
        CODI_UNUSED(adjointsManagement);

        Base::verifyTagAndProperties(identifier.tag, 0.0, identifier.properties);

        return tempGradient;
      }

      /// @}
      /*******************************************************************************/
      /// @name IdentifierInformationTapeInterface interface implementation
      /// @{

      /// Behave as linear index handler.
      static bool constexpr LinearIndexHandling = true;

      /// Zero tag.
      Identifier getPassiveIndex() const {
        return Identifier(Base::PassiveTag);
      }

      /// -1 tag.
      Identifier getInvalidIndex() const {
        return Identifier(Base::InvalidTag);
      }

      /// Verify tag.
      bool isIdentifierActive(Identifier const& index) const {
        return index.tag != Base::PassiveTag;
      }

      /// Set tag to passive.
      template<typename Lhs>
      void deactivateValue(LhsExpressionInterface<Real, Gradient, TagTapeReverse, Lhs>& value) {
        value.getIdentifier() = getPassiveIndex();
      }

      /// @}
      /*******************************************************************************/
      /// @name InternalStatementRecordingTapeInterface interface implementation
      /// @{

      /// Do not allow Jacobian optimization.
      static bool constexpr AllowJacobianOptimization = false;

      /// Do nothing.
      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);
        identifier = Identifier();
      }

      /// Do nothing.
      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value, identifier);
      }

      /// Verify all tags of the rhs and the lhs properties.
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, TagTapeReverse, Lhs>& lhs,
                             ExpressionInterface<Real, Rhs> const& rhs) {
        typename Base::ValidateTags validate;
        ValidationIndicator<Real, Tag> vi;

        validate.eval(rhs, vi, *this);

        Base::checkLhsError(lhs, rhs.cast().getValue());

        Base::handleError(vi);

        if (vi.isActive) {
          Base::setTag(lhs.cast().getIdentifier().tag);
        } else {
          Base::resetTag(lhs.cast().getIdentifier().tag);
        }
        lhs.cast().value() = rhs.cast().getValue();
      }

      /// Verify all tags of the rhs and the lhs properties.
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, TagTapeReverse, Lhs>& lhs,
                             LhsExpressionInterface<Real, Gradient, TagTapeReverse, Rhs> const& rhs) {
        store<Lhs, Rhs>(lhs, static_cast<ExpressionInterface<Real, Rhs> const&>(rhs));
      }

      /// Verify the lhs properties.
      template<typename Lhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, TagTapeReverse, Lhs>& lhs, PassiveReal const& rhs) {
        Base::checkLhsError(lhs, rhs);

        Base::resetTag(lhs.cast().getIdentifier().tag);

        lhs.cast().value() = rhs;
      }

      /// @}
      /*******************************************************************************/
      /// @name ManualStatementPushTapeInterface interface implementation
      /// @{

      /// Do nothing.
      void pushJacobiManual(Real const& jacobian, Real const& value, Identifier const& index) {
        CODI_UNUSED(jacobian, value, index);
      }

      /// Set tag on lhs.
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue, size);

        Base::checkLhsError(lhsValue, lhsIndex, lhsValue);
        setTag(lhsIndex.tag);
      }

      /// @}
      /*******************************************************************************/
      /// @name PositionalEvaluationTapeInterface interface implementation
      /// @{

      /// Do nothing.
      void evaluate(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);
      }

      /// Do nothing.
      void clearAdjoints(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);
      }

      /// Do nothing.
      Position getPosition() const {
        return Position();
      }

      /// Do nothing.
      Position getZeroPosition() const {
        return Position();
      }

      /// Do nothing.
      void resetTo(Position const& pos, bool resetAdjoints = true) {
        CODI_UNUSED(pos, resetAdjoints);
      }

      /// @}
      /*******************************************************************************/
      /// @name PreaccumulationEvaluationTapeInterface interface implementation
      /// @{

      /// Do nothing.
      void evaluateKeepState(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);
      }
      /// Do nothing.
      void evaluateForwardKeepState(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);
      }

      /// @}
      /*******************************************************************************/
      /// @name PrimalEvaluationTapeInterface interface implementation
      /// @{

      static bool constexpr HasPrimalValues = false;        ///< No primal values.
      static bool constexpr RequiresPrimalRestore = false;  ///< No primal values.

      /// Do nothing.
      void evaluatePrimal(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);
      }

      /// Do nothing.
      void evaluatePrimal() {}

      /// Do nothing.
      void setPrimal(Identifier const& identifier, Real const& gradient) {
        CODI_UNUSED(identifier, gradient);
      }

      /// Do nothing.
      Real const& getPrimal(Identifier const& identifier) const {
        CODI_UNUSED(identifier);
        return tempPrimal;
      }

      /// Do nothing.
      Real& primal(Identifier const& identifier) {
        CODI_UNUSED(identifier);
        return tempPrimal;
      }

      /// Do nothing.
      Real const& primal(Identifier const& identifier) const {
        CODI_UNUSED(identifier);
        return tempPrimal;
      }

      /// Do nothing.
      void revertPrimals(Position const& pos) {
        CODI_UNUSED(pos);
      }

      /// @}
      /*******************************************************************************/
      /// @name ReverseTapeInterface interface implementation
      /// @{

      /// Verify value properties.
      template<typename Lhs>
      void registerInput(LhsExpressionInterface<Real, Gradient, TagTapeReverse, Lhs>& value) {
        Base::setTag(value.cast().getIdentifier().tag);
        Base::verifyRegisterValue(value, value.cast().getIdentifier());  // verification is mainly for the properties
      }

      /// Verify tag.
      template<typename Lhs>
      void registerOutput(LhsExpressionInterface<Real, Gradient, TagTapeReverse, Lhs>& value) {
        Base::verifyRegisterValue(value, value.cast().getIdentifier());
      }

      /// Set tape to active.
      void setActive() {
        active = true;
      }

      /// Set tape to passive.
      void setPassive() {
        active = false;
      }

      /// Check if tape is active.
      bool isActive() const {
        return active;
      }

      /// Default check.
      bool isActive(Identifier const& identifier) const {
        return identifier.tag != Base::PassiveTag;
      }

      void evaluate() {}  ///< Do nothing.

      /// Do nothing.
      void clearAdjoints() {}

      /// Do nothing.
      void reset(bool resetAdjoints = true) {
        CODI_UNUSED(resetAdjoints);
      }

      /// Do nothing.
      template<typename Stream = std::ostream>
      void printStatistics(Stream& out = std::cout) const {
        CODI_UNUSED(out);
      }

      /// Do nothing.
      template<typename Stream = std::ostream>
      void printTableHeader(Stream& out = std::cout) const {
        CODI_UNUSED(out);
      }

      /// Do nothing.
      template<typename Stream = std::ostream>
      void printTableRow(Stream& out = std::cout) const {
        CODI_UNUSED(out);
      }

      /// Do nothing.
      TapeValues getTapeValues() const {
        return TapeValues("TagTapeReverse");
      }

      /// @}
  };
}
