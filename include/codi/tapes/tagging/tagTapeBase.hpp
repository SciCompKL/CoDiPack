/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../expressions/logic/helpers/forEachLeafLogic.hpp"
#include "../../expressions/aggregate/aggregatedActiveType.hpp"
#include "../../misc/enumBitset.hpp"
#include "../indices/indexManagerInterface.hpp"
#include "../interfaces/fullTapeInterface.hpp"
#include "../misc/adjointVectorAccess.hpp"
#include "tagData.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper class for statement validation.
  template<typename Real, typename Tag>
  struct ValidationEntry {
      bool isTagError;        ///< true if the tag on the value is wrong.
      bool isPropertyError;   ///< true if a property on the value is wrong.
      Tag tagError;           ///< Wrong tag.
      TagFlags propertyError; ///< Wrong property.
      Real const* value;      ///< Address of the offending value.
      Real newValue;          ///< New value for property errors.
  };

  /// Helper class for statement validation.
  template<typename Real, typename Tag>
  struct ValidationIndicator {
      /// Stores all problematic entries
      std::vector<ValidationEntry<Real, Tag>> errorEntries = {};

      /// Add a new entry.
      void addEntry(ValidationEntry<Real, Tag> const& entry) {
        errorEntries.push_back(entry);
      }

      /// Clear all errors.
      void reset() {
        errorEntries.resize(0);
      }

      /// Test if there is an error.
      bool hasError() {
        return !errorEntries.empty();
      }
  };

  /**
   * @brief Base implementation for tagging tapes.
   *
   * Provides all basic management routines for the tag.
   *
   * See tests/functional/src/testTagging.cpp for an example.
   *
   * @tparam T_Real  The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_tag   The type of the tag, usually int.
   */
  template<typename T_Real, typename T_Tag, typename T_Gradient, typename T_Impl>
  struct TagTapeBase {
      using Real = CODI_DD(T_Real, double);       ///< See TagTapeBase.
      using Tag = CODI_DD(T_Tag, int);            ///< See TagTapeBase.
      using Gradient = CODI_DD(T_Gradient, int);  ///< See TagTapeBase.
      using Impl = CODI_DD(T_Impl, TagTapeBase);  ///< See TagTapeBase.

      using Identifier = int;                   ///< See TapeTypesInterface.
      using ActiveTypeTapeData = TagData<Tag>;  ///< See TapeTypesInterface.

      /// Callback for a change in a lhs value.
      using TagPropertyErrorCallback = void (*)(Real const& currentValue, Real const& newValue, TagFlags flag,
                                                void* userData);

      /// Callback for a tag error.
      using TagErrorCallback = void (*)(Tag const& correctTag, Tag const& wrongTag, void* userData);

      /// Callback for all errors.
      using ErrorCallback = void (*)(ValidationEntry<Real, Tag> const& data, Tag const& curTag, void* userData);

      static Tag constexpr PassiveTag = Tag(0);   ///< Tag indicating an inactive value.
      static Tag constexpr InvalidTag = Tag(-1);  ///< Tag indicating an invalid value.

    protected:

      Tag curTag;  ///< Current tag for new values.

      TagPropertyErrorCallback tagPropertyErrorCallback;  ///< User defined callback for lhs value errors.
      void* tagPropertyErrorUserData;                     ///< User data in call to callback for lhs value errors.

      TagErrorCallback tagErrorCallback;  ///< User defined callback for tag errors.
      void* tagErrorUserData;             ///< User data in call to callback for tag errors.

      ErrorCallback errorCallback;  ///< User defined callback for tag and property errors. Has priority over the other callbacks.
      void* errorUserData;          ///< User data in call to callback for errors.

      bool preaccumulationHandling;  ///< Parameter to enable/disable preaccumulation handling.
      Tag preaccumulationTag;        ///< Tag used for preaccumulation specialized handling.

      Identifier tempIdentifier = 0;           ///< Temporary for identifier values.
      Identifier const activeIdentifier = 1;   ///< Temporary for active identifier.
      Identifier const passiveIdentifier = 0;  ///< Temporary for passive identifier.

      ValidationIndicator<Real, Tag> vi = {}; ///< Helper for error detection.

    public:

      /// Constructor.
      TagTapeBase()
          : curTag(),
            tagPropertyErrorCallback(defaultPropertyErrorCallback),
            tagPropertyErrorUserData(nullptr),
            tagErrorCallback(defaultTagErrorCallback),
            tagErrorUserData(this),
            errorCallback(nullptr),
            errorUserData(nullptr),
            preaccumulationHandling(true),
            preaccumulationTag(1337) {}

      /// Looks at the tags for the expression.
      struct ValidateTags : public ForEachLeafLogic<ValidateTags> {
        public:

          /// \copydoc codi::ForEachLeafLogic::handleActive
          template<typename Node>
          CODI_INLINE void handleActive(Node const& node, ValidationIndicator<Real, Tag>& vi, Impl& tape, bool& isActive) {
            ActiveTypeTapeData tagData = node.getTapeData();

            tape.verifyRhsAccess(node.getValue(), tagData, vi);

            isActive |= tagData.tag != 0;
          }
      };

      /// Swap members.
      void swap(Impl& other) {
        std::swap(curTag, other.curTag);
        std::swap(tagPropertyErrorCallback, other.tagPropertyErrorCallback);
        std::swap(tagPropertyErrorUserData, other.tagPropertyErrorUserData);
        std::swap(tagErrorCallback, other.tagErrorCallback);
        std::swap(tagErrorUserData, other.tagErrorUserData);
        std::swap(preaccumulationHandling, other.preaccumulationHandling);
        std::swap(preaccumulationTag, other.preaccumulationTag);
        std::swap(vi, other.vi);
      }

      /*******************************************************************************/
      /// @name IdentifierInformationTapeInterface interface partial implementation
      /// @{

      /// \copydoc codi::IdentifierInformationTapeInterface::getIdentifier()
      CODI_INLINE Identifier const& getIdentifier(ActiveTypeTapeData const& data) {
        if (data.tag != 0) {
          return activeIdentifier;
        } else {
          return passiveIdentifier;
        }
      }

      /// \copydoc codi::IdentifierInformationTapeInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier(ActiveTypeTapeData& data) {
        tempIdentifier = getIdentifier(std::as_const(data));
        return tempIdentifier;
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of InternalStatementRecordingTapeInterface
      /// @{

      static bool constexpr AllowJacobianOptimization = false;  ///< Do not allow Jacobian optimization.

      /// Do nothing.
      template<typename Real>
      void initTapeData(Real& value, ActiveTypeTapeData& data) {
        CODI_UNUSED(value);
        data = ActiveTypeTapeData();
      }

      /// Do nothing.
      template<typename Real>
      void destroyTapeData(Real& value, ActiveTypeTapeData& data) {
        CODI_UNUSED(value, data);
      }

      /// Verify all tags of the rhs and properties from the lhs and rhs.
      template<typename Aggregated, typename Type, typename Lhs, typename Rhs>
      CODI_INLINE void store(AggregatedActiveType<Aggregated, Type, Lhs>& lhs,
                             ExpressionInterface<Aggregated, Rhs> const& rhs) {
        using AggregatedTraits = RealTraits::AggregatedTypeTraits<Aggregated>;

        int constexpr Elements = AggregatedTraits::Elements;

        std::array<bool, Elements> isActive = {};

        Aggregated real = rhs.cast().getValue();

        static_for<Elements>([this, &isActive, &lhs, &rhs, &real](auto i) CODI_LAMBDA_INLINE {
          ValidateTags validate;

          validate.eval(ArrayAccessExpression<Aggregated, i.value, Rhs>(rhs), vi, cast(), isActive[i.value]);
          verifyLhsWrite(lhs.values[i.value], AggregatedTraits::template arrayAccess<i.value>(real));

          handleError(vi);
        });

        static_for<Elements>([this, &isActive, &lhs, &real](auto i) CODI_LAMBDA_INLINE {
          if (isActive[i.value]) {
            setTag(lhs.values[i.value].getTapeData().tag);
          } else {
            resetTag(lhs.values[i.value].getTapeData().tag);
          }
          lhs.values[i.value].getTapeData().properties.reset(TagFlags::DoNotUse);
          lhs.values[i.value].value() = AggregatedTraits::template arrayAccess<i.value>(real);
        });
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Optimization for copy statements of aggregated types.
      template<typename Aggregated, typename Type, typename Lhs, typename Rhs>
      CODI_INLINE void store(AggregatedActiveType<Aggregated, Type, Lhs>& lhs,
                             AggregatedActiveType<Aggregated, Type, Rhs> const& rhs) {
        store<Aggregated, Type, Lhs, Rhs>(lhs, static_cast<ExpressionInterface<Aggregated, Rhs> const&>(rhs));
      }

      /// Verify all tags of the rhs and properties from the lhs and rhs.
      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, ExpressionInterface<Real, Rhs> const& rhs) {
        ValidateTags validate;
        bool isActive = false;
        validate.eval(rhs, vi, cast(), isActive);

        verifyLhsWrite(lhs.cast().getValue(), lhs.cast().getTapeData(), rhs.cast().getValue(), vi);

        handleError(vi);

        if (isActive) {
          setTag(lhs.cast().getTapeData().tag);
        } else {
          resetTag(lhs.cast().getTapeData().tag);
        }
        lhs.cast().getTapeData().properties.reset(TagFlags::DoNotUse);

        lhs.cast().value() = rhs.cast().getValue();
      }

      /// Verify all tags of the rhs and the lhs properties.
      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, Impl, Rhs> const& rhs) {
        store<Lhs, Rhs>(lhs, static_cast<ExpressionInterface<Real, Rhs> const&>(rhs));
      }

      /// Verify the lhs properties.
      template<typename Lhs>
      void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, Real const& rhs) {
        verifyLhsWrite(lhs.cast().getValue(), lhs.cast().getTapeData(), rhs, vi);
        handleError(vi);

        resetTag(lhs.cast().getTapeData().tag);

        lhs.cast().value() = rhs;
      }

      /// @}
      /*******************************************************************************/
      /// @name Tagging specific functions.
      /// @{

      /// Set the current tag of the tape.
      void setCurTag(const Tag& tag) {
        this->curTag = tag;
      }

      /// Get the current tag of the tape.
      Tag getCurTag() {
        return this->curTag;
      }

      /// Get tag of a CoDiPack active type.
      template<typename Lhs>
      Tag getTagFromVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        return value.cast().getTapeData().tag;
      }

      /// Set tag on a CoDiPack active type.
      template<typename Lhs>
      void setTagOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs> const& value) {
        value.cast().getTapeData().tag = this->curTag;
      }

      /// Clear tag on a CoDiPack active type.
      template<typename Lhs>
      void clearTagOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        value.cast().getTapeData().tag = Tag();
      }

      /// Clear properties on a CoDiPack active type.
      template<typename Lhs>
      void clearTagPropertiesOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        value.cast().getTapeData().properties.reset();
      }

      /// Set properties on a CoDiPack active type.
      template<typename Lhs>
      void setTagPropertyOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value, TagFlags flag) {
        value.cast().getTapeData().properties.set(flag);
      }

      /// Check properties on a CoDiPack active type.
      template<typename Lhs>
      bool hasTagPropertyOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value, TagFlags flag) {
        return value.cast().getTapeData().properties.test(flag);
      }

      /// Set the callback and user data for a property error error.
      void setTagPropertyErrorCallback(TagPropertyErrorCallback const& callback, void* userData) {
        tagPropertyErrorCallback = callback;
        tagPropertyErrorUserData = userData;
      }

      /// Set the callback and user data for a tag error.
      void setTagErrorCallback(TagErrorCallback const& callback, void* userData) {
        tagErrorCallback = callback;
        tagErrorUserData = userData;
      }

      /// Set a general error handle for both, tag and property, errors. If this error handler is set,
      /// the others are no longer called.
      void setErrorCallback(ErrorCallback const& callback, void* userData) {
        errorCallback = callback;
        errorUserData = userData;
      }

      /// @brief Enable or disable specialized handling for preaccumulation. Default: true
      /// Uses a special tag to sanitize preaccumulation regions.
      CODI_INLINE void setPreaccumulationHandlingEnabled(bool enabled) {
        preaccumulationHandling = enabled;
      }

      /// Set the special tag for preaccumulation regions. See setPreaccumulationHandlingEnabled().
      CODI_INLINE void setPreaccumulationHandlingTag(Tag tag) {
        preaccumulationTag = tag;
      }

      /// If handling for preaccumulation is enabled.
      CODI_INLINE bool isPreaccumulationHandlingEnabled() {
        return preaccumulationHandling;
      }

      /// The special tag for preaccumulation.
      CODI_INLINE Tag getPreaccumulationHandlingTag() {
        return preaccumulationTag;
      }

    protected:

      /// Create a tag error entry.
      ValidationEntry<Real, Tag> makeTagError(Tag tag, Real const* value) {
        return {true, false, tag, TagFlags::MaxElement, value, Real()};
      }

      /// Create a property error entry.
      ValidationEntry<Real, Tag> makeProeprtyError(TagFlags flag, Real const* value, Real newValue = {}) {
        return {false, true, 0, flag, value, newValue};
      }

      /// Check tag and property errors.
      void verifyRhsAccess(Real const& value, ActiveTypeTapeData const& data, ValidationIndicator<Real, Tag>& vi) {
        verifyRhsAccessForTag(value, data, vi);
        verifyRhsAccessForProperty(value, data, vi);
      }

      /// Check tag errors. Raises an error for tags that are not equal to the current tag. Ignores passive tags.
      void verifyRhsAccessForTag(Real const& value, ActiveTypeTapeData const& data, ValidationIndicator<Real, Tag>& vi) {
        if (PassiveTag != data.tag && InvalidTag != data.tag) {
          if (data.tag != curTag) {
            vi.addEntry(makeTagError(data.tag, &value));
          }
        }
      }

      /// Check properties for rhs values. Validates the 'DoNotUse' property.
      void verifyRhsAccessForProperty(Real const& value, ActiveTypeTapeData const& data, ValidationIndicator<Real, Tag>& vi) {
        if (data.properties.test(TagFlags::DoNotUse)) {
          vi.addEntry(makeProeprtyError(TagFlags::DoNotUse, &value));
        }
      }

      /// Validates the lhs properties.
      CODI_INLINE void verifyLhsWrite(Real const& lhsValue, ActiveTypeTapeData const& lhsData, Real const& rhsValue, ValidationIndicator<Real, Tag>& vi) {
        if (lhsData.properties.test(TagFlags::DoNotChange)) {
          if (lhsValue != rhsValue) {
            vi.addEntry(makeProeprtyError(TagFlags::DoNotChange, &lhsValue, rhsValue));
          }
        } else if (lhsData.properties.test(TagFlags::DoNotWrite)) {
          vi.addEntry(makeProeprtyError(TagFlags::DoNotWrite, &lhsValue));
        }
      }

      /// Verify tag, properties and lhs error.
      template<typename Lhs>
      CODI_INLINE void verifyRegisterValue(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value,
                                           ActiveTypeTapeData const& tag) {
        verifyRhsAccess(value.cast().getValue(), value.cast().getTapeData(), vi);
        verifyLhsWrite(value.cast().getValue(), value.cast().getTapeData(), value.cast().getValue(), vi);

        handleError(vi);
      }

      /// Default callback for TagPropertyErrorCallback.
      static void defaultPropertyErrorCallback(Real const& currentValue, Real const& newValue, TagFlags flag,
                                               void* userData) {
        CODI_UNUSED(userData);

        std::cerr << "Property error '" << std::to_string(flag) << "' on value. current value: " << currentValue
                  << " new value: " << newValue << "" << std::endl;
      }

      /// Default callback for TagErrorCallback.
      static void defaultTagErrorCallback(Tag const& correctTag, Tag const& wrongTag, void* userData) {
        TagTapeBase& impl = *static_cast<TagTapeBase*>(userData);

        // output default warning if no handle is defined.
        std::cerr << "Use of variable with bad tag '" << wrongTag << "', should be '" << correctTag << "'.";
        if (wrongTag == impl.preaccumulationTag) {
          std::cerr << " The value seems to be a preaccumulation output.";
        } else if (correctTag == impl.preaccumulationTag) {
          std::cerr << " The value seems to be used during a preaccumulation but is not declared as an input.";
        }
        std::cerr << std::endl;
      }

      /// Default callback for ErrorCallback. (Currently not used)
      static void defaultErrorCallback(ValidationEntry<Real, Tag> const& data, Tag const& curTag, void* userData) {
        if (data.isTagError) {
          defaultTagErrorCallback(curTag, data.tagError, userData);
        }
        if (data.isPropertyError) {
          defaultPropertyErrorCallback(*data.value, data.newValue, data.propertyError, userData);
        }
      }

      /// Call tag error callback.
      CODI_INLINE void handleError(ValidationIndicator<Real, Tag>& vi) const {
        if (vi.hasError()) {
          for(ValidationEntry<Real, Tag> const& entry : vi.errorEntries) {
            if(nullptr != errorCallback) {
              errorCallback(entry, curTag, errorUserData);
            } else {
              if (entry.isTagError) {
                tagErrorCallback(curTag, entry.tagError, tagErrorUserData);
              }
              if (entry.isPropertyError) {
                tagPropertyErrorCallback(*entry.value, entry.newValue, entry.propertyError, tagPropertyErrorUserData);
              }
            }
          }
        }

        vi.reset();
      }

      /// Set tag on value.
      CODI_INLINE void setTag(Tag& tag) const {
        tag = curTag;
      }

      /// Reset tag on value.
      CODI_INLINE void resetTag(Tag& tag) const {
        tag = Tag();
      }
      /// @}

    private:
      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
