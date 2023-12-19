/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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

#include "../../expressions/logic/helpers/forEachLeafLogic.hpp"
#include "../../misc/enumBitset.hpp"
#include "../indices/indexManagerInterface.hpp"
#include "../interfaces/fullTapeInterface.hpp"
#include "../misc/adjointVectorAccess.hpp"
#include "tagData.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper class for statement validation.
  template<typename Tag>
  struct ValidationIndicator {
      bool isActive;     ///< true if an active rhs is detected. tag != 0
      bool hasError;     ///< true if an error is detected.
      bool hasTagError;  ///< true if a tag not the current required tag.
      bool hasUseError;  ///< true if a value is used in the wrong way.
      Tag errorTag;      ///< Value of the wrong tag.

      /// Constructor.
      ValidationIndicator() : isActive(false), hasError(false), hasTagError(false), hasUseError(false), errorTag() {}
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
      using Impl = CODI_DD(T_Impl, int);          ///< See TagTapeBase.

      using Identifier = TagData<Tag>;  ///< See TapeTypesInterface.

      /// Callback for a change in a lhs value.
      using TagLhsChangeErrorCallback = void (*)(Real const& currentValue, Real const& newValue, void* userData);

      /// Callback for a tag error.
      using TagErrorCallback = void (*)(Tag const& correctTag, Tag const& wrongTag, bool tagError, bool useError,
                                        void* userData);

      static Tag constexpr PassiveTag = Tag(0);  ///< Tag indicating an inactive value.

    protected:

      Tag curTag;  ///< Current tag for new values.

      TagLhsChangeErrorCallback tagLhsChangeErrorCallback;  ///< User defined callback for lhs value errors.
      void* tagChangeErrorUserData;                         ///< User data in call to callback for lhs value errors.

      TagErrorCallback tagErrorCallback;  ///< User defined callback for tag errors.
      void* tagErrorUserData;             ///< User data in call to callback for tag errors.

      bool preaccumulationHandling;  ///< Parameter to enable/disable preaccumulation handling.
      Tag preaccumulationTag;        ///< Tag used for preaccumulation specialized handling.

    public:

      /// Constructor.
      TagTapeBase()
          : curTag(),
            tagLhsChangeErrorCallback(defaultTagLhsChangeErrorCallback),
            tagChangeErrorUserData(nullptr),
            tagErrorCallback(defaultTagErrorCallback),
            tagErrorUserData(this),
            preaccumulationHandling(true),
            preaccumulationTag(1337) {}

      /// Looks at the tags for the expression.
      struct ValidateTags : public ForEachLeafLogic<ValidateTags> {
        public:

          /// \copydoc codi::ForEachLeafLogic::handleActive
          template<typename Node>
          CODI_INLINE void handleActive(Node const& node, ValidationIndicator<Tag>& vi, Impl& tape) {
            Identifier tagData = node.getIdentifier();
            tape.verifyTag(vi, tagData.tag);
            tape.verifyProperties(vi, tagData.properties);
          }
      };

      /// Swap members.
      void swap(Impl& other) {
        std::swap(curTag, other.curTag);
        std::swap(tagLhsChangeErrorCallback, other.tagLhsChangeErrorCallback);
        std::swap(tagChangeErrorUserData, other.tagChangeErrorUserData);
        std::swap(tagErrorCallback, other.tagErrorCallback);
        std::swap(tagErrorUserData, other.tagErrorUserData);
        std::swap(preaccumulationHandling, other.preaccumulationHandling);
        std::swap(preaccumulationTag, other.preaccumulationTag);
      }

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
        return value.cast().getIdentifier().tag;
      }

      /// Set tag on a CoDiPack active type.
      template<typename Lhs>
      void setTagOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        value.cast().getIdentifier().tag = this->curTag;
      }

      /// Clear tag on a CoDiPack active type.
      template<typename Lhs>
      void clearTagOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        value.cast().getIdentifier().tag = Tag();
      }

      /// Clear properties on a CoDiPack active type.
      template<typename Lhs>
      void clearTagPropertiesOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        value.cast().getIdentifier().properties.reset();
      }

      /// Set properties on a CoDiPack active type.
      template<typename Lhs>
      void setTagPropertyOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value, TagFlags flag) {
        value.cast().getIdentifier().properties.set(flag);
      }

      /// Check properties on a CoDiPack active type.
      template<typename Lhs>
      bool hasTagPropertyOnVariable(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value, TagFlags flag) {
        return value.cast().getIdentifier().properties.test(flag);
      }

      /// Set the callback and user data for a lhs error.
      void setTagLhsChangeErrorCallback(TagLhsChangeErrorCallback const& callback, void* userData) {
        tagLhsChangeErrorCallback = callback;
        tagChangeErrorUserData = userData;
      }

      /// Set the callback and user data for a tag error.
      void setTagErrorCallback(TagErrorCallback const& callback, void* userData) {
        tagErrorCallback = callback;
        tagErrorUserData = userData;
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

      /// Checks if the tag is correct. Errors are set on the ValidationIndicator object.
      CODI_INLINE void verifyTag(ValidationIndicator<Tag>& vi, Tag const& tag) const {
        if (PassiveTag != tag) {
          vi.isActive = true;
          if (tag != curTag) {
            vi.hasError = true;
            vi.hasTagError = true;
            vi.errorTag = tag;
          }
        }
      }

      /// Checks if the tag is correct and creates an error.
      CODI_INLINE void verifyTag(Tag const& tag) const {
        ValidationIndicator<Tag> vi;

        verifyTag(vi, tag);
        handleError(vi);
      }

      /// Checks if the tag properties are correct.
      CODI_INLINE void verifyProperties(ValidationIndicator<Tag>& vi, const EnumBitset<TagFlags>& properties) const {
        if (properties.test(TagFlags::DoNotUse)) {
          vi.hasError = true;
          vi.hasUseError = true;
        }
      }

      /// Checks if the tag and the properties are correct.
      CODI_INLINE void verifyTagAndProperties(Tag const& tag, const EnumBitset<TagFlags>& properties) const {
        ValidationIndicator<Tag> vi;

        verifyTag(vi, tag);
        verifyProperties(vi, properties);
        handleError(vi);
      }

      /// Default callback for TagLhsChangeErrorCallback.
      static void defaultTagLhsChangeErrorCallback(Real const& currentValue, Real const& newValue, void* userData) {
        CODI_UNUSED(userData);

        std::cerr << "DoNotChange variable changes value from '" << currentValue << "' to '" << newValue << "'." << std::endl;
      }

      /// Default callback for TagErrorCallback.
      static void defaultTagErrorCallback(Tag const& correctTag, Tag const& wrongTag, bool tagError, bool useError,
                                          void* userData) {

        TagTapeBase& impl = *static_cast<TagTapeBase*>(userData);

        // output default warning if no handle is defined.
        if (useError) {
          std::cerr << "DoNotUse variable is used." << std::endl;
        }
        if (tagError) {
          std::cerr << "Use of variable with bad tag '" << wrongTag << "', should be '" << correctTag << "'.";
          if(wrongTag == impl.preaccumulationTag) {
            std::cerr << " The value seems to be a preaccumulation output.";
          } else if(correctTag == impl.preaccumulationTag) {
            std::cerr << " The value seems to be used during a preaccumulation but is not declared as an input.";
          }
          std::cerr << std::endl;
        }
      }

      /// Check if the lhs value is changed.
      CODI_INLINE void checkLhsError(Real& lhsValue, Identifier& lhsIdentifier, const Real& rhs) const {
        if (lhsIdentifier.properties.test(TagFlags::DoNotChange)) {
          if (lhsValue != rhs) {
            tagLhsChangeErrorCallback(lhsValue, rhs, tagChangeErrorUserData);
          }
        }
      }

      /// Check if the lhs value is changed.
      template<typename Lhs>
      CODI_INLINE void checkLhsError(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, const Real& rhs) const {
        checkLhsError(lhs.cast().value(), lhs.cast().getIdentifier(), rhs);
      }

      /// Call tag error callback.
      CODI_INLINE void handleError(ValidationIndicator<Tag>& vi) const {
        if (vi.hasError) {
          tagErrorCallback(curTag, vi.errorTag, vi.hasTagError, vi.hasUseError, tagErrorUserData);
        }
      }

      /// Verify tag, properties and lhs error.
      template<typename Lhs>
      CODI_INLINE void verifyRegisterValue(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value,
                                           const Identifier& tag) {
        ValidationIndicator<Tag> vi;

        verifyTag(vi, tag.tag);
        verifyProperties(vi, tag.properties);
        handleError(vi);

        checkLhsError(value, value.cast().getValue());
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
  };
}
