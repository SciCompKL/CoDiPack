
#include <codi.hpp>
#include <iostream>
#include <Eigen/Eigen>

//! [Tutorial 7 - Aggregated type implementation]

//! [Base declaration]
// Declare all types without the CoDiPack active type.
using Number = double;
using Vector4 = Eigen::Matrix<Number, 4, 1>;
using Vector4Transpose = Eigen::Matrix<Number, 1, 4>;

// Declare all types with the CoDiPack active type.
using ActiveNumber = codi::RealReverse;
//! [Base declaration]

// Start with the implementation of aggregated types.

// 1. Specialize the necessary traits.
//! [Vector4 specializations]
template<>
struct codi::RealTraits::AggregatedTypeTraits<Vector4>
    : public codi::RealTraits::ArrayAggregatedTypeTraitsBase<Vector4, Number, Vector4, 4> {};

namespace codi {
  CODI_CREATE_TRANSPOSE(Vector4, Vector4Transpose, jacobian.transpose()); ///< Transpose operation for the vector type.
}
//! [Vector4 specializations]

// 2. Implement CoDiPack expression operations for the vector type.

//! [Operation implementation]
/// Operation for scalar * vector.
template<typename T> // Template argument not used, we use the direct type specification of Vector4.
struct ScalarVectorMultiplicationOperation : public codi::BinaryJacobianOperation<Vector4, ScalarVectorMultiplicationOperation<T>> {

    /// Primal evaluation.
    static CODI_INLINE Vector4 primal(Number const& s, Vector4 const& v) {
      return s * v;
    }

    /// Gradient with respect to s.
    static CODI_INLINE Vector4 gradientA(Number const& s, Vector4 const& v, Vector4 const& result) {
      return v;
    }

    /// Gradient with respect to v.
    static CODI_INLINE Number gradientB(Number const& s, Vector4 const& v, Vector4 const& result) {
      return s;
    }
};

// Define scalar * vector overload.
template<typename ArgS, typename ArgV>
auto operator*(codi::ExpressionInterface<Number, ArgS> const& s, codi::ExpressionInterface<Vector4, ArgV> const& v) {
  return codi::ComputeExpression<Vector4, ScalarVectorMultiplicationOperation, ArgS, ArgV>(s.cast(), v.cast());
}

/// Operation for vector + vector.
template<typename T> // Template argument not used, we use the direct type specification of Vector4.
struct VectorAdditionOperation : public codi::BinaryJacobianOperation<Vector4, VectorAdditionOperation<T>> {

    /// Primal evaluation.
    static CODI_INLINE Vector4 primal(Vector4 const& v1, Vector4 const& v2) {
      return v1 + v2;
    }

    /// Gradient with respect to v1.
    static CODI_INLINE Number gradientA(Vector4 const& v1, Vector4 const& v2, Vector4 const& result) {
      return 1.0;
    }

    /// Gradient with respect to v2.
    static CODI_INLINE Number gradientB(Vector4 const& v1, Vector4 const& v2, Vector4 const& result) {
      return 1.0;
    }
};

// Define vector + vector overload.
template<typename ArgV1, typename ArgV2>
auto operator+(codi::ExpressionInterface<Vector4, ArgV1> const& v1, codi::ExpressionInterface<Vector4, ArgV2> const& v2) {
  return codi::ComputeExpression<Vector4, VectorAdditionOperation, ArgV1, ArgV2>(v1.cast(), v2.cast());
}
//! [Operation implementation]

// 3. Add member CoDiPack expression operations on the vector type.

//! [Member operation implementation]
// Operation for vector.norm().
template<typename T> // Template argument not used, we use the direct type specification of Vector4.
struct VectorNormOperation : public codi::UnaryJacobianOperation<Number, VectorNormOperation<T>> {

    /// Primal evaluation.
    static CODI_INLINE Number primal(Vector4 const& v) {
      return v.norm();
    }

    /// Gradient with respect to v.
    static CODI_INLINE Vector4 gradient(Vector4 const& v, Number const& result) {
      return v / result;
    }
};

/// Specialization of ExpressionMemberOperations for the vector type. This injects the member operations into all
/// expressions.
template<typename Impl>
struct codi::ExpressionMemberOperations<Vector4, Impl> {
  public:
    using Real = Vector4; ///< The expression value is the vector type.

    /// Returns the expression with the norm operation.
    auto norm() const {
      return codi::ComputeExpression<Number, VectorNormOperation, Impl>(this->cast());
    }

  private:
    CODI_INLINE Impl const& cast() const {
      return static_cast<Impl const&>(*this);
    }
};
//! [Member operation implementation]

// 4. Implement the ActiveVector4 with the help of the aggregated type class in CoDiPack.
//! [ActiveVector4 definition]
struct ActiveVector4 : public codi::AggregatedActiveType<Vector4, ActiveNumber, ActiveVector4> {

    using Base = codi::AggregatedActiveType<Vector4, ActiveNumber, ActiveVector4>; ///< Base abbreviation.

    using Base::Base; // Include constructors from the base class.
    /// Constructor with 4 arguments.
    template<typename Arg1, typename Arg2, typename Arg3>
    ActiveVector4(
        codi::ExpressionInterface<Number, Arg1> const& arg1,
        codi::ExpressionInterface<Number, Arg2> const& arg2,
        codi::ExpressionInterface<Number, Arg3> const& arg3,
        codi::ExpressionInterface<Number, Arg3> const& arg4) : Base() {
      Base::values[0] = arg1;
      Base::values[1] = arg2;
      Base::values[2] = arg3;
      Base::values[3] = arg4;
    }

    /// Array access operation.
    ActiveNumber& operator[](size_t index) {
      return Base::values[index];
    }
};
//! [ActiveVector4 definition]

//! [Function]
template<typename VectorType>
ActiveNumber func(ActiveNumber const& s1, VectorType const& v1, ActiveNumber const& s2, VectorType const& v2) {
  return (s1 * v1 + s2 * v2).norm();
}
//! [Function]

// Test the implementation with either ActiveVector4 or Vector4WithActiveType.
template<typename VectorType>
void test() {
  ActiveNumber::Tape& tape = ActiveNumber::getTape();

  VectorType v1(ActiveNumber(1.0), ActiveNumber(2.0), ActiveNumber(4.0), ActiveNumber(8.0));
  VectorType v2(ActiveNumber(0.1), ActiveNumber(0.2), ActiveNumber(0.4), ActiveNumber(0.8));

  ActiveNumber s1(5.0);
  ActiveNumber s2(0.5);

  tape.setActive();
  tape.registerInput(s1);
  tape.registerInput(s2);
  for(int i = 0; i < 4; i += 1) {
    tape.registerInput(v1[i]);
    tape.registerInput(v2[i]);
  }

  codi::TapeValues before = tape.getTapeValues();
  ActiveNumber res = func(s1, v1, s2, v2);
  codi::TapeValues diff = tape.getTapeValues().subtract(before);

  tape.registerOutput(res);

  tape.setPassive();

  res.gradient() = 1.0;

  tape.evaluate();

  std::cout << "d f/d s1 = " << s1.getGradient() << std::endl;
  std::cout << "d f/d s2 = " << s2.getGradient() << std::endl;
  for(int i = 0; i < 4; i += 1) {
    std::cout << "d f/d v1[" << i << "] = " << v1[i].getGradient() << std::endl;
  }
  for(int i = 0; i < 4; i += 1) {
    std::cout << "d f/d v2[" << i << "] = " << v2[i].getGradient() << std::endl;
  }

  diff.formatDefault();

  tape.reset();
}

int main() {
  using Vector4WithActiveType = Eigen::Matrix<ActiveNumber, 4, 1>;

  std::cout << "Running example with 'Vector4WithActiveType' vector type. No specialization are used for Eigen vector." << std::endl;
  test<Vector4WithActiveType>();

  std::cout << std::endl;
  std::cout << "Running example with 'ActiveVector4' vector type. The specializations are used for Eigen vector." << std::endl;
  test<ActiveVector4>();

  return 0;
}
//! [Tutorial 7 - Aggregated type implementation]
