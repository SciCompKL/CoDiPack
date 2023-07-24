//! [Simple Tape]

#include <codi.hpp>

//! [Storing - Operator codes]
enum class OperatorCode {
  ADD,
  SUB,
  MUL,
  DIV,
  SIN,
  COS,
  COPY
};

struct OperatorCodeLookup {
  public:
    template<typename Op>
    static OperatorCode get();
};

template<typename Op>
OperatorCode OperatorCodeLookup::get() {
  std::cerr << "Missing specialization for operator code lookup." << std::endl;
  exit(-1);
}

#define SPECIALIZE_LOOKUP(NAME, CODE) \
  template<> \
  OperatorCode OperatorCodeLookup::get<codi:: NAME<double>>() { \
    return OperatorCode:: CODE; \
  }

SPECIALIZE_LOOKUP(OperationAdd, ADD)
SPECIALIZE_LOOKUP(OperationSubstract, SUB)
SPECIALIZE_LOOKUP(OperationMultiply, MUL)
SPECIALIZE_LOOKUP(OperationDivide, DIV)
SPECIALIZE_LOOKUP(OperationSin, SIN)
SPECIALIZE_LOOKUP(OperationCos, COS)

#undef SPECIALIZE_LOOKUP
//! [Storing - Operator codes]

struct SimpleTape : public codi::ReverseTapeInterface<double, double, int> {
  public:
    using Real = double;
    using Gradient = double;
    using Identifier = int;

//! [Data stream - Type definition]
    using OperatorData = codi::ChunkedData<codi::Chunk1<OperatorCode>>;
    using IdentifierData = codi::ChunkedData<codi::Chunk1<int>, OperatorData>;
    using PrimalData = codi::ChunkedData<codi::Chunk1<double>, IdentifierData>;
//! [Data stream - Type definition]

    using Position = PrimalData::Position;

  private:
    bool active;

//! [Identifiers - Member definition]
    std::vector<double> adjointVec;

    int maxIdentifier;
//! [Identifiers - Member definition]

//! [Data stream - Member definition]
    codi::EmptyData emptyData; // Required for the termination.
    OperatorData operatorData;
    IdentifierData identifierData;
    PrimalData primalData;
//! [Data stream - Member definition]

  public:

    SimpleTape() :
//! [Identifiers - Member initialization]
      active(false),
      adjointVec(1), // Reserve one for out of bounds gradient access.
      maxIdentifier(0),
//! [Identifiers - Member initialization]
//! [Data stream - Member creation]
      emptyData(),
      operatorData(1024),
      identifierData(1024),
      primalData(1024)
    {
      operatorData.setNested(&emptyData);
      identifierData.setNested(&operatorData);
      primalData.setNested(&identifierData);
    }
//! [Data stream - Member creation]

    /*******************************************************************************/
    // ReverseTapeInterface implementation

//! [Identifiers - Registration]
    template<typename Lhs> void registerInput(codi::LhsExpressionInterface<Real, Gradient, SimpleTape, Lhs>& value) {
      if (active) {
        value.cast().getIdentifier() = generateIdentifier();
      } else {
        value.cast().getIdentifier() = 0;
      }
    }


    template<typename Lhs> void registerOutput(codi::LhsExpressionInterface<Real, Gradient, SimpleTape, Lhs>& value) {
      // Do nothing, every identifier is unique.
    }
//! [Identifiers - Registration]

//! [Other - Activity]
    void setActive()      {active = true;}
    void setPassive()     {active = false;}
    bool isActive() const {return active;}
//! [Other - Activity]

//! [Evaluation - Entry]
    void evaluate() {
      primalData.evaluateReverse(primalData.getPosition(), primalData.getZeroPosition(), SimpleTape::evaluateStack, adjointVec);
    }
//! [Evaluation - Entry]

//! [Other - Misc]
    void clearAdjoints() {
      for (double& adj : adjointVec) {
        adj = 0.0;
      }
    }
    void reset(bool resetAdjoints = true) {
      if (resetAdjoints) {
        clearAdjoints();
      }

      maxIdentifier = 0;

      primalData.reset();
    }

    template<typename Stream = std::ostream> void printStatistics(Stream& out = std::cout) const {
      getTapeValues().formatDefault(out);
    }

    template<typename Stream = std::ostream> void printTableHeader(Stream& out = std::cout) const {
      getTapeValues().formatHeader(out);
    }

    template<typename Stream = std::ostream> void printTableRow(Stream& out = std::cout) const {
      getTapeValues().formatRow(out);
    }

    codi::TapeValues getTapeValues() const {
      codi::TapeValues values("Example tape");

      values.addSection("Adjoint vector");
      values.addLongEntry("Number of adjoints", (1 + maxIdentifier));
      values.addDoubleEntry("Memory allocated", sizeof(double) * (1 + maxIdentifier), true, true);

      values.addSection("Index manager");
      values.addLongEntry("Max. live indices", (1 + maxIdentifier));

      values.addSection("Primal data");
      primalData.addToTapeValues(values);
      values.addSection("Identifier data");
      identifierData.addToTapeValues(values);
      values.addSection("Operator data");
      operatorData.addToTapeValues(values);

      return values;
    }
//! [Other - Misc]

    /*******************************************************************************/
    // InternalStatementRecordingInterface implementation

    static bool constexpr AllowJacobianOptimization = false; // If certain operations can be hidden from the tape.

//! [Identifiers - Initialization]
    template<typename Real>
    void initIdentifier(Real& value, Identifier& identifier) {
      identifier = 0; // Initialize with zero we perform an online activity analysis.
    }

    template<typename Real>
    void destroyIdentifier(Real& value, Identifier& identifier) {
      // Do nothing: Identifiers are not reused.
    }
//! [Identifiers - Initialization]

//! [Storing - Entry]
    template<typename Lhs, typename Rhs>
    void store(Lhs& lhs, Rhs const& rhs) {
      if (active) {
        storeOperator(rhs, lhs.value(), lhs.getIdentifier(), true);
      } else {
        lhs.value() = rhs.getValue();
        lhs.getIdentifier() = 0;
      }
    }
//! [Storing - Entry]

    /*******************************************************************************/
    // GradientAccessTapeInterface implementation

//! [Adjoint - Access]
    void setGradient(Identifier const& identifier, Gradient const& grad,
                     AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
      gradient(identifier, adjointsManagement) = grad;
    }
    Gradient const& getGradient(Identifier const& identifier) const {
      return gradient(identifier);
    }

    Gradient& gradient(Identifier const& identifier,
                       AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
      if (AdjointsManagement::Automatic == adjointsManagement) {
        checkAndResizeAdjoints(identifier);
      }

      return adjointVec[identifier];
    }

    Gradient const& gradient(Identifier const& identifier,
                             AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const {
      if (AdjointsManagement::Automatic == adjointsManagement && identifier >= (int)adjointVec.size()) {
        return adjointVec[0];
      } else {
        return adjointVec[identifier];
      }
    }
//! [Adjoint - Access]

  private:

//! [Identifiers - Helper]
    void checkAndResizeAdjoints(int const& identifier) {
      if (identifier > maxIdentifier) {
        std::cerr << "Error: Tryinig to access an identifier which was not distributed." << std::endl;
      }
      if (identifier >= (int)adjointVec.size()) { // Only resize if necessary.
        adjointVec.resize(maxIdentifier + 1);
      }
    }

    int generateIdentifier() {
      maxIdentifier += 1;

      return maxIdentifier;
    }
//! [Identifiers - Helper]

//! [Storing - Helper class]
    template<typename Operation, typename = void>
    struct StoreOperator_Impl {
      public:
        template<typename Exp>
        static void store(
            Exp const& exp,
            SimpleTape& tape,
            double& resultValue,
            int& resultIdentifier,
            bool copy);
    };
//! [Storing - Helper class]

//! [Storing - Unary operator]
    template<typename Arg, template<typename> class Op>
    struct StoreOperator_Impl<codi::UnaryExpression<double, Arg, Op>> {
      public:
        static void store(
            codi::UnaryExpression<double, Arg, Op> const& exp, SimpleTape& tape,
            double& resultValue, int& resultIdentifier, bool copy)
        {

          double argValue;
          int argIdentifier;

          tape.storeOperator(exp.arg, argValue, argIdentifier, false);

          if (argIdentifier != 0) {
            // Active argument or branch => store the operator.
            tape.operatorData.reserveItems(1);
            tape.identifierData.reserveItems(2);
            tape.primalData.reserveItems(1);

            resultIdentifier = tape.generateIdentifier();

            tape.operatorData.pushData(OperatorCodeLookup::get<Op<double>>());
            tape.identifierData.pushData(argIdentifier);
            tape.identifierData.pushData(resultIdentifier);
            tape.primalData.pushData(argValue);
          } else {
            // Passive argument or branch => do not store anything.
            resultIdentifier = 0;
          }

          resultValue = exp.getValue();
        }
    };
//! [Storing - Unary operator]

//! [Storing - Other operators]
    template<typename Arg1, typename Arg2, template<typename> class Op>
    struct StoreOperator_Impl<codi::BinaryExpression<double, Arg1, Arg2, Op>> {
      public:
        static void store(
            codi::BinaryExpression<double, Arg1, Arg2, Op> const& exp, SimpleTape& tape,
            double& resultValue, int& resultIdentifier, bool copy)
        {

          double argAValue;
          double argBValue;
          int argAIdentifier;
          int argBIdentifier;

          tape.storeOperator(exp.argA, argAValue, argAIdentifier, false);
          tape.storeOperator(exp.argB, argBValue, argBIdentifier, false);

          if (argAIdentifier != 0 || argBIdentifier != 0) {
            // Active argument or branch => store the operator.
            tape.operatorData.reserveItems(1);
            tape.identifierData.reserveItems(3);
            tape.primalData.reserveItems(2);

            resultIdentifier = tape.generateIdentifier();

            tape.operatorData.pushData(OperatorCodeLookup::get<Op<double>>());
            tape.identifierData.pushData(argAIdentifier);
            tape.identifierData.pushData(argBIdentifier);
            tape.identifierData.pushData(resultIdentifier);
            tape.primalData.pushData(argAValue);
            tape.primalData.pushData(argBValue);

          } else {
            // Passive argument or branch => do not store anything.
            resultIdentifier = 0;
          }

          resultValue = exp.getValue();
        }
    };

    template<typename Exp>
    struct StoreOperator_Impl<Exp, codi::ExpressionTraits::EnableIfConstantExpression<Exp>> {
      public:
        static void store(
            codi::ConstantExpression<double> const& exp, SimpleTape& tape,
            double& resultValue, int& resultIdentifier, bool copy)
        {
          resultValue = exp.getValue();
          resultIdentifier = 0;
        }
    };

    template<typename Exp>
    struct StoreOperator_Impl<Exp, codi::ExpressionTraits::EnableIfLhsExpression<Exp>> {
      public:
        static void store(
            Exp const& exp, SimpleTape& tape,
            double& resultValue, int& resultIdentifier, bool copy)
        {
          if (copy && 0 != exp.getIdentifier()) {
            // Active argument and a copy operation => store the operator.
            tape.operatorData.reserveItems(1);
            tape.identifierData.reserveItems(2);
            tape.primalData.reserveItems(1);

            resultIdentifier = tape.generateIdentifier();

            tape.operatorData.pushData(OperatorCode::COPY);
            tape.identifierData.pushData(exp.getIdentifier());
            tape.identifierData.pushData(resultIdentifier);
            tape.primalData.pushData(exp.getValue());
          } else {
            // No copy operation or passive value => just pass the data.
            resultIdentifier = exp.getIdentifier();
          }

          resultValue = exp.getValue();
        }
    };
//! [Storing - Other operators]

    template<typename Exp>
    void storeOperator(Exp const& exp, double& value, int& identifier, bool copy) {
      StoreOperator_Impl<Exp>::store(exp, *this, value, identifier, copy);
    }

//![Evaluation - Stack]
    static void evaluateStack(
        /* data from call */
        std::vector<double>& adjointVector,
        /* primal data */
        size_t& curPrimalPos, size_t const& endPrimalPos, double const* const primalData,
        /* identifier data */
        size_t& curIdentifierPos, size_t const& endIdentifierPos, int const* const identifierData,
        /* operator data */
        size_t& curOperatorPos, size_t const& endOperatorPos, OperatorCode const* const operatorData) {

      while (curOperatorPos > endOperatorPos) {
        curOperatorPos -= 1;

        int resultIdentifier = 0;
        int arg1Identifier = 0;
        int arg2Identifier = 0;
        double arg1Value = 0.0;
        double arg2Value = 0.0;

        // Get the data from the stacks.
        switch (operatorData[curOperatorPos]) {
          // Binary operations.
          case OperatorCode::ADD:
          case OperatorCode::SUB:
          case OperatorCode::MUL:
          case OperatorCode::DIV:
            // Get the data.
            resultIdentifier = identifierData[curIdentifierPos - 1];
            arg2Identifier = identifierData[curIdentifierPos - 2];
            arg1Identifier = identifierData[curIdentifierPos - 3];
            arg2Value = primalData[curPrimalPos - 1];
            arg1Value = primalData[curPrimalPos - 2];

            // Adjust positions.
            curIdentifierPos -= 3;
            curPrimalPos -= 2;
            break;

          // Unary operations.
          case OperatorCode::COS:
          case OperatorCode::SIN:
          case OperatorCode::COPY:
            // Get the data
            resultIdentifier = identifierData[curIdentifierPos - 1];
            arg1Identifier = identifierData[curIdentifierPos - 2];
            arg1Value = primalData[curPrimalPos - 1];

            // Adjust positions.
            curIdentifierPos -= 2;
            curPrimalPos -= 1;
            break;
        default:
            std::cerr << "Error: unhandled pop operator '" << (int)operatorData[curOperatorPos] << "'." << std::endl;
            exit(-1);
          break;
        }

        double resultAdjoint = adjointVector[resultIdentifier];
        adjointVector[resultIdentifier] = 0.0; // Perform the reset of the lhs.
        double& arg1Adjoint = adjointVector[arg1Identifier];
        double& arg2Adjoint = adjointVector[arg2Identifier];

        switch (operatorData[curOperatorPos]) {
          // Binary operations
          case OperatorCode::ADD:
            arg1Adjoint += resultAdjoint;
            arg2Adjoint += resultAdjoint;
            break;
          case OperatorCode::SUB:
            arg1Adjoint += resultAdjoint;
            arg2Adjoint -= resultAdjoint;
            break;
          case OperatorCode::MUL:
            arg1Adjoint += arg2Value * resultAdjoint;
            arg2Adjoint += arg1Value * resultAdjoint;
            break;
          case OperatorCode::DIV:
            arg1Adjoint += resultAdjoint / arg2Value;
            arg2Adjoint += - arg1Value * resultAdjoint / arg2Value / arg2Value;
            break;

          // Unary operations
          case OperatorCode::COS:
            arg1Adjoint +=  -std::sin(arg1Value) * resultAdjoint;
            break;
          case OperatorCode::SIN:
            arg1Adjoint +=  std::cos(arg1Value) * resultAdjoint;
            break;
          case OperatorCode::COPY:
            arg1Adjoint +=  resultAdjoint;
            break;
        default:
            std::cerr << "Error: unhandled adjoint operator '" << (int)operatorData[curOperatorPos] << "'." << std::endl;
            exit(-1);
          break;
        }
      }
    }
//![Evaluation - Stack]
};

//![Example]
template<typename Real>
void eval() {
  using Tape = typename Real::Tape;
  Tape& tape = Real::getTape();

  Real a = 3.0;
  Real b = 4.0;

  tape.setActive();
  tape.registerInput(a);
  tape.registerInput(b);

  Real c = sin(a + b) * cos(a - b); // Single statement.

  tape.registerOutput(c);

  tape.setPassive();

  c.gradient() = 1.0;

  tape.evaluate();

  std::cout << "c = " << c.getValue() << std::endl;
  std::cout << "d c/d a = " << a.getGradient() << std::endl;
  std::cout << "d c/d b = " << b.getGradient() << std::endl;

  tape.printStatistics();

  tape.reset();
}

int main(int nargs, char** args) {

  std::cout << "Simple tape:" << std::endl;
  eval<codi::ActiveType<SimpleTape>>();

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "codi::RealReverse:" << std::endl;
  eval<codi::RealReverse>();


  return 0;
}
//![Example]
//! [Simple Tape]
