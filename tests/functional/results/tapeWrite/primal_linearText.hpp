#include <codi.hpp>

template <typename Tape>
std::vector<typename Tape::EvalHandle> primal_linearTextCreateEvalHandles(){

  std::vector<typename Tape::EvalHandle> evalHandles;
  using Impl = codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> >;

  evalHandles.resize(4);
  evalHandles[0] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >>();
  evalHandles[1] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::UnaryExpression<double, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::OperationSin>>();
  evalHandles[2] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::BinaryExpression<double, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::OperationAdd>>();
  evalHandles[3] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::BinaryExpression<double, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::OperationMultiply>>();

  return evalHandles;
}