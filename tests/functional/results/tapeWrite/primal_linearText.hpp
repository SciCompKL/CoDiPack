#include <codi.hpp>

template <typename Tape>
std::vector<typename Tape::EvalHandle> primal_linearTextCreateEvalHandles(){

  std::vector<typename Tape::EvalHandle> evalHandles;
  using Impl = codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> >;

  evalHandles.resize(4);
  evalHandles[0] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >>();
  evalHandles[1] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::ComputeExpression<double, codi::OperationSin, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > > >>();
  evalHandles[2] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::ComputeExpression<double, codi::OperationAdd, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > > >>();
  evalHandles[3] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::ComputeExpression<double, codi::OperationMultiply, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::ActiveType<codi::PrimalValueLinearTape<codi::PrimalValueTapeTypes<double, double, codi::LinearIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > > >>();

  return evalHandles;
}