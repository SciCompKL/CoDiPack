#include <codi.hpp>

template <typename Tape>
std::vector<typename Tape::EvalHandle> primal_multiuseBinaryCreateEvalHandles(){

  std::vector<typename Tape::EvalHandle> evalHandles;
  using Impl = codi::PrimalValueReuseTape<codi::PrimalValueTapeTypes<double, double, codi::MultiUseIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> >;

  evalHandles.resize(4);
  evalHandles[0] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::UnaryExpression<double, codi::ActiveType<codi::PrimalValueReuseTape<codi::PrimalValueTapeTypes<double, double, codi::MultiUseIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::OperationSin>>();
  evalHandles[1] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::BinaryExpression<double, codi::ActiveType<codi::PrimalValueReuseTape<codi::PrimalValueTapeTypes<double, double, codi::MultiUseIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::ActiveType<codi::PrimalValueReuseTape<codi::PrimalValueTapeTypes<double, double, codi::MultiUseIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::OperationAdd>>();
  evalHandles[2] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::BinaryExpression<double, codi::ActiveType<codi::PrimalValueReuseTape<codi::PrimalValueTapeTypes<double, double, codi::MultiUseIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::ActiveType<codi::PrimalValueReuseTape<codi::PrimalValueTapeTypes<double, double, codi::MultiUseIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >, codi::OperationMultiply>>();
  evalHandles[3] = Tape::StatementEvaluator::template createHandle<Impl, Impl, codi::ActiveType<codi::PrimalValueReuseTape<codi::PrimalValueTapeTypes<double, double, codi::MultiUseIndexManager<int>, codi::InnerStatementEvaluator, codi::DefaultChunkedData> > >>();

  return evalHandles;
}