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

#include "basic/testCopy.hpp"
#include "basic/testCopyHigherOrder.hpp"
#include "basic/testExpr.hpp"
#include "basic/testExprHigherOrder.hpp"
#include "basic/testIndices.hpp"
#include "basic/testOutput.hpp"
#include "exceptions/testOneArgumentExceptions.hpp"
#include "exceptions/testTwoArgumentExceptions.hpp"
#include "expressions/testAssignOperators1.hpp"
#include "expressions/testAssignOperators2.hpp"
#include "expressions/testBigExpressions.hpp"
#include "expressions/testIncrementOperators.hpp"
#include "expressions/testOneArgumentExpr1.hpp"
#include "expressions/testOneArgumentExpr2.hpp"
#include "expressions/testOneArgumentExpr3.hpp"
#include "expressions/testTwoArgumentExpr1.hpp"
#include "expressions/testTwoArgumentExpr2.hpp"
#include "externalFunctions/testExtFunctionCall.hpp"
#include "externalFunctions/testExtFunctionCallMultiple.hpp"
#include "io/testIO.hpp"
#include "io/testSwap.hpp"
#include "tools/helpers/testEigenLinearSystemSolverHandler.hpp"
#include "tools/helpers/testEigenSparseLinearSystemSolverHandler.hpp"
#include "tools/helpers/testEnzymeExternalFunctionHelper.hpp"
#include "tools/helpers/testExternalFunctionHelper.hpp"
#include "tools/helpers/testExternalFunctionHelperPassive.hpp"
#include "tools/helpers/testPreaccumulation.hpp"
#include "tools/helpers/testPreaccumulationForward.hpp"
#include "tools/helpers/testPreaccumulationForwardInvalidAdjoint.hpp"
#include "tools/helpers/testPreaccumulationLargeStatement.hpp"
#include "tools/helpers/testPreaccumulationPassiveValue.hpp"
#include "tools/helpers/testPreaccumulationZeroJacobi.hpp"
#include "tools/helpers/testReset.hpp"
#include "tools/helpers/testStatementPushHelper.hpp"
#include "tools/testReferenceActiveType.hpp"
#include "traits/testNumericLimits.hpp"
