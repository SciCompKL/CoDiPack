/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#include <toolDefines.h>

#include <vector>

IN(2)
OUT(2)
POINTS(1) =
{
  {  1.0,     0.5}
};


void evalFunc(NUMBER* x, NUMBER* y) {
  y[0] = x[0];
  y[1] = x[1];
  for(int i = 0; i < 5; ++i) {
   NUMBER xTemp = y[0];
   NUMBER yTemp = y[1];

   y[0] = xTemp * xTemp - yTemp * yTemp - 0.65;
   y[1] = 2.0 * yTemp * xTemp;
  }
}

void func(NUMBER* x, NUMBER* y) {

  std::vector<NUMBER::GradientData> inputData;
  std::vector<NUMBER*> outputData;
  inputData.push_back(x[0].getGradientData());
  inputData.push_back(x[1].getGradientData());

  NUMBER::TapeType& tape = NUMBER::getGlobalTape();
  NUMBER::TapeType::Position startPos = tape.getPosition();

  evalFunc(x, y);

  outputData.push_back(&y[0]);
  outputData.push_back(&y[1]);

  NUMBER::TapeType::Position endPos = tape.getPosition();

  unsigned short nVarIn  = inputData.size();
  unsigned short nVarOut = outputData.size();
	double* jacobi     = new double[nVarOut*nVarIn];
	unsigned short* nNonzero        = new unsigned short[nVarOut];

	for (unsigned short iVarOut = 0; iVarOut < nVarOut; iVarOut++) {
		nNonzero[iVarOut] = 0;
		NUMBER::GradientData index_out = outputData[iVarOut]->getGradientData();

		tape.setGradient(index_out, 1.0);
		tape.evaluate(endPos, startPos);

		for (unsigned short iVarIn= 0; iVarIn < nVarIn; iVarIn++) {
			NUMBER::GradientData index_in =  inputData[iVarIn];
			jacobi[iVarOut*nVarIn+iVarIn] = tape.getGradient(index_in);
			if (jacobi[iVarOut*nVarIn+iVarIn] != 0.0) {
				nNonzero[iVarOut]++;
			}
			tape.setGradient(index_in, 0.0);
		}
		tape.setGradient(index_out, 0.0);
		tape.clearAdjoints(endPos, startPos);
	}

	if (nVarOut > 0) {
		tape.reset(startPos);
	}

	for (unsigned short iVarOut = 0; iVarOut < nVarOut; iVarOut++) {
		if (nNonzero[iVarOut] != 0){
			tape.store(outputData[iVarOut]->getValue(), outputData[iVarOut]->getGradientData(), nNonzero[iVarOut]);
			for (unsigned short iVarIn = 0; iVarIn < nVarIn; iVarIn++) {
				NUMBER::GradientData index_in =  inputData[iVarIn];
				tape.pushJacobi(jacobi[iVarOut*nVarIn+iVarIn], jacobi[iVarOut*nVarIn+iVarIn], 0.0, index_in);
			}
		}
	}
}
