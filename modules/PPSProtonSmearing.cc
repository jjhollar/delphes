/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class PPSProtonSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PPSProtonSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

PPSProtonSmearing::PPSProtonSmearing() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

PPSProtonSmearing::~PPSProtonSmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void PPSProtonSmearing::Init()
{
    fXiSmear = GetDouble("XiResolution", 0.05);
    fXiMin = GetDouble("XiMin",0.02);
    fXiMax = GetDouble("XiMax",0.20);
    fTimeRes = GetDouble("TimeRes",0.1);
    
  // read resolution formula
  
  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ProtonFilter/protons"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void PPSProtonSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void PPSProtonSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t z, e, pz, t, pt, eta, phi, res, enew, znew, xi, xinew;
  Double_t ttrue, tnew;

  const Double_t c_light = 2.99792458E8;
  
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    z = candidatePosition.Z();
    e = candidateMomentum.E();
    pz = candidateMomentum.Pz();
    
    t = candidatePosition.T();
    pt = candidateMomentum.Pt();
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();

    //    res = fFormula->Eval(pt, eta, phi, e, candidate);
    xi = 1.0 - (e/6800.0);
    xinew = gRandom->Gaus(xi,fXiSmear*xi);

    if(pz>0)
      ttrue = std::abs((21650.0 - (z/10))/30.0);
    if(pz<0)
      ttrue = std::abs((-21650.0 - (z/10))/30.0);

    //    std::cout << "ttrue = " << ttrue << " (pz = " << pz << ", z = " << z/10 << "), ";
    
    // apply smearing formula
    //pt = gRandom->Gaus(pt, fFormula->Eval(pt, eta, phi, e) * pt);
    enew = 6800.0*(1.0 -xinew);
    //    znew = gRandom->Gaus(z,0.5);
    tnew = gRandom->Gaus(ttrue,fTimeRes);

    //    std::cout << "tsmear = " << tnew << std::endl;
    
    //    res = (res > 1.0) ? 1.0 : res;

    //if(pt <= 0.0) continue;

    if(e > 3000 && (xinew>=fXiMin) && (xinew<=fXiMax))
      {
	candidate->Momentum.SetPtEtaPhiE(pt, eta, phi, enew);
	candidate->Position.SetT(tnew * 1.0E3 * c_light);	
	//	candidate->Position.SetXYZT(0,0,z,tnew);
	//	std::cout << "\ttnew = " << tnew << ", Candidate T = " << candidate->Position.T() << std::endl;
	fOutputArray->Add(candidate);
      }
  }
}
//----------------------------------------------------------------

//------------------------------------------------------------------------------
