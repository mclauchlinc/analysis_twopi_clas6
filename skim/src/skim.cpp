/**************************************/
/*										*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "branches.hpp"
#include "skim.hpp"

Skim::Skim(std::string input, std::string output) {
  fin = input;
  fout = output;

  if (getenv("GPART") != NULL) {
    GPART = atof(getenv("GPART"));
    std::cout << "Cutting: gpart >= " << GPART << std::endl;
  }

  chain = new TChain("h10");
  chain->AddFile(fin.c_str());
  getBranches(chain);

  RootOutputFile = new TFile(fout.c_str(), "RECREATE");
  RootOutputFile->SetCompressionLevel(5);
}
Skim::~Skim() {}

void Skim::Basic() {
  bool electron_cuts;
  std::cout << BOLDWHITE << "Skimming file " << BOLDGREEN << fin << BOLDBLUE << " ==> " << BOLDGREEN << fout << RESET
            << std::endl;

  int num_of_events = (int)chain->GetEntries();

  // Make a copy of the chain we created earlier
  TTree *skim = chain->CloneTree(0);

  // Loop through all the events
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    electron_cuts = true;

    // Charge of 0th particle is negative
    electron_cuts &= (q[0] == -1);
    // The event has more than GPART particles
    electron_cuts &= (gpart >= GPART);
    // 0th particle hit EC
    electron_cuts &= (ec[0] > 0);
    // 0th particle hit CC
    electron_cuts &= (cc[0] > 0);
    // 0th particle hit SC
    electron_cuts &= (sc[0] > 0);
    // 0th particle hit DC
    electron_cuts &= (dc[0] > 0);

    if (electron_cuts) {
      // If the event passes all the cuts then fill the skim branch
      skim->Fill();
    }
  }

  // delete Chain object
  chain->Reset();
  delete chain;

  // CD into output file and write
  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
  delete RootOutputFile;
}