#include "TChain.h"
#include "TStopwatch.h"
#include "skim.hpp"

using namespace std;

int main(int argc, char **argv) {
  std::string infile;
  std::string outfile;

  if (argc == 2) {
    infile = argv[1];
    outfile = infile.substr(0, infile.size() - 5) + "_skim.root";
  } else if (argc == 3) {
    infile = argv[1];
    outfile = argv[2];
  } else {
    std::cerr << BOLDRED << "Error: \n";
    std::cerr << BOLDWHITE << "\tNeed input file and output file\n";
    std::cerr << BOLDBLUE << "Usage:\n\t";
    std::cerr << RESET << argv[0] << " /path/to/unskimmed.root\n";
    std::cerr << BOLDBLUE << "or:\n\t";
    std::cerr << RESET << argv[0] << " /path/to/unskimmed.root /path/to/skimmed.root\n";
    std::cerr << std::endl;
    return 1;
  }

  Skim *s = new Skim(infile, outfile);
  s->Basic();

  return 0;
}