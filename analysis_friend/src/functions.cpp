#include "functions.hpp"

double* fun::Vector_Array(std::vector<double> vec){
	double arr[vec.size()];
	for(int i=0; i<vec.size(); i++){
		arr[i] = vec.at(i);
	}
	return arr;
}

int* fun::Vector_Array(std::vector<int> vec){
	int arr[vec.size()];
	for(int i=0; i<vec.size(); i++){
		arr[i] = vec.at(i);
	}
	return arr;
}

bool fun::replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::shared_ptr<TFile> fun::Name_File(std::string file_name_){
    return std::make_shared<TFile>(file_name_.c_str(),"RECREATE");
}

double fun::nSparseIntegral(THnSparseD* nhist_)
{
    //std::cout<<"\tTaking an Integral\n";
  //Taken from most of the ComputeIntegral() for THnSparse, but got rid of some checks and the normalization
   // Calculate the integral of the histogram
  /*int dim = nhist_->GetNdimensions();
  int nbins = nhist_->GetNbins();
  double integral = NAN;
  for(int i=0; i<nbins; i++){
    integral+= nhist_->ComputeIntegral()[i];
  }
  return integral;
  */
   // check number of bins
   if (nhist_->GetNbins() == 0) {
      //Error("ComputeIntegral", "The histogram must have at least one bin.");
      return 0.;
   }
   // allocate integral array
   Double_t* fIntegral = new  Double_t [nhist_->GetNbins() + 1];
   fIntegral[0] = 0.;
   // fill integral array with contents of regular bins (non over/underflow)
   Int_t* coord = new Int_t[nhist_->GetNdimensions()];
   for (Long64_t i = 0; i < nhist_->GetNbins(); ++i) {
      Double_t v = nhist_->GetBinContent(i, coord);
      // check whether the bin is regular
      bool regularBin = true;
      for (Int_t dim = 0; dim < nhist_->GetNdimensions(); dim++){
         if (coord[dim] < 1 || coord[dim] > nhist_->GetAxis(dim)->GetNbins()) {
            regularBin = false;
            break;
         }
      }
      // if outlayer, count it with zero weight
      if (!regularBin) v = 0.;
      fIntegral[i + 1] = fIntegral[i] + v;
   }
   delete [] coord;
   // check sum of weights
   if (fIntegral[nhist_->GetNbins()] == 0.) {
      //Error("ComputeIntegral", "No hits in regular bins (non over/underflow).");
      delete [] fIntegral;
      return 0.;
   }
   // normalize the integral array
   /*for (Long64_t i = 0; i <= nhist_->GetNbins(); ++i){
      fIntegral[i] = fIntegral[i] / fIntegral[nhist_->GetNbins()];
   }*/

   // set status to valid
   //fIntegralStatus = kValidInt;
   return fIntegral[nhist_->GetNbins()];
}
