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


THnSparseD* fun::Add_THnSparse(THnSparseD* hist1_, THnSparseD* hist2_, int sign_, std::vector<int> num_bins_){
   long idx = 0;
   std::vector<long> space_dims;
   for(int i=0; i<num_bins_.size(); i++){
      space_dims.push_back(num_bins_[i]);
   }
   CartesianGenerator cart(space_dims);
   int bin[num_bins_.size()];
   THnSparseD* output = (THnSparseD*)hist1_->Clone();
   while(cart.GetNextCombination()){
      for(int j = 0; j<num_bins_.size(); j++){
         bin[j] = cart[j];
      }
      output->AddBinContent(bin,hist2_->GetBinContent(bin)*sign_);
   }
   return output; 
}

THnSparseD* fun::Localized_Holes(THnSparseD* exp_hist_, THnSparseD* sim_hist_, THnSparseD* sim_hole_hist_, std::vector<int> num_bins_){
   long idx = 0;
   std::vector<long> space_dims;
   for(int i=0; i<num_bins_.size(); i++){
      space_dims.push_back(num_bins_[i]);
   }
   CartesianGenerator cart(space_dims);
   int bin[num_bins_.size()];
   int bin2[num_bins_.size()];
   THnSparseD* output = (THnSparseD*)sim_hole_hist_->Clone();
   double top = NAN; 
   double bot = NAN;
   double prev_val = NAN;
   bool look_further = true;
   int dist = 0; 
   std::vector<std::vector<int>> surr_bins; 
   while(cart.GetNextCombination()){
      for(int j = 0; j<num_bins_.size(); j++){
         bin[j] = cart[j];
      }
      while(look_further){
         dist++;
         surr_bins = fun::Surrounding_Bin(bin,dist,num_bins_);
         for(int i=0; i<surr_bins.size(); i++){
            for(int k=0; k<num_bins_.size(); k++){
               bin2[k] = surr_bins[i][k];
            }
            top += exp_hist_->GetBinContent(bin2);
            bot += sim_hist_->GetBinContent(bin2);
         }
         if(top!=0 && bot!=0){
            look_further=false;
            prev_val = sim_hole_hist_->GetBinContent(bin);
            output->SetBinContent(bin,(top/bot)*prev_val);
         }
         surr_bins.clear();
      }
   }
   return output; 
}

std::vector<std::vector<int>> fun::Surrounding_Bin(int* bin_, int dist_, std::vector<int> num_bins_){
   std::vector<std::vector<int>> surr_bin;
   std::vector<int> surr_bin1; 
   std::vector<long> space_dims;
   for(int i=0; i<num_bins_.size(); i++){
      space_dims.push_back(1+2*dist_);
   }
   CartesianGenerator cart(space_dims);
   int bin[num_bins_.size()];
   int bin2[num_bins_.size()];
   //int not_the_bin = 0;
   int inside_dist_range = 0;
   bool in_range = true;
   while(cart.GetNextCombination()){
      for(int j=0; j<dist_; j++){
         //not_the_bin = 0;
         inside_dist_range = 0;
         for(int i= 0; i<num_bins_.size(); i++){
            surr_bin1.push_back(bin_[i] + (-(j+1)+cart[i]));
            //if((-(j+1)+cart[i])==0){//Not counting the bin itself
            //   not_the_bin++;
            //}
            if(abs(-(j+1)+cart[i])<(j+1)){//Eliminate double counting from previous bouts
               inside_dist_range++;
            }
            if(surr_bin1[i] >= num_bins_[i]){//Stay within the bins that are actually accessible
               in_range &=false;
            }         
         }
         if(inside_dist_range<num_bins_.size() && in_range){//&&not_the_bin<num_bins_.size() 
            surr_bin.push_back(surr_bin1);
            surr_bin1.clear();
         }
      }
   }
   return surr_bin;
}
