# analysis_twopi_clas6
 <p>**PhD Analysis of Two Pion Electroproduction off the Proton with CLAS**<br>
 Utilizing Datasets E1-6 and E1F</p>

## How to Run
 This analysis is separated into various steps, which I'll go through each individually
### Event Selection
 <p>This is located in the main directory and all the src files refer to it<br>
 The purpose of this program is to read in root files from CLAS6 runs<br>
 Inputs are given in the form of flags</p>

#### Event Selection Flags
 <p>Flags are handled flags.(c/h)pp and Class "Flags"<br>
 Information about what is to be performed is put in several different channels<br>
 **target information (-t=)**<br>
 <ul>   
    <li>e16 => Looking at data from run group E1-6</li>
    <li>e1f => Looking at data from run group E1F</li>
    <li>filled => Target was filled</li>
    <li>sim => This is simulated Data</li>
    <li>hel => Helicity and Half Wave Plate status will be utilized</li>
 </ul>
 **Data Files (-loc=)**<br>
 Data is expected to be in root files with a separate text file having all their full paths separted by line. The pull path of this text file should be input as: -loc=<Full Path to Path File><br>
 **Cuts (-cut= || -ncut=)**<br>
 {particle} => {ele,pro,pip,pim}<br>
 {topology} => {pro,pip,pim,zero}<br>
 <ul>
    <li>all</li>
    <li>fid_{particle}</li>
    <li>dt_{particle}</li>
    <li>sf</li>
    <li>cc</li>
    <li>ec</li>
    <li>id</li>
    <li>beta_{particle}</li>
    <li>vertex</li>
    <li>mm_{topology}</li>
 </ul>
 **Efficiency Cuts (-eff=)** (not usable yet)<br>
 <ul>
    <li>cc</li>
    <li>dc</li>
    <li>ec</li>
    <li>sc</li>
 </ul>
 **Corrections (-corr=)** <br>
 <ul>   
    <li>p_corr  => Momentum Correction</li>
    <li>e_corr  => Energy Correction</li>
    <li>v_corr  => Vertex Correction</li>
 </ul>
 **Plots (-plot= || -nplot=)**<br>
 <ul>   
    <li>all</li>
    <li>wq2</li>
    <li>fid_{particle}</li>
    <li>sf</li>
    <li>cc</li>
    <li>dt_{particle}</li>
    <li>cc_eff</li>
    <li>sc_eff</li>
    <li>ec_eff</li>
    <li>dc_eff</li>
    <li>cc_geo</li>
    <li>sc_geo</li>
    <li>ec_geo</li>
    <li>mm_{topology}</li>
    <li>beta_{particle}</li>
    <li>vertex</li>
 </ul>
 **THnSparse Friend (-friend=)**<br>
 If wanting to output a THnSparse full of events, simply name the full path of the intended friend<br>
 **Image (-image=)** (not usable yet)<br>
  If wanting to output image files of histograms, then name the full path of the intended image file<br>
 **Output file (-name=)**<br>
 Name the full path of the intended output root file with event selection histograms in it<br>
 **Other**<br>
 <ul>
    <li>-n=  => number of files to be used. use -1 if you want all of them in the given path file</li>
    <li>-cores=  => Used for multithreading. Input the number of threads you wish to utilize</li>
 </ul>
 **Example Running**<br>
 ./bin/analysis -t=e16 -t=filled -loc=/folder/path.txt -cut=dt_pro -cut=fid_ele -cut=sf -cut=mm_pim -plot=all -nplot=cc_eff -nplot=mm_pro -name=/outdir/name.root -friend=/outdir/name_friend.root -n=30 -cores=6</p>

### Cross Section Extraction
 <p>This takes in THnSparse histograms as output by the Event Selection above and outputs relevant cross sections and other yields</p>

#### Cross Section Extraction Flags
 <p>**Topologies (-t= and -r=)**<br>
 {topologies} => {mpro,mpip,mpim,mzero,mall}<br>
 <ul>
    <li>-t= determines which THnSparse is being extracted based on topology. This is usually "mall"</li>
    <li>-r= determines what the actual topology given is. Within the context of the analysis, there may be times when it really is a missing pi-, but we want to use the "mall" file</li>
 </ul>
 **Histograms (-h=)**<br>
  <ul>
    <li>all</li>
    <li>pol => Polarization Observables</li>
    <li>single_diff => Single Differential Cross Sections</li>
    <li>beam_spin => Beam Spin Asymmetry</li>
    <li>eff => Efficiencies</li>
    <li>err => Errors</li>
    <li>wq2 => W Q2 distributions</li>
    <li>accept => Acceptance</li>
 </ul>
 **Output (-name=)**<br>
 Full path of the output file<br>
 **Files to put in**<br>
 Give full path for these input files<br>
 <ul>
    <li>-sim=</li>
    <li>-exp=</li>
    <li>-weight=</li>
    <li>-empty=</li>
    <li>-sim2=</li>
    <li>-exp2=</li>
    <li>-weight2=</li>
    <li>-empty2=</li>
 </ul>
 the 2s are used for running over both e16 and e1f data sets simultaneously for a combined measured cross section
 **Run Information (-i=)**<br>
 <ul>
    <li>e16</li>
    <li>e1f</li>
    <li>both</li>
   <li> var_pim => Specific variable set where pim is a focused product hadron</li>
    <li>var_pro => Specific variable set where pro is a focused product hadron</li>
    <li>var_pip => Specific variable set where pip is a focused product hadron</li>
 </ul>
 **Image File (-image=)** (not used yet)<br>
 Include full path for intended image file to output images of all plots made. </p>

### Golden Run
 <p>Used to determine a golden run list based on Integrated Faraday Cup Charge normalized by number of Events <br>
 **Example Run**<br>
 all files are full paths<br>
 ./golden_run <data path file> <number of files> <output name> <FRONT> <MID> <LOW> <TOP> <RUN><br>
 **Explanation of some parameters**
 FRONT<br>
    The front string of a data file before getting to its run number. For example if file is usually:<br>
    /directory/e16_run_34950_pass2.a.07.root where that 34950 is the run number, then <FRONT> would be /directory/e16_run_<br>
 MID<br>
    The string between the run number and the iteration within the run. Using the previous example, MID =_pass2.a.<br>
 TOP/LOW<br>
    These are float numbers that give the acceptable bounds through which to cut on the normalized Integrated Faraday Cup Charge. This is only useful when outputting the list, but needs to utilize fitting programs elsewhere to themselves be determined. Until determined they can be whatever you wish<br>
 RUN<br>
    This is either "e16" or "e1f"</p>