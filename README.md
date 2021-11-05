# analysis_twopi_clas6
 **PhD Analysis of Two Pion Electroproduction off the Proton with CLAS**
 Utilizing Datasets E1-6 and E1F

## How to Run
 This analysis is separated into various steps, which I'll go through each individuall
### Event Selection
 This is located in the main directory and all the src files refer to it
 The purpose of this program is to read in root files from CLAS6 runs
 Inputs are given in the form of flags
#### Event Selection Flags
 Flags are handled flags.(c/h)pp and Class "Flags"
 Information about what is to be performed is put in several different channels
 **target information (-t=)**
    e16 => Looking at data from run group E1-6
    e1f => Looking at data from run group E1F
    filled => Target was filled
    sim => This is simulated Data
    hel => Helicity and Half Wave Plate status will be utilized
 **Data Files (-loc=)**
 Data is expected to be in root files with a separate text file having all their full paths separted by line. The pull path of this text file should be input as: -loc=<Full Path to Path File>
 **Cuts (-cut= || -ncut=)**
 all
 {particle} => {ele,pro,pip,pim}
 {topology} => {pro,pip,pim,zero}
    fid_{particle}
    dt_{particle}
    sf
    cc
    ec
    id
    beta_{particle}
    vertex
    mm_{topology}
 **Efficiency Cuts (-eff=)** (not usable yet)
    cc
    dc
    ec
    sc
 **Corrections (-corr=)** 
    p_corr  => Momentum Correction
    e_corr  => Energy Correction
    v_corr  => Vertex Correction
 **Plots (-plot= || -nplot=)**
    all
    wq2
    fid_{particle}
    sf
    cc
    dt_{particle}
    cc_eff
    sc_eff
    ec_eff
    dc_eff
    cc_geo
    sc_geo
    ec_geo
    mm_{topology}
    beta_{particle}
    vertex
 **THnSparse Friend (-friend=)**
 If wanting to output a THnSparse full of events, simply name the full path of the intended friend
 **Image (-image=)** (not usable yet)
  If wanting to output image files of histograms, then name the full path of the intended image file
 **Output file (-name=)**
 Name the full path of the intended output root file with event selection histograms in it
 **Other**
 -n=  => number of files to be used. use -1 if you want all of them in the given path file
 -cores=  => Used for multithreading. Input the number of threads you wish to utilize
 **Example Running**
 ./bin/analysis -t=e16 -t=filled -loc=/folder/path.txt -cut=dt_pro -cut=fid_ele -cut=sf -cut=mm_pim -plot=all -nplot=cc_eff -nplot=mm_pro -name=/outdir/name.root -friend=/outdir/name_friend.root -n=30 -cores=6
### Cross Section Extraction
 This takes in THnSparse histograms as output by the Event Selection above and outputs relevant cross sections and other yields
#### Cross Section Extraction Flags
 **Topologies (-t= and -r=)**
 {topologies} => {mpro,mpip,mpim,mzero,mall}
 -t= determines which THnSparse is being extracted based on topology. This is usually "mall"
 -r= determines what the actual topology given is. Within the context of the analysis, there may be times when it really is a missing pi-, but we want to use the "mall" file
 **Histograms (-h=)**
    all
    pol => Polarization Observables
    single_diff => Single Differential Cross Sections
    beam_spin => Beam Spin Asymmetry
    eff => Efficiencies
    err => Errors
    wq2 => W Q2 distributions
    accept => Acceptance
 **Output (-name=)**
 Full path of the output file
 **Files to put in**
 Give full path for these input files
    -sim=
    -exp=
    -weight=
    -empty=
    -sim2=
    -exp2=
    -weight2=
    -empty2=
 the 2s are used for running over both e16 and e1f data sets simultaneously for a combined measured cross section
 **Run Information (-i=)**
    e16
    e1f
    both
    var_pim => Specific variable set where pim is a focused product hadron
    var_pro => Specific variable set where pro is a focused product hadron
    var_pip => Specific variable set where pip is a focused product hadron
 **Image File (-image=)** (not used yet)
 Include full path for intended image file to output images of all plots made. 
### Golden Run
 Used to determine a golden run list based on Integrated Faraday Cup Charge normalized by number of Events 
 **Example Run**
 all files are full paths
 ./golden_run <data path file> <number of files> <output name> <FRONT> <MID> <LOW> <TOP> <RUN>
 **Explanation of some parameters**
 FRONT
    The front string of a data file before getting to its run number. For example if file is usually:
    /directory/e16_run_34950_pass2.a.07.root where that 34950 is the run number, then <FRONT> would be /directory/e16_run_
 MID
    The string between the run number and the iteration within the run. Using the previous example, MID =_pass2.a.
 TOP/LOW
    These are float numbers that give the acceptable bounds through which to cut on the normalized Integrated Faraday Cup Charge. This is only useful when outputting the list, but needs to utilize fitting programs elsewhere to themselves be determined. Until determined they can be whatever you wish
 RUN
    This is either "e16" or "e1f"