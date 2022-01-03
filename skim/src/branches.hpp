/**************************************/
/*										*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TChain.h"
#include "constants.hpp"

int npart;
int evstat;
int intt;
int evntid;
int evtype;
int evntclas;
int evthel;
int evntclas2;
float q_l;
float t_l;
float tr_time;
float rf_time1;
float rf_time2;
int gpart;
int id[MAX_PARTS];    //[gpart]
int stat[MAX_PARTS];  //[gpart]
int dc[MAX_PARTS];    //[gpart]
int cc[MAX_PARTS];    //[gpart]
int sc[MAX_PARTS];    //[gpart]
int ec[MAX_PARTS];    //[gpart]
int lec[MAX_PARTS];   //[gpart]
int ccst[MAX_PARTS];  //[gpart]
float p[MAX_PARTS];   //[gpart]
float m[MAX_PARTS];   //[gpart]
int q[MAX_PARTS];     //[gpart]
float b[MAX_PARTS];   //[gpart]
float cx[MAX_PARTS];  //[gpart]
float cy[MAX_PARTS];  //[gpart]
float cz[MAX_PARTS];  //[gpart]
float vx[MAX_PARTS];  //[gpart]
float vy[MAX_PARTS];  //[gpart]
float vz[MAX_PARTS];  //[gpart]
int dc_part;
int dc_sect[MAX_PARTS];    //[dc_part]
int dc_trk[MAX_PARTS];     //[dc_part]
int dc_stat[MAX_PARTS];    //[dc_part]
float dc_vx[MAX_PARTS];    //[dc_part]
float dc_vy[MAX_PARTS];    //[dc_part]
float dc_vz[MAX_PARTS];    //[dc_part]
float dc_vr[MAX_PARTS];    //[dc_part]
float dc_xsc[MAX_PARTS];   //[dc_part]
float dc_ysc[MAX_PARTS];   //[dc_part]
float dc_zsc[MAX_PARTS];   //[dc_part]
float dc_cxsc[MAX_PARTS];  //[dc_part]
float dc_cysc[MAX_PARTS];  //[dc_part]
float dc_czsc[MAX_PARTS];  //[dc_part]
float dc_c2[MAX_PARTS];    //[dc_part]
int ec_part;
int ec_stat[MAX_PARTS];  //[ec_part]
int ec_sect[MAX_PARTS];  //[ec_part]
int ec_whol[MAX_PARTS];  //[ec_part]
int ec_inst[MAX_PARTS];  //[ec_part]
int ec_oust[MAX_PARTS];  //[ec_part]
float etot[MAX_PARTS];   //[ec_part]
float ec_ei[MAX_PARTS];  //[ec_part]
float ec_eo[MAX_PARTS];  //[ec_part]
float ec_t[MAX_PARTS];   //[ec_part]
float ec_r[MAX_PARTS];   //[ec_part]
float ech_x[MAX_PARTS];  //[ec_part]
float ech_y[MAX_PARTS];  //[ec_part]
float ech_z[MAX_PARTS];  //[ec_part]
float ec_m2[MAX_PARTS];  //[ec_part]
float ec_m3[MAX_PARTS];  //[ec_part]
float ec_m4[MAX_PARTS];  //[ec_part]
float ec_c2[MAX_PARTS];  //[ec_part]
int sc_part;
int sc_sect[MAX_PARTS];  //[sc_part]
int sc_hit[MAX_PARTS];   //[sc_part]
int sc_pd[MAX_PARTS];    //[sc_part]
int sc_stat[MAX_PARTS];  //[sc_part]
float edep[MAX_PARTS];   //[sc_part]
float sc_t[MAX_PARTS];   //[sc_part]
float sc_r[MAX_PARTS];   //[sc_part]
float sc_c2[MAX_PARTS];  //[sc_part]
int cc_part;
int cc_sect[MAX_PARTS];  //[cc_part]
int cc_hit[MAX_PARTS];   //[cc_part]
int cc_segm[MAX_PARTS];  //[cc_part]
int nphe[MAX_PARTS];     //[cc_part]
float cc_t[MAX_PARTS];   //[cc_part]
float cc_r[MAX_PARTS];   //[cc_part]
float cc_c2[MAX_PARTS];  //[cc_part]

void getBranches(TTree *myTree) {
  myTree->SetBranchStatus("*", 0);
  myTree->SetBranchAddress("npart", &npart);
  myTree->SetBranchAddress("evstat", &evstat);
  myTree->SetBranchAddress("intt", &intt);
  myTree->SetBranchAddress("evntid", &evntid);
  myTree->SetBranchAddress("evtype", &evtype);
  myTree->SetBranchAddress("evntclas", &evntclas);
  myTree->SetBranchAddress("evthel", &evthel);
  myTree->SetBranchAddress("evntclas2", &evntclas2);
  myTree->SetBranchAddress("q_l", &q_l);
  myTree->SetBranchAddress("t_l", &t_l);
  myTree->SetBranchAddress("tr_time", &tr_time);
  myTree->SetBranchAddress("rf_time1", &rf_time1);
  myTree->SetBranchAddress("rf_time2", &rf_time2);
  myTree->SetBranchAddress("gpart", &gpart);
  myTree->SetBranchAddress("id", id);
  myTree->SetBranchAddress("stat", stat);
  myTree->SetBranchAddress("dc", dc);
  myTree->SetBranchAddress("cc", cc);
  myTree->SetBranchAddress("sc", sc);
  myTree->SetBranchAddress("ec", ec);
  myTree->SetBranchAddress("lec", lec);
  myTree->SetBranchAddress("ccst", ccst);
  myTree->SetBranchAddress("p", p);
  myTree->SetBranchAddress("m", m);
  myTree->SetBranchAddress("q", q);
  myTree->SetBranchAddress("b", b);
  myTree->SetBranchAddress("cx", cx);
  myTree->SetBranchAddress("cy", cy);
  myTree->SetBranchAddress("cz", cz);
  myTree->SetBranchAddress("vx", vx);
  myTree->SetBranchAddress("vy", vy);
  myTree->SetBranchAddress("vz", vz);
  myTree->SetBranchAddress("dc_part", &dc_part);
  myTree->SetBranchAddress("dc_sect", dc_sect);
  myTree->SetBranchAddress("dc_trk", dc_trk);
  myTree->SetBranchAddress("dc_stat", dc_stat);
  myTree->SetBranchAddress("dc_vx", dc_vx);
  myTree->SetBranchAddress("dc_vy", dc_vy);
  myTree->SetBranchAddress("dc_vz", dc_vz);
  myTree->SetBranchAddress("dc_vr", dc_vr);
  myTree->SetBranchAddress("dc_xsc", dc_xsc);
  myTree->SetBranchAddress("dc_ysc", dc_ysc);
  myTree->SetBranchAddress("dc_zsc", dc_zsc);
  myTree->SetBranchAddress("dc_cxsc", dc_cxsc);
  myTree->SetBranchAddress("dc_cysc", dc_cysc);
  myTree->SetBranchAddress("dc_czsc", dc_czsc);
  myTree->SetBranchAddress("dc_c2", dc_c2);
  myTree->SetBranchAddress("ec_part", &ec_part);
  myTree->SetBranchAddress("ec_stat", ec_stat);
  myTree->SetBranchAddress("ec_sect", ec_sect);
  myTree->SetBranchAddress("ec_whol", ec_whol);
  myTree->SetBranchAddress("ec_inst", ec_inst);
  myTree->SetBranchAddress("ec_oust", ec_oust);
  myTree->SetBranchAddress("etot", etot);
  myTree->SetBranchAddress("ec_ei", ec_ei);
  myTree->SetBranchAddress("ec_eo", ec_eo);
  myTree->SetBranchAddress("ec_t", ec_t);
  myTree->SetBranchAddress("ec_r", ec_r);
  myTree->SetBranchAddress("ech_x", ech_x);
  myTree->SetBranchAddress("ech_y", ech_y);
  myTree->SetBranchAddress("ech_z", ech_z);
  myTree->SetBranchAddress("ec_m2", ec_m2);
  myTree->SetBranchAddress("ec_m3", ec_m3);
  myTree->SetBranchAddress("ec_m4", ec_m4);
  myTree->SetBranchAddress("ec_c2", ec_c2);
  myTree->SetBranchAddress("sc_part", &sc_part);
  myTree->SetBranchAddress("sc_sect", sc_sect);
  myTree->SetBranchAddress("sc_hit", sc_hit);
  myTree->SetBranchAddress("sc_pd", sc_pd);
  myTree->SetBranchAddress("sc_stat", sc_stat);
  myTree->SetBranchAddress("edep", edep);
  myTree->SetBranchAddress("sc_t", sc_t);
  myTree->SetBranchAddress("sc_r", sc_r);
  myTree->SetBranchAddress("sc_c2", sc_c2);
  myTree->SetBranchAddress("cc_part", &cc_part);
  myTree->SetBranchAddress("cc_sect", cc_sect);
  myTree->SetBranchAddress("cc_hit", cc_hit);
  myTree->SetBranchAddress("cc_segm", cc_segm);
  myTree->SetBranchAddress("nphe", nphe);
  myTree->SetBranchAddress("cc_t", cc_t);
  myTree->SetBranchAddress("cc_r", cc_r);
  myTree->SetBranchAddress("cc_c2", cc_c2);
}

#endif