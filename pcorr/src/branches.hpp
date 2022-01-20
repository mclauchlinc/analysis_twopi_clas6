
#ifndef BRANCHES_HPP
#define BRANCHES_HPP
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "TChain.h"
#include "constants.hpp"


//Create a data chain for a given data set
//Optional number of files to load

class Branches {
 private:
  std::shared_ptr<TChain> _tree;
  bool _MC = false;
  bool _include_weights = false; 
  int _npart = 0;
  int _evstat = 0;
  int _intt = 0;
  int _evntid = 0;
  int _evtype = 0;
  int _evntclas= 0;
  int _evthel= 0;
  int _evntclas2 = 0;
  float _q_l = NAN;
  float _t_l= NAN;
  float _tr_time= NAN;
  float _rf_time1= NAN;
  float _rf_time2= NAN;
  int _gpart = 0;
  int _id[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[gpart]
  int _stat[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[gpart]
  int _dc[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[gpart]
  int _cc[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[gpart]
  int _sc[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[gpart]
  int _ec[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[gpart]
  int _lec[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};   //[gpart]
  int _ccst[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[gpart]
  float _p[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[gpart]
  float _m[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[gpart]
  int _q[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};     //[gpart]
  float _b[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[gpart]
  float _cx[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[gpart]
  float _cy[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[gpart]
  float _cz[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[gpart]
  float _vx[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[gpart]
  float _vy[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[gpart]
  float _vz[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[gpart]
  int _dc_part = 0;
  int _dc_sect[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[dc_part]
  int _dc_trk[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};     //[dc_part]
  int _dc_stat[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[dc_part]
  float _dc_vx[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[dc_part]
  float _dc_vy[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[dc_part]
  float _dc_vz[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[dc_part]
  float _dc_vr[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[dc_part]
  float _dc_xsc[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[dc_part]//x y z of SC as determined by DC
  float _dc_ysc[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[dc_part]
  float _dc_zsc[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[dc_part]
  float _dc_cxsc[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[dc_part]
  float _dc_cysc[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[dc_part]
  float _dc_czsc[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[dc_part]
  float _dc_c2[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[dc_part]
  int _ec_part = 0;
  int _ec_stat[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[ec_part]
  int _ec_sect[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[ec_part]
  int _ec_whol[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[ec_part]
  int _ec_inst[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[ec_part]
  int _ec_oust[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[ec_part]
  float _etot[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[ec_part]
  float _ec_ei[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  float _ec_eo[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  float _ec_t[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[ec_part]
  float _ec_r[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[ec_part]
  float _ech_x[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part] //x y z of EC determined from cross hatched hits of different layers
  float _ech_y[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  float _ech_z[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  float _ec_m2[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  float _ec_m3[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  float _ec_m4[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  float _ec_c2[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[ec_part]
  int _sc_part = 0;
  int _sc_sect[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[sc_part]
  int _sc_hit[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};   //[sc_part]
  int _sc_pd[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[sc_part]
  int _sc_stat[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[sc_part]
  float _edep[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};   //[sc_part]
  float _sc_t[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};   //[sc_part]
  float _sc_r[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};   //[sc_part]
  float _sc_c2[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[sc_part]
  int _cc_part = 0;
  int _cc_sect[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[cc_part]
  int _cc_hit[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};   //[cc_part]
  int _cc_segm[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[cc_part]
  int _nphe[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};     //[cc_part]
  float _cc_t[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[cc_part]
  float _cc_r[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[cc_part]// distance from vertex for hit in CC
  float _cc_c2[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};  //[cc_part]// Cos squared for angle out towards CC
  ////////////
  int _nprt = 0;
  int _pidpart[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    //[nprt]
  float _xpart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[nprt]
  float _ypart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[nprt]
  float _zpart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[nprt]
  float _epart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[nprt]
  float _pxpart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[nprt]
  float _pypart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[nprt]
  float _pzpart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};   //[nprt]
  float _qpart[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};    //[nprt]
  int _flagspart[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //[nprt]

  float _mcp[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
  float _mctheta[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
  float _mcphi[MAX_PARTS] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
  int _mcid[MAX_PARTS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  float _weight = NAN;

 public:
  Branches(std::shared_ptr<TChain> tree, bool sim_);
  ~Branches(){};
  bool MC(); 
  void init();
  int npart();
  int evstat();
  int intt();
  int evntid();
  int evtype();
  int evntclas();
  int evthel();
  int evntclas2();
  float q_l();
  float t_l();
  float tr_time();
  float rf_time1();
  float rf_time2();
  int gpart();
  int dc_part();
  int ec_part();
  int sc_part();
  int cc_part();
  int lac_part();
  int st_part();
  int nprt();

  int id(int i);
  int stat(int i);
  int dc(int i);
  int cc(int i);
  int sc(int i);
  int ec(int i);
  int lec(int i);
  int ccst(int i);
  float p(int i);
  float px(int i);
  float py(int i);
  float pz(int i);
  float m(int i);
  int q(int i);
  float b(int i);
  float cx(int i);
  float cy(int i);
  float cz(int i);
  float vx(int i);
  float vy(int i);
  float vz(int i);
  int dc_sect(int i);
  int dc_trk(int i);
  int dc_stat(int i);
  float dc_vx(int i);
  float dc_vy(int i);
  float dc_vz(int i);
  float dc_vr(int i);
  float dc_xsc(int i);
  float dc_ysc(int i);
  float dc_zsc(int i);
  float dc_cxsc(int i);
  float dc_cysc(int i);
  float dc_czsc(int i);
  float dc_c2(int i);
  int ec_stat(int i);
  int ec_sect(int i);
  int ec_whol(int i);
  int ec_inst(int i);
  int ec_oust(int i);
  float etot(int i);
  float ec_ei(int i);
  float ec_eo(int i);
  float ec_t(int i);
  float ec_r(int i);
  float ech_x(int i);
  float ech_y(int i);
  float ech_z(int i);
  float ec_m2(int i);
  float ec_m3(int i);
  float ec_m4(int i);
  float ec_c2(int i);
  int sc_sect(int i);
  int sc_hit(int i);
  int sc_pd(int i);
  int sc_stat(int i);
  float edep(int i);
  float sc_t(int i);
  float sc_r(int i);
  float sc_c2(int i);
  int cc_sect(int i);
  int cc_hit(int i);
  int cc_segm(int i);
  int nphe(int i);
  float cc_t(int i);
  float cc_r(int i);
  float cc_c2(int i);
  //////////
  int pidpart(int i);
  float xpart(int i);
  float ypart(int i);
  float zpart(int i);
  float epart(int i);
  float pxpart(int i);
  float pypart(int i);
  float pzpart(int i);
  float qpart(int i);
  float weight();

  float mcp(int i);
  float mctheta(int i);
  float mcphi(int i);
  int mcid(int i);
};

#endif
