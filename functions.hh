
//#define DEBUG



typedef std::numeric_limits< double > dbl;


enum cases { Anc, Al2Al3A, Al2Bl3B, Al2Al3nc, Al2Bl3nc, Al2ncl3A, Al2ncl3B, Al2Al3B, Al2Bl3A,
             Bnc, Bl2Al3A, Bl2Bl3B, Bl2Al3nc, Bl2Bl3nc, Bl2ncl3A, Bl2ncl3B, Bl2Al3B, Bl2Bl3A};


//A is 1
//B is 0
template<typename T>
class levelkworld {
public:
  levelkworld(T l0ratio, T l1ratio, T l2ratio, T l3ratio, T q_) :
    lkr{l0ratio, l1ratio, l2ratio, l3ratio}, q(q_) {

      updatep();
    }

  T lkr[4]; //lkratio
  T q;        //prob q that signal is correct
  //calculate probability p that A comes for each situation
  //in the same situation probability for B is of course 1-p
  T p[18];
  //No cascade
  T nc() {
    return 0.5*lkr[0] + q*(lkr[1]+lkr[2]+lkr[3]);
  }
  //l2+l3 A cascade
  T l2l3A() {
    return 0.5*lkr[0] + q*lkr[1] + lkr[2]+lkr[3];
  }
  //l2+l3 B cascade
  T l2l3B() {
    return 0.5*lkr[0] + q*lkr[1];
  }
  //only l2 A cascade
  T l2A() {
     return 0.5*lkr[0] + q*(lkr[1]+lkr[3]) + lkr[2];
  }
  //only l2 B cascade
  T l2B() {
     return 0.5*lkr[0] + q*(lkr[1]+lkr[3]);
  }
  //only l3 A cascade
  T l3A() {
    return 0.5*lkr[0] + q*(lkr[1]+lkr[2]) + lkr[3];
  }
  //only l3 B cascade
  T l3B() {
    return 0.5*lkr[0] + q*(lkr[1]+lkr[2]);
  }

  // l2 A cascade l3 B cascade
  T l2Al3B() {
    return 0.5*lkr[0] + q*lkr[1] + lkr[2];
  }
  // l2 B l3 A
  T l2Bl3A() {
    return 0.5*lkr[0]+q*lkr[1]+lkr[3];
  }



  //update probabilities
  void updatep() {
    p[Anc] = nc();
    p[Al2Al3A] = l2l3A();
    p[Al2Bl3B] = l2l3B();
    p[Al2Al3nc] = l2A();
    p[Al2Bl3nc] = l2B();
    p[Al2ncl3A] = l3A();
    p[Al2ncl3B] = l3B();
    p[Al2Al3B] = l2Al3B();
    p[Al2Bl3A] = l2Bl3A();
    for(size_t i=0; i<9; ++i)
      p[i+9] = 1. - p[i];

  }
  //print probabilities for A
  void print(){
    std::cout<<"NC  "<<"ACR "<<"BCR "<<"ACF2 "<<"BCF2 "<<"ACF3 "<<"BCF3 "<<std::endl;
    std::cout<<nc()<<" "<<l2l3A()<<" "<<l2l3B()<<" "<<l2A()<<" "<<l2B()<<" "<<l3A()<<" "<<l3B()<<std::endl;
  }

};




// int signal(double q = 0.6) {
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   std::discrete_distribution<> distrib({ 1-q, q });
//   return distrib(gen);
// }


// //randomly return 0 or 1 with probability 0.5
// int level0(boost::dynamic_bitset<> seq, int signal){
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   std::uniform_int_distribution<> dis(0, 1);
//   return dis(gen);
// }
// //return signal
// int level1(boost::dynamic_bitset<> seq, int signal){
//   return signal;
// }




//detect cascade
//@param cascade =0 NC, 1 AC, 2 BC
size_t level2(boost::dynamic_bitset<> seq, size_t signal, size_t t, int& cascade){
  int A=0;
  int B=0;
  //just observe until t-1
  //for (boost::dynamic_bitset<>::size_type i = 0; i < t-1; ++i) {
    for (size_t i = 0; i < t; ++i) {
    if(seq[i] == 1) {
      A+=1;
    }
    else {
      B+=1;
    }
  }

  if(A-B>1) {
    cascade = 1;
    //std::cout<<"L2 A cascade"<<std::endl;
    return 1; //A cascade => return A
  }
  if(B-A>1) {
    cascade = 2;
    //std::cout<<"L2 B cascade"<<std::endl;
    return 0; //B cascade => return B
  }
  //no cascade => return signal
  cascade = 0;
  //std::cout<<"L2 NO cascade"<<std::endl;
  return signal;
}

//detect cascade and detect fake cascade
//@param cascade =0 NC, 1 AC, 2 BC
size_t level3(boost::dynamic_bitset<> seq, size_t signal, size_t t, int& cascade){
  int A=0;
  int B=0;
  int level3A=0;
  int level3B=0;


  //just observe until t-1
  //for (boost::dynamic_bitset<>::size_type i = 0; i < t-1; ++i) {
  for (size_t i = 0; i < t; ++i) {

    //check for level2 cascade
    if(std::abs(A-B)>1) {
      //check if already in A cascade => fake signal => do not increment A
      if(A-B > 1) {
        //std::cout<<"Already in l2 A cascade omit l3 A"<<std::endl;
        if(seq[i] == 0) {
          level3B+=1;
        }
      }
      //check if already in B cascade => fake signal => do not increment B
      if(B-A > 1) {
        //std::cout<<"Already in l2 B cascade omit l3 B"<<std::endl;
        if(seq[i] == 1) {
          level3A+=1;
        }
      }
    }
    // no cascade => increment
    else {
      if(seq[i] == 1) {
        level3A+=1;
      }
      else {
        level3B+=1;
      }
    }
    //increment to check for level2 cascade
    if(seq[i] == 1) {
      A+=1;
    }
    else {
      B+=1;
    }

  }

  if(level3A-level3B>1) {
    cascade = 1;
    //std::cout<<"L3 A cascade"<<std::endl;
    return 1; //A cascade => return A
  }
  if(level3B-level3A>1) {
    cascade = 2;
    //std::cout<<"L3 B cascade"<<std::endl;
    return 0; //B cascade => return B
  }

    // no l3 cascade => return signal
  cascade=0;
  //std::cout<<"L3 NO cascade"<<std::endl;
  return signal;

}




template <typename T>
void tab2(levelkworld<T> w, std::vector<size_t>& ct, size_t i, std::vector<double>& res){

  double tmp=1.;
  //check for all 9 cases, respectively 18 (if one counts B cases too)
  for(size_t k=0; k<9; ++k) {
    if( (ct[i] == k) || (ct[i] == k+9) ) {
      for(size_t j=0; j<i+1; ++j){
        tmp *= w.p[ct[j]];
      }
      res[k] = tmp;
      //speed up computation, because these others results must be zero
      for(size_t l=k+1; l<9; ++l){
        res[l] = 0.;
      }
      break;
    }
    else{
      res[k] = 0.;
    }
  }

  return;
}






//frequency of cascades
template <typename T>
void freqswitchover(levelkworld<T> w, std::vector<size_t>& ct,  std::vector<double>& res){

  //cascades can only emerge if size > 2
  if(ct.size() < 3)
    return;


  double pseq=1.;
  for(size_t j=0; j<ct.size(); ++j){
    pseq *= w.p[ct[j]];
  }


  //int tmp=  -1;// last value in ct 0 A 1 B, -1 anything else
  int lastcas= -1;  //-1 unset or no cascade, 0 A cascade 1 B cascade
  int blastcas= -1;  // cascade before lastcascade
  //the first two values are always no cascade
  for(size_t i=2; i < ct.size(); ++i){

#ifdef DEBUG
    std::cout<<"ct["<<i<<"]: "<<ct[i]<<std::endl;
#endif



    //A cascade
    if((ct[i] == Al2Al3A) || (ct[i] == Al2Al3B ) || (ct[i] == Al2Al3nc)) {
    acas:
      for(;  i<ct.size(); ++i) {
#ifdef DEBUG
        std::cout<<"A cascade"<<std::endl;
#endif
        //broken A cascade
        if((ct[i] == Bl2Al3A) || (ct[i] == Bl2Al3B ) || (ct[i] == Bl2Al3nc)) {
          for(; i<ct.size(); ++i){
#ifdef DEBUG
            std::cout<<"broken A cascade"<<std::endl;
#endif
            if( (ct[i] == Al2Al3A) || (ct[i] == Al2Al3B ) || (ct[i] == Al2Al3nc) ) {
              // switch over A to A cascade
              res[0] += 1;
              res[4] += pseq;
#ifdef DEBUG
              std::cout<<"switch over A to A cascade"<<std::endl;
#endif
              goto acas;
            }
            // switch over A to B cascade
            if( (ct[i] == Bl2Bl3A) || (ct[i] == Bl2Bl3B ) || (ct[i] == Bl2Bl3nc) ) {
              res[1] += 1;
              res[5] += pseq;
#ifdef DEBUG
              std::cout<<"switch over A to B cascade"<<std::endl;
#endif
              goto bcas;
            }
          }
        }
      }
    }
    //B cascade
    if((ct[i] == Bl2Bl3A) || (ct[i] == Bl2Bl3B ) || (ct[i] == Bl2Bl3nc)) {
    bcas:
      for(;  i<ct.size(); ++i) {
#ifdef DEBUG
        std::cout<<"B cascade"<<std::endl;
#endif
        //broken B cascade
        if((ct[i] == Al2Bl3A) || (ct[i] == Al2Bl3B ) || (ct[i] == Al2Bl3nc) ) {
          for(; i<ct.size(); ++i) {
#ifdef DEBUG
            std::cout<<"broken B cascade"<<std::endl;
#endif
            // switch over B to A cascade
            if((ct[i] == Al2Al3A) || (ct[i] == Al2Al3B ) || (ct[i] == Al2Al3nc))  {
              res[2] += 1;
              res[6] += pseq;
#ifdef DEBUG
              std::cout<<"switch over B to A cascade"<<std::endl;
#endif
              goto acas;
            }
            // switch over B to B cascade
            if((ct[i] == Bl2Bl3A) || (ct[i] == Bl2Bl3B ) || (ct[i] == Bl2Bl3nc)) {
              res[3] += 1;
              res[7] += pseq;
#ifdef DEBUG
              std::cout<<"switch over B to B cascade"<<std::endl;
#endif
              goto bcas;
            }
          }
        }
      }
    }

  }



    // // switch over A to A cascade
    // if( lastcas == 0 && blastcas == 0 && ((ct[i] == Al2Al3A) || (ct[i] == Al2Al3B ) || (ct[i] == Al2Al3nc)) ) {
    //     res[0] += 1;
    //     res[4] += pseq;
    //     blastcas = -1;
    //     continue;
    // }
    // // switch over A to B cascade
    // if(lastcas == 0 && blastcas == 0  && ((ct[i] == Bl2Bl3A) || (ct[i] == Bl2Bl3B ) || (ct[i] == Bl2Bl3nc)) ) {
    //     res[1] += 1;
    //     res[5] += pseq;
    //     lastcas = 1;
    //     blastcas = -1;
    //     continue;
    // }

    // // switch over B to A cascade
    // if( lastcas == 1 && blastcas == 1  && ((ct[i] == Al2Al3A) || (ct[i] == Al2Al3B ) || (ct[i] == Al2Al3nc)) ) {
    //     res[2] += 1;
    //     res[6] += pseq;
    //     blastcas = -1;
    //     continue;
    // }
    // // switch over B to B cascade
    // if( lastcas == 1 && blastcas == 1  && ((ct[i] == Bl2Bl3A) || (ct[i] == Bl2Bl3B ) || (ct[i] == Bl2Bl3nc)) ) {
    //     res[3] += 1;
    //     res[7] += pseq;
    //     lastcas = 1;
    //     blastcas = -1;
    //     continue;
    // }

    // // broken A cascade
    // if(lastcas == 0 && ((ct[i] == Bl2Al3A) || (ct[i] == Bl2Al3B ) || (ct[i] == Bl2Al3nc)) ) {
    //   blastcas = 0;
    //   continue;
    // }
    // // broken B cascade
    // if (lastcas == 1 && (ct[i] == Al2Bl3A) || (ct[i] == Al2Bl3B ) || (ct[i] == Al2Bl3nc) ) {
    //   blastcas = 1;
    //   continue;
    // }



}




//output
template <typename T>
void outputsw(std::string filename, std::vector<T>& pt){


  std::ofstream output;
  output.precision(dbl::digits10);
  output.setf( std::ios::fixed);
  output.open (filename);



  std::vector<std::string> rows{"A to A", "A to B", "B to A", "B to B"};


  output<<std::endl;

  output <<"Switch over:";
  output << "\n";
  for(size_t j=0; j<4; ++j){
    output<<rows[j];
    output<<","<<pt[j]<<","<<pt[j+4];
    output<<"\n";
  }

  output.close();

}





//length of cascades
template <typename T>
void lencas(levelkworld<T> w, std::vector<size_t>& ct,  std::vector<std::vector<double> >& res,  std::vector<std::vector<double> >& res2){

  //cascades can only emerge if size > 2
  if(ct.size() < 3)
    return;


  double pseq=1.;

  for(size_t j=0; j<ct.size(); ++j){
    pseq *= w.p[ct[j]];
  }




  std::vector<int> tmp{Al2Al3A, Al2Al3nc, Al2Al3B, Bl2Bl3B,  Bl2Bl3nc, Bl2Bl3A};



  //the first two values are always no cascade
  for(int i=2; static_cast<size_t>(i) < ct.size(); ++i){

    //check for all relevant cases
    for(size_t k=0; k<tmp.size(); ++k) {
      int len=0;
      // for( ; ( static_cast<size_t>(i+len)<ct.size() && ( (ct[static_cast<size_t>(i+len)] == k) || (ct[static_cast<size_t>(i+len)] == k+9 ))  )   ; ++len ){
      // }
      for( ; ( static_cast<size_t>(i+len)<ct.size() && ( (ct[static_cast<size_t>(i+len)] == tmp[k]) )  )   ; ++len ){
      }
      if(len > 0){


        res[tmp[k]][len-1] += pseq;
        res2[tmp[k]][len-1] += 1;
        i += len;
        --i;
        //std::cout<<"res[k][len-1]="<<"res["<<k<<"]"<<"["<<len-1<<"]: "<<res[tmp[k]][len-1]<<std::endl;
        break; //only 1 case is possible
      }
    }
  }



}




//output
void outputlcas(size_t tn, std::string filename, std::vector<std::vector<double> >& pt){


  std::ofstream output;
  output.precision(dbl::digits10);
  output.setf( std::ios::fixed);
  output.open (filename);

  // //level k world parameters
  // std::vector<std::string> lkwrows{"l0ratio", "l1ratio", "l2ratio", "l3ratio"};
  // output << "level k world parameters\n";
  // for(size_t i=0; i<4; ++i){
  //   output << lkwrows[i]<<","<<w.lkr[i] <<"\n";
  // }
  // output <<"signal,"<<w.q<<"\n\n";

  //std::vector<std::string> rows{"NC", "l2Al3A", "l2Bl3B", "l2Al3NC", "l2Bl3NC", "l3Al2NC", "l3Bl2NC", "l2Al3B", "l2Bl3A"};
  std::vector<std::string> rows{"Al2Al3A", "Al2Al3nc", "Al2Al3B", "Bl2Bl3B",  "Bl2Bl3nc", "Bl2Bl3A"};

  std::vector<int> tmp{Al2Al3A, Al2Al3nc, Al2Al3B, Bl2Bl3B,  Bl2Bl3nc, Bl2Bl3A};

  // output<<"P(at|w=A)";
  // for(size_t i=0; i<9; ++i){
  //   output<<","<<rows[i];
  // }
  // output<<std::endl;
  // output<<"A";
  // for(size_t i=0; i<9; ++i){
  //   output<<","<<w.p[i];
  // }
  // output<<std::endl;
  // output<<"B";
  // for(size_t i=0; i<9; ++i){
  //   output<<","<<1.- w.p[i];
  // }
  // output<<std::endl;
  output<<std::endl;

  output <<"length";
  for(size_t i=0; i<tn;++i) {
    output<<","<<i+1;
  }
  output << "\n";
  for(size_t j=0; j<rows.size(); ++j){
    output<<rows[j];
    for(size_t i=0; i < pt[j].size(); ++i){
      output<<","<<pt[tmp[j]][i];
    }
    output<<"\n";
  }

  output.close();

}







//compute tab1 coloured table
size_t coltab(boost::dynamic_bitset<>& x, std::vector<size_t>& res){

  // p[Anc] = nc();
  //   p[Al2Al3A] = l2l3A();
  //   p[Al2Bl3B] = l2l3B();
  //   p[Al2Al3nc] = l2A();
  //   p[Al2Bl3nc] = l2B();
  //   p[Al2ncl3A] = l3A();
  //   p[Al2ncl3B] = l3B();
  //   p[Al2Al3B] = l2Al3B();
  //   p[Al2Bl3A] = l2Bl3A();

  for(size_t i=0; i<x.size(); ++i){
    int l2c=-1;
    int l3c=-1;

    //check in which state the levelks would be at bit i
    //@param cascade =0 NC, 1 AC, 2 BC
    level2(x, 0, i, l2c);
    level3(x, 0, i, l3c);
    // std::cout<<"l2c: "<<l2c<<std::endl;
    // std::cout<<"l3c: "<<l3c<<std::endl;

    //no cascade
    if(l2c==0 && l3c==0) {
      // A
      if(x[i] == 1) {
        res[i] = Anc;
      }
      else{ // B
        res[i] = Bnc;
      }
      continue;
    }

    //both A cascade
    if(l2c==1 && l3c==1) {
      // A
      if(x[i] == 1) {
        res[i] = Al2Al3A;
      }
      else{ // B
        res[i] = Bl2Al3A;
      }
      continue;
    }
    //both B cascade
    if(l2c==2 && l3c==2) {
      // A
      if(x[i] == 1) {
        res[i] = Al2Bl3B;
      }
      else{ // B
        res[i] = Bl2Bl3B;
      }
      continue;
    }
    // l2 A, l3 nc
    if(l2c==1 && l3c==0) {
      // A
      if(x[i] == 1) {
        res[i] = Al2Al3nc;
      }
      else{ // B
        res[i] = Bl2Al3nc;
      }
      continue;
    }
    //l2B, l3 nc
    if(l2c==2 && l3c==0) {
      // A
      if(x[i] == 1) {
        res[i] = Al2Bl3nc;
      }
      else{ // B
        res[i] = Bl2Bl3nc;
      }
      continue;
    }
    //l3A, l2 nc
    if(l2c==0 && l3c==1) {
      // A
      if(x[i] == 1) {
        res[i] = Al2ncl3A;
      }
      else{ // B
        res[i] = Bl2ncl3A;
      }
      continue;
    }
    //l3B, l2nc
    if(l2c==0 && l3c==2) {
      // A
      if(x[i] == 1) {
        res[i] = Al2ncl3B;
      }
      else{ // B
        res[i] = Bl2ncl3B;
      }
      continue;
    }

    // l2A l3 B
    if(l2c==1 && l3c==2) {
      // A
      if(x[i] == 1) {
        res[i] = Al2Al3B;
      }
      else{ // B
        res[i] = Bl2Al3B;
      }
      continue;
    }
    // l2B l3 A
    if(l2c==2 && l3c==1) {
      // A
      if(x[i] == 1) {
        res[i] = Al2Bl3A;
      }
      else{ // B
        res[i] = Bl2Bl3A;
      }
      continue;
    }
    else{
      std::cout<<"ERROR: unspecified case in col tab"<<std::endl;
    }
  }

  return 0.;
}




//calculate efficiency
template <typename T>
void efficiency(levelkworld<T>& w, std::vector<std::vector<double> >& pt) {
	for (size_t i = 0; i<pt[0].size(); ++i) {
		pt[10][i] = 0.5*w.lkr[0] + w.q * w.lkr[1] + w.lkr[2] * ( w.q*(pt[0][i] + pt[5][i] + pt[6][i]) + pt[3][i] + pt[1][i] + pt[7][i] )
			+ w.lkr[3] * (w.q*(pt[0][i] + pt[3][i] + pt[4][i]) + pt[5][i] + pt[1][i] + pt[8][i]);
	}
}



//calculate public belief
template <typename T>
void publicbelief(levelkworld<T>& w, std::vector<std::vector<double> >& pt) {
  double pbnc = 0.5*w.lkr[0] + w.q * (w.lkr[1] + w.lkr[2] + w.lkr[3]);
  double pbcr = 0.5*w.lkr[0] + w.q * w.lkr[1] + w.lkr[2] + w.lkr[3];
  double pbl2cl3nc = 0.5*w.lkr[0] + w.q * (w.lkr[1] + w.lkr[3]) + w.lkr[2];
  double pbl2ncl3c = 0.5*w.lkr[0] + w.q * (w.lkr[1] + w.lkr[2]) + w.lkr[3];
  for(size_t i=0; i<pt[0].size(); ++i){
    pt[11][i] = pbnc*pt[0][i] + pbcr*( pt[1][i] + pt[2][i] + pt[7][i] + pt[8][i] )
      + pbl2cl3nc * ( pt[3][i] + pt[4][i]  )
      + pbl2ncl3c * ( pt[5][i] + pt[6][i]   );
  }
}



//output
template <typename T>
void output(levelkworld<T>& w, size_t tn, std::vector<std::vector<double> >& pt){


  std::ofstream output;
  output.precision(dbl::digits10);
  output.setf( std::ios::fixed);
  output.open ("result.csv");

  //level k world parameters
  std::vector<std::string> lkwrows{"l0ratio", "l1ratio", "l2ratio", "l3ratio"};
  output << "level k world parameters\n";
  for(size_t i=0; i<4; ++i){
    output << lkwrows[i]<<","<<w.lkr[i] <<"\n";
  }
  output <<"signal,"<<w.q<<"\n\n";

  std::vector<std::string> rows{"NC", "l2Al3A", "l2Bl3B", "l2Al3NC", "l2Bl3NC", "l3Al2NC", "l3Bl2NC", "l2Al3B", "l2Bl3A",
      "CHKSUM", "efficiency", "publicbelief"};

  output<<"P(at|w=A)";
  for(size_t i=0; i<9; ++i){
    output<<","<<rows[i];
  }
  output<<std::endl;
  output<<"A";
  for(size_t i=0; i<9; ++i){
    output<<","<<w.p[i];
  }
  output<<std::endl;
  output<<"B";
  for(size_t i=9; i<18; ++i){
    output<<","<<w.p[i];
  }
  output<<std::endl;
  output<<std::endl;

  output <<"t";
  for(size_t i=0; i<tn;++i) {
    output<<","<<i+1;
  }
  output << "\n";
  for(size_t j=0; j<pt.size(); ++j){
    output<<rows[j];
    for(size_t i=0; i < pt[j].size(); ++i){
      output<<","<<pt[j][i];
    }
    output<<"\n";
  }

  output.close();

}
