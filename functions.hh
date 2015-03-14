
typedef std::numeric_limits< double > dbl;

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
  T p[9];
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
    p[0] = nc();
    p[1] = l2l3A();
    p[2] = l2l3B();
    p[3] = l2A();
    p[4] = l2B();
    p[5] = l3A();
    p[6] = l3B();
    p[7] = l2Al3B();
    p[8] = l2Bl3A();

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
int level2(boost::dynamic_bitset<> seq, int signal, int t, int& cascade){
  int A=0;
  int B=0;
  //just observe until t-1
  for (boost::dynamic_bitset<>::size_type i = 0; i < t-1; ++i) {
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
int level3(boost::dynamic_bitset<> seq, int signal, int t, int& cascade){
  int A=0;
  int B=0;
  int level3A=0;
  int level3B=0;


  //just observe until t-1
  for (boost::dynamic_bitset<>::size_type i = 0; i < t-1; ++i) {

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
void tab2(levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, int i, std::vector<double>& res){

  double tmp=1.;
  //check for all 9 cases
  for(int k=0; k<9; ++k) {
    if( (ct[i] == w.p[k]) || (ct[i] == 1.-w.p[k]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res[k] = tmp;
      //speed up computation, because these others results must be zero
      for(int l=k+1; l<9; ++l){
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
void freqswitchover(levelkworld<T> w, std::vector<double>& ct,  std::vector<int>& res){

  //cascades can only emerge if size > 2
  if(ct.size() < 3)
    return;

  int tmp=  -1;// last value in ct 0 A 1 B , -1 anything else
  int lastcas= -1;  //-1 unset or no cascade, 0 A cascade 1 B cascade
  //the first two values are always no cascade
  for(int i=2; i < ct.size(); ++i){

    //A cascade
    if(  ((ct[i] == w.p[1]) || (ct[i] == 1. - w.p[1] )) ) {

      // tmp contains informaiton about cascade of bit before
      if ( tmp == 0 )  {//already in A cascade=> do nothing
        continue;
      }
      if(lastcas == 0) { //switch over A to A cascade
        res[0] += 1;
      }
      if(lastcas == 1) { //switch over B to A cascade
        res[1] += 1;
      }
      lastcas = 0;
      tmp = 0;
      continue;
    }

    //B cascade
    if(  ((ct[i] == w.p[2]) || (ct[i] == 1. - w.p[2] )) ) {

      // tmp contains informaiton about cascade of bit before
      if ( tmp == 1 )  {//already in B cascade=> do nothing
        continue;
      }
      if(lastcas == 0) { //switch over A to B cascade
        res[2] += 1;
      }
      if(lastcas == 1) { //switch over B to B cascade
        res[3] += 1;
      }
      lastcas = 1;
      tmp = 1;
      continue;
    }
    //anything else in ct[i] but A or B cascade
    tmp = -1;

  }

}




//output
template <typename T>
void outputsw(levelkworld<T>& w, int tn, std::string filename, std::vector<int>& pt){


  std::ofstream output;
  output.precision(dbl::digits10);
  output.setf( std::ios::fixed);
  output.open (filename);



  std::vector<std::string> rows{"A to A", "A to B", "B to A", "B to B"};


  output<<std::endl;

  output <<"Switch over:";
  output << "\n";
  for(int j=0; j<pt.size(); ++j){
    output<<rows[j];
    output<<","<<pt[j];
    output<<"\n";
  }

  output.close();

}








//frequency of cascades
template <typename T>
void freqcas(levelkworld<T> w, std::vector<double>& ct,  std::vector<std::vector<double>>& res){

  //cascades can only emerge if size > 2
  if(ct.size() < 3)
    return;


  //the first two values are always no cascade
  for(int i=2; i < ct.size(); ++i){

    //check for all 9 cases
    for(int k=0; k<9; ++k) {
      int len=0;
      for( ; ( len<ct.size() && ( (ct[i+len] == w.p[k]) || (ct[i+len] == 1. - w.p[k] ))  )   ; ++len ){
      }
      if(len > 0){

        double tmp=1.;
        // for(int j=0; j<i+len-1; ++j){
        //     tmp *= ct[j];
        // }
        for(int j=0; j<ct.size(); ++j){
            tmp *= ct[j];
        }
        res[k][len-1] += tmp;  //weighted frequency
        //res[k][len-1] += 1; //absolute frequency
        i += len-1;
        //std::cout<<"k "<<k<<" len  "<<len<<"res[k][len-1]   "<<res[k][len-1]<<std::endl;
        break; //only 1 case is possible
      }
    }
  }



}







//length of cascades
template <typename T>
void lencas(levelkworld<T> w, std::vector<double>& ct,  std::vector<std::vector<double>>& res){

  //cascades can only emerge if size > 2
  if(ct.size() < 3)
    return;


  //the first two values are always no cascade
  for(int i=2; i < ct.size(); ++i){

    //check for all 9 cases
    for(int k=0; k<9; ++k) {
      int len=0;
      for( ; ( len<ct.size() && ( (ct[i+len] == w.p[k]) || (ct[i+len] == 1. - w.p[k] ))  )   ; ++len ){
      }
      if(len > 0){

        double tmp=1.;
        // for(int j=0; j<i+len-1; ++j){
        //     tmp *= ct[j];
        // }
        for(int j=0; j<ct.size(); ++j){
            tmp *= ct[j];
        }
        res[k][len-1] += tmp*len;  //weighted len
        //res[k][len-1] += len; //absolute length
        i += len-1;
        //std::cout<<"k "<<k<<" len  "<<len<<"res[k][len-1]   "<<res[k][len-1]<<std::endl;
        break; //only 1 case is possible
      }
    }
  }



}


//output
template <typename T>
void outputlcas(levelkworld<T>& w, int tn, std::string filename, std::vector<std::vector<double>>& pt){


  std::ofstream output;
  output.precision(dbl::digits10);
  output.setf( std::ios::fixed);
  output.open (filename);

  // //level k world parameters
  // std::vector<std::string> lkwrows{"l0ratio", "l1ratio", "l2ratio", "l3ratio"};
  // output << "level k world parameters\n";
  // for(int i=0; i<4; ++i){
  //   output << lkwrows[i]<<","<<w.lkr[i] <<"\n";
  // }
  // output <<"signal,"<<w.q<<"\n\n";

  std::vector<std::string> rows{"NC", "l2Al3A", "l2Bl3B", "l2Al3NC", "l2Bl3NC", "l3Al2NC", "l3Bl2NC", "l2Al3B", "l2Bl3A"};

  // output<<"P(at|w=A)";
  // for(int i=0; i<9; ++i){
  //   output<<","<<rows[i];
  // }
  // output<<std::endl;
  // output<<"A";
  // for(int i=0; i<9; ++i){
  //   output<<","<<w.p[i];
  // }
  // output<<std::endl;
  // output<<"B";
  // for(int i=0; i<9; ++i){
  //   output<<","<<1.- w.p[i];
  // }
  // output<<std::endl;
  output<<std::endl;

  output <<"length";
  for(int i=0; i<tn;++i) {
    output<<","<<i+1;
  }
  output << "\n";
  for(int j=0; j<pt.size(); ++j){
    output<<rows[j];
    for(int i=0; i < pt[j].size(); ++i){
      output<<","<<pt[j][i];
    }
    output<<"\n";
  }

  output.close();

}











//compute tab1 coloured table
template <typename T>
T coltab(boost::dynamic_bitset<>& x, levelkworld<T> w, std::vector<double>& res){



  for(int i=0; i<x.size(); ++i){
    int l2c=-1;
    int l3c=-1;

    //check in which state the levelks would be at bit i
    //@param cascade =0 NC, 1 AC, 2 BC
    level2(x, 0, i+1, l2c);
    level3(x, 0, i+1, l3c);
    //std::cout<<"l2c: "<<l2c<<std::endl;
    //std::cout<<"l3c: "<<l3c<<std::endl;
    //both A cascade
    if(l2c==1 && l3c==1) {
      //res *= w.p[1];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[1];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[1];
      }
      continue;
    }
    //both B cascade
    if(l2c==2 && l3c==2) {
      //res *= w.p[2];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[2];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[2];
      }
      continue;
    }
    // l2 A, l3 nc
    if(l2c==1 && l3c==0) {
      //res *= w.p[3];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[3];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[3];
      }
      continue;
    }
    //l2B, l3 nc
    if(l2c==2 && l3c==0) {
      //res *= w.p[4];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[4];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[4];
      }
      continue;
    }
    //l3A, l2 nc
    if(l2c==0 && l3c==1) {
      //res *= w.p[5];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[5];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[5];
      }
      continue;
    }
    //l3B, l2nc
    if(l2c==0 && l3c==2) {
      //res *= w.p[6];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[6];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[6];
      }
      continue;
    }
    //no cascade
    if(l2c==0 && l3c==0) {
      //res *= w.p[0];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[0];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[0];
      }
      continue;
    }
    // l2A l3 B
    if(l2c==1 && l3c==2) {
      //res *= w.p[0];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[7];
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[7];
      }
      continue;
    }
    // l2B l3 A
    if(l2c==2 && l3c==1) {
      //res *= w.p[0];
      // A push back prob for A
      if(x[i] == 1) {
        res[i] = w.p[8] ;
      }
      else{ // B push back prob for B
        res[i] = 1.- w.p[8];
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
void efficiency(levelkworld<T>& w, std::vector<std::vector<double>>& pt) {
  for(int i=0; i<pt[0].size(); ++i){
    pt[10][i] = 0.5*w.lkr[0] + w.q * w.lkr[1] + w.lkr[2] * (  w.q*( pt[0][i] + pt[5][i] + pt[6][i] ) + pt[4][i] + pt[1][i] + pt[7][i]  )
      + w.lkr[3] * ( w.q*( pt[0][i] + pt[3][i] + pt[4][i] ) + pt[5][i] + pt[1][i] + pt[8][i]     );
  }
}

//calculate public belief
template <typename T>
void publicbelief(levelkworld<T>& w, std::vector<std::vector<double>>& pt) {
  double pbnc = 0.5*w.lkr[0] + w.q * (w.lkr[1] + w.lkr[2] + w.lkr[3]);
  double pbcr = 0.5*w.lkr[0] + w.q * w.lkr[1] + w.lkr[2] + w.lkr[3];
  double pbl2cl3nc = 0.5*w.lkr[0] + w.q * (w.lkr[1] + w.lkr[3]) + w.lkr[2];
  double pbl2ncl3c = 0.5*w.lkr[0] + w.q * (w.lkr[1] + w.lkr[2]) + w.lkr[3];
  for(int i=0; i<pt[0].size(); ++i){
    pt[11][i] = pbnc*pt[0][i] + pbcr*( pt[1][i] + pt[2][i] + pt[7][i] + pt[8][i] )
      + pbl2cl3nc * ( pt[3][i] + pt[4][i]  )
      + pbl2ncl3c * ( pt[5][i] + pt[6][i]   );
  }
}



//output
template <typename T>
void output(levelkworld<T>& w, int tn, std::vector<std::vector<double>>& pt){


  std::ofstream output;
  output.precision(dbl::digits10);
  output.setf( std::ios::fixed);
  output.open ("result.csv");

  //level k world parameters
  std::vector<std::string> lkwrows{"l0ratio", "l1ratio", "l2ratio", "l3ratio"};
  output << "level k world parameters\n";
  for(int i=0; i<4; ++i){
    output << lkwrows[i]<<","<<w.lkr[i] <<"\n";
  }
  output <<"signal,"<<w.q<<"\n\n";

  std::vector<std::string> rows{"NC", "l2Al3A", "l2Bl3B", "l2Al3NC", "l2Bl3NC", "l3Al2NC", "l3Bl2NC", "l2Al3B", "l2Bl3A",
      "CHKSUM", "efficiency", "publicbelief"};

  output<<"P(at|w=A)";
  for(int i=0; i<9; ++i){
    output<<","<<rows[i];
  }
  output<<std::endl;
  output<<"A";
  for(int i=0; i<9; ++i){
    output<<","<<w.p[i];
  }
  output<<std::endl;
  output<<"B";
  for(int i=0; i<9; ++i){
    output<<","<<1.- w.p[i];
  }
  output<<std::endl;
  output<<std::endl;

  output <<"t";
  for(int i=0; i<tn;++i) {
    output<<","<<i+1;
  }
  output << "\n";
  for(int j=0; j<pt.size(); ++j){
    output<<rows[j];
    for(int i=0; i < pt[j].size(); ++i){
      output<<","<<pt[j][i];
    }
    output<<"\n";
  }

  output.close();

}








// //compute probability of sequence
// template <typename T>
// T seqprob(boost::dynamic_bitset<>& x, T q, levelkworld<T> lkw){
//   T res=1.;

//   for(int i=0; i<x.size(); ++i){
//   int l2c=-1;
//   int l3c=-1;
//     //create copy of bitset up to bit i
//     boost::dynamic_bitset<> y(i);
//     for(int j=0; j<i; ++j){
//       y[j] = x[j];
//     }
//     //check in which state the levelks would be at bit i
//     //@param cascade =0 NC, 1 AC, 2 BC
//     level2(y, 0, l2c);
//     level3(y, 0, l3c);

//     //std::cout<<"l2c: "<<l2c<<std::endl;
//     //std::cout<<"l3c: "<<l3c<<std::endl;

//     //both A cascade
//     if(l2c==1 && l3c==1) {
//       res *= lkw.p[1];
//       continue;
//     }
//     //both B cascade
//     if(l2c==2 && l3c==2) {
//       res *= lkw.p[2];
//       continue;
//     }
//     // l2 A, l3 nc
//     if(l2c==1 && l3c==0) {
//        res *= lkw.p[3];
//        continue;
//     }
//     //l2B, l3 nc
//     if(l2c==2 && l3c==0) {
//       res *= lkw.p[4];
//       continue;
//     }
//     //l3A, l2 nc
//     if(l2c==0 && l3c==1) {
//       res *= lkw.p[5];
//       continue;
//     }
//     //l3B, l2nc
//     if(l2c==0 && l3c==2) {
//        res *= lkw.p[6];
//        continue;
//     }
//     //no cascade
//     if(l2c==0 && l3c==0) {
//       res *= lkw.p[0];
//       continue;
//     }
//     else{
//       std::cout<<"ERROR: unspecified case in seq prob"<<std::endl;
//     }

//   }

//   return res;
// }
