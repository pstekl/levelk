//A is 1
//B is 0
template<typename T>
class levelkworld {
  T lkr[4]; //lkratio
  T q;        //prob q that signal is correct
public:
  levelkworld(T l0ratio, T l1ratio, T l2ratio, T l3ratio, T q_) :
    lkr{l0ratio, l1ratio, l2ratio, l3ratio}, q(q_) {

      updatep();
    }

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

//randomly return 0 or 1 with probability 0.5
int level0(boost::dynamic_bitset<> seq, int signal){
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, 1);
  return dis(gen);
}
//return signal
int level1(boost::dynamic_bitset<> seq, int signal){
  return signal;
}
//detect cascade
//@param cascade =0 NC, 1 AC, 2 BC
int level2(boost::dynamic_bitset<> seq, int signal, int& cascade){
  int A=0;
  int B=0;
  //just observe until t-1
  for (boost::dynamic_bitset<>::size_type i = 0; i < seq.size()-1; ++i) {
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
int level3(boost::dynamic_bitset<> seq, int signal, int& cascade){
  int A=0;
  int B=0;
  int level3A=0;
  int level3B=0;


  //just observe until t-1
  for (boost::dynamic_bitset<>::size_type i = 0; i < seq.size()-1; ++i) {

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
void tab2nc(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

    // p[0] = nc();
    // p[1] = l2l3A();
    // p[2] = l2l3B();
    // p[3] = l2A();
    // p[4] = l2B();
    // p[5] = l3A();
    // p[6] = l3B();

  double tmp=1.;

  //check for no cascade
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[0]) || (ct[i] == 1.-w.p[0]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }

  }

  // std::cout<<"Table2 No cascade:"<<std::endl;
  // for(std::vector<double>::iterator it = res.begin(); it != res.end(); ++it){
  //   std::cout<<*it<<std::endl;
  // }

  return;
}
template <typename T>
void tab2l2l3A(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

    // p[0] = nc();
    // p[1] = l2l3A();
    // p[2] = l2l3B();
    // p[3] = l2A();
    // p[4] = l2B();
    // p[5] = l3A();
    // p[6] = l3B();

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[1]) || (ct[i] == 1.-w.p[1]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }

  }

}

template <typename T>
void tab2l2l3B(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

    // p[0] = nc();
    // p[1] = l2l3A();
    // p[2] = l2l3B();
    // p[3] = l2A();
    // p[4] = l2B();
    // p[5] = l3A();
    // p[6] = l3B();

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[2]) || (ct[i] == 1.-w.p[2]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }

  }

}

template <typename T>
void tab2l2A(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

    // p[0] = nc();
    // p[1] = l2l3A();
    // p[2] = l2l3B();
    // p[3] = l2A();
    // p[4] = l2B();
    // p[5] = l3A();
    // p[6] = l3B();

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[3]) || (ct[i] == 1.-w.p[3]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }

  }

}

template <typename T>
void tab2l2B(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

    // p[0] = nc();
    // p[1] = l2l3A();
    // p[2] = l2l3B();
    // p[3] = l2A();
    // p[4] = l2B();
    // p[5] = l3A();
    // p[6] = l3B();

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[4]) || (ct[i] == 1.-w.p[4]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }

  }

}


template <typename T>
void tab2l3A(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

    // p[0] = nc();
    // p[1] = l2l3A();
    // p[2] = l2l3B();
    // p[3] = l2A();
    // p[4] = l2B();
    // p[5] = l3A();
    // p[6] = l3B();

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[5]) || (ct[i] == 1.-w.p[5]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }

  }

}


template <typename T>
void tab2l3B(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

    // p[0] = nc();
    // p[1] = l2l3A();
    // p[2] = l2l3B();
    // p[3] = l2A();
    // p[4] = l2B();
    // p[5] = l3A();
    // p[6] = l3B();

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[6]) || (ct[i] == 1.-w.p[6]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }

  }

}



template <typename T>
void tab2l2Al3B(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[7]) || (ct[i] == 1.-w.p[7]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }
  }

}


template <typename T>
void tab2l2Bl3A(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

  double tmp=1.;

  //check for l2l3A
  for(int i=0; i < ct.size(); ++i){
    tmp=1.;
    if( (ct[i] == w.p[8]) || (ct[i] == 1.-w.p[8]) ) {
      for(int j=0; j<i+1; ++j){
        tmp *= ct[j];
      }
      res.push_back(tmp);
    }
    else{
      res.push_back(0.);
    }
  }

}






//compute tab1 coloured table
template <typename T>
T coltab(boost::dynamic_bitset<>& x, T q, levelkworld<T> lkw, std::vector<T>& res){
  //T res=1.;

  for(int i=0; i<x.size(); ++i){
  int l2c=-1;
  int l3c=-1;
    //create copy of bitset up to bit i
    boost::dynamic_bitset<> y(i+1);
    for(int j=0; j<i+1; ++j){
      y[j] = x[j];
    }
    //check in which state the levelks would be at bit i
    //@param cascade =0 NC, 1 AC, 2 BC
    level2(y, 0, l2c);
    level3(y, 0, l3c);

    //std::cout<<"l2c: "<<l2c<<std::endl;
    //std::cout<<"l3c: "<<l3c<<std::endl;

    //both A cascade
    if(l2c==1 && l3c==1) {
      //res *= lkw.p[1];

      // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[1]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[1]);
      }
      continue;
    }
    //both B cascade
    if(l2c==2 && l3c==2) {
      //res *= lkw.p[2];
            // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[2]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[2]);
      }
      continue;
    }
    // l2 A, l3 nc
    if(l2c==1 && l3c==0) {
      //res *= lkw.p[3];
             // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[3]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[3]);
      }
       continue;
    }
    //l2B, l3 nc
    if(l2c==2 && l3c==0) {
      //res *= lkw.p[4];
            // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[4]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[4]);
      }
      continue;
    }
    //l3A, l2 nc
    if(l2c==0 && l3c==1) {
      //res *= lkw.p[5];
            // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[5]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[5]);
      }
      continue;
    }
    //l3B, l2nc
    if(l2c==0 && l3c==2) {
      //res *= lkw.p[6];
            // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[6]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[6]);
      }
       continue;
    }
    //no cascade
    if(l2c==0 && l3c==0) {
      //res *= lkw.p[0];
            // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[0]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[0]);
      }
      continue;
    }

   // l2A l3 B
    if(l2c==1 && l3c==2) {
      //res *= lkw.p[0];
            // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[7]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[7]);
      }
      continue;
    }
   // l2B l3 A
    if(l2c==2 && l3c==1) {
      //res *= lkw.p[0];
            // A push back prob for A
      if(y[i] == 1) {
        res.push_back(lkw.p[8]);
      }
      else{ // B push back prob for B
        res.push_back(1.- lkw.p[8]);
      }
      continue;
    }
    else{
      std::cout<<"ERROR: unspecified case in col tab"<<std::endl;
    }

  }
  // for(std::vector<double>::iterator it = res.begin(); it != res.end(); ++it){
  //   std::cout<<*it<<std::endl;
  // }
  return 0.;
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
