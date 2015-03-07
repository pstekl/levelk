#include <iostream>
#include <random>
#include <string>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include "functions.hh"

//A is 1
//B is 0



int signal(double q = 0.6) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> distrib({ 1-q, q });
  return distrib(gen);
}


//new test

// //new testing
// template <typename T>
// T prob(boost::dynamic_bitset<>& x, levelkworld<T> lkw){
//   T tmp=1.;

//   std::vector<int> cases;
//   std::vector<T> res;

//   for(int i=0; i<x.size(); ++i){
//     int l2c=-1;
//     int l3c=-1;
//     //create copy of bitset up to bit i
//     boost::dynamic_bitset<> y(i+1);
//     for(int j=0; j<i+1; ++j){
//       y[j] = x[j];
//     }
//     //check in which state the levelks would be at bit i
//     //@param cascade =0 NC, 1 AC, 2 BC
//     level2(y, 0, l2c);
//     level3(y, 0, l3c);

//     //std::cout<<"l2c: "<<l2c<<std::endl;
//     //std::cout<<"l3c: "<<l3c<<std::endl;

//     //both A cascade : case 0
//     if(l2c==1 && l3c==1) {
//       // A push back prob for A
//       if(y[i] == 1) {
//         tmp*= lkw.p[1];
//         res.push_back(lkw.p[1]);
//       }
//       else{ // B push back prob for B
//         tmp*= (1.-lkw.p[1]);
//         res.push_back(1.- lkw.p[1]);
//       }
//       continue;
//     }
//     //both B cascade : case 1
//     if(l2c==2 && l3c==2) {
//       //res *= lkw.p[2];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[2]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[2]);
//       }
//       continue;
//     }
//     // l2 A, l3 nc : case 2
//     if(l2c==1 && l3c==0) {
//       //res *= lkw.p[3];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[3]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[3]);
//       }
//       continue;
//     }
//     //l2B, l3 nc : case 3
//     if(l2c==2 && l3c==0) {
//       //res *= lkw.p[4];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[4]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[4]);
//       }
//       continue;
//     }
//     //l3A, l2 nc : case 4
//     if(l2c==0 && l3c==1) {
//       //res *= lkw.p[5];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[5]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[5]);
//       }
//       continue;
//     }
//     //l3B, l2nc : case 5
//     if(l2c==0 && l3c==2) {
//       //res *= lkw.p[6];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[6]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[6]);
//       }
//       continue;
//     }
//     //no cascaed : case 6
//     if(l2c==0 && l3c==0) {
//       //res *= lkw.p[0];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[0]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[0]);
//       }
//       continue;
//     }

//     // l2A l3 B : case 7
//     if(l2c==1 && l3c==2) {
//       //res *= lkw.p[0];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[7]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[7]);
//       }
//       continue;
//     }
//     // l2B l3 A : case 8
//     if(l2c==2 && l3c==1) {
//       //res *= lkw.p[0];
//       // A push back prob for A
//       if(y[i] == 1) {
//         res.push_back(lkw.p[8]);
//       }
//       else{ // B push back prob for B
//         res.push_back(1.- lkw.p[8]);
//       }
//       continue;
//     }
//     else{
//       std::cout<<"ERROR: unspecified case in col tab"<<std::endl;
//     }
//   }








// template <typename T>
// void tab2nc(T q, levelkworld<T> w, boost::dynamic_bitset<>& x, std::vector<double>& ct, std::vector<double>& res){

//     // p[0] = nc();
//     // p[1] = l2l3A();
//     // p[2] = l2l3B();
//     // p[3] = l2A();
//     // p[4] = l2B();
//     // p[5] = l3A();
//     // p[6] = l3B();

//   double tmp=1.;

//   //check for no cascade
//   for(int i=0; i < ct.size(); ++i){
//     tmp=1.;
//     if( (ct[i] == w.p[0]) || (ct[i] == 1.-w.p[0]) ) {
//       for(int j=0; j<i+1; ++j){
//         tmp *= ct[j];
//       }
//       res.push_back(tmp);
//     }
//     else{
//       res.push_back(0.);
//     }

//   }

//   // std::cout<<"Table2 No cascade:"<<std::endl;
//   // for(std::vector<double>::iterator it = res.begin(); it != res.end(); ++it){
//   //   std::cout<<*it<<std::endl;
//   // }

//   return;
// }

















//   // for(std::vector<double>::iterator it = res.begin(); it != res.end(); ++it){
//   //   std::cout<<*it<<std::endl;
//   // }
//   return 0.;
// }


















template <typename T>
void table1(T q, levelkworld<T> lkw){

  int n=20;


  //probabilities for the 9 cases
  std::vector<std::vector<double>> pt(9, std::vector<double>(n));


  //create sequences
  for(unsigned long i=0; i < n ; ++i) {
    pt[0][i] = 0.;
    for(unsigned long j=0; j < (1<<(i+1)); ++j){
      std::vector<double> tab2ncres;
      std::vector<double> tab2l2l3Ares;
      std::vector<double> tab2l2l3Bres;
      std::vector<double> tab2l2Ares;
      std::vector<double> tab2l2Bres;
      std::vector<double> tab2l3Ares;
      std::vector<double> tab2l3Bres;
      std::vector<double> tab2l2Al3Bres;
      std::vector<double> tab2l2Bl3Ares;
      boost::dynamic_bitset<> x(i+1, j);
      std::vector<double> ct;
      // for (boost::dynamic_bitset<>::size_type k = 0; k < x.size(); ++k)
      //   std::cout << x[k];
      // std::cout << "\n";
      coltab(x, q, lkw, ct);


      // p[0] = nc();
      // p[1] = l2l3A();
      // p[2] = l2l3B();
      // p[3] = l2A();
      // p[4] = l2B();
      // p[5] = l3A();
      // p[6] = l3B();

      tab2nc(q, lkw, x, ct, tab2ncres);
      tab2l2l3A(q, lkw, x, ct, tab2l2l3Ares);
      tab2l2l3B(q, lkw, x, ct, tab2l2l3Bres);
      tab2l2A(q, lkw, x, ct, tab2l2Ares);
      tab2l2B(q, lkw, x, ct, tab2l2Bres);
      tab2l3A(q, lkw, x, ct, tab2l3Ares);
      tab2l3B(q, lkw, x, ct, tab2l3Bres);
      tab2l2Al3B(q, lkw, x, ct, tab2l2Al3Bres);
      tab2l2Bl3A(q, lkw, x, ct, tab2l2Bl3Ares);
      pt[0][i] += tab2ncres[i];
      pt[1][i] += tab2l2l3Ares[i];
      pt[2][i] += tab2l2l3Bres[i];
      pt[3][i] += tab2l2Ares[i];
      pt[4][i] += tab2l2Bres[i];
      pt[5][i] += tab2l3Ares[i];
      pt[6][i] += tab2l3Bres[i];
      pt[7][i] += tab2l2Al3Bres[i];
      pt[8][i] += tab2l2Bl3Ares[i];


    }
  }

  for(int j=0; j<pt.size(); ++j){
    std::cout<<"table1 "<<j<<":"<<std::endl;
    for(int i=0; i < pt[j].size(); ++i){
      std::cout<<pt[j][i]<<"\t";
    }
    std::cout<<"\n";
  }
  return;
}



// template <typename T>
// void table1(T q, levelkworld<T> lkw){

//   int n=3;


//   //probabilities for the 7 cases
//   std::vector<std::vector<double>> pt(7, std::vector<double>(n));


//  //create sequences
//   for(unsigned long i=0; i < n; ++i) {
//     pt[0][i] = 0.;
//     for(unsigned long j=0; j < (1<<(i+1)); ++j){
//       std::vector<double> tab2ncres;
//       boost::dynamic_bitset<> x(i+1, j);
//       std::vector<double> ct;
//       // for (boost::dynamic_bitset<>::size_type k = 0; k < x.size(); ++k)
//       //   std::cout << x[k];
//       // std::cout << "\n";
//       coltab(x, q, lkw, ct);
//       tab2nc(q, lkw, x, ct, tab2ncres);
//       //sum up for table1
//       for(int l=0; l< tab2ncres.size(); ++l){
//         pt[0][l] += tab2ncres[l];
//       }
//     }
//   }

//   // std::cout<<"table1 nc:"<<std::endl;
//   // for(std::vector<double>::iterator it = pt[0].begin(); it != pt[0].end(); ++it){
//   //   std::cout<<*it<<std::endl;
//   // }

// }


int main() {


  // //create input sequence
  // boost::dynamic_bitset<> seq(5); // all 0's by default
  // seq[0] = 1;
  // seq[1] = 1;
  // seq[4] = 1;
  // std::cout<<"sequence: ";
  // for (boost::dynamic_bitset<>::size_type i = 0; i < seq.size(); ++i)
  //   std::cout << seq[i];
  // std::cout << "\n";

  // //create signal
  // int q = signal(0.7);

  // std::cout<<"SIGNAL Q: "<<q<<std::endl;

  // //return decisions of levelk players with same signal and same sequence
  // std::cout<<"LEVEL0: "<<level0(seq,q)<<std::endl;
  // std::cout<<"LEVEL1: "<<level1(seq,q)<<std::endl;
  // std::cout<<"LEVEL2: "<<level2(seq,q)<<std::endl;
  // std::cout<<"LEVEL3: "<<level3(seq,q)<<std::endl;

  levelkworld<double> w(0.1,0.3,0.4,0.2,0.66666);
  //std::cout<<"p[0]: "<<w.p[0]<<std::endl;


  // // line 41 of ulf table
  // boost::dynamic_bitset<> seq(10);
  // seq[0] = 1;
  // seq[1] = 1;
  // seq[2] = 1;
  // seq[3] = 1;
  // seq[4] = 1;
  // seq[5] = 0;
  // seq[6] = 1;
  // seq[7] = 1;
  // seq[8] = 1;
  // seq[9] = 1;

  // double q = 0.6;

  // std::cout<<"probability seq AA: "<<seqprob(seq, 0.6, w)<<std::endl;
  //genseq(6, 3.);

  table1(0.6, w);
    // std::vector<double> ct;
    // coltab(seq, q, w, ct);

  return EXIT_SUCCESS;
}
