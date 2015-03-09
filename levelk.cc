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





template <typename T>
void table1(T q, levelkworld<T> lkw){

  int n=15;


  //probabilities for the 9 cases
  std::vector<std::vector<double>> pt(9, std::vector<double>(n));

  std::vector<double> chksum(n,0.);


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

      for(int k=0; k<9; ++k)
        chksum[i] += pt[k][i];
  }

  for(int j=0; j<pt.size(); ++j){
    std::cout<<"table1 "<<j<<":"<<std::endl;
    for(int i=0; i < pt[j].size(); ++i){
      std::cout<<pt[j][i]<<"\t";
    }
    std::cout<<"\n";
  }
  std::cout<<"CHK SUM:"<<std::endl;
  for(int i=0; i < chksum.size(); ++i){
    std::cout<<chksum[i]<<"\t";
  }
  std::cout<<"\n";
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
