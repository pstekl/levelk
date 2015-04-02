#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "functions.hh"

//A is 1
//B is 0





void cascadestats(levelkworld<double> w, size_t tn){

  //length
  std::vector<std::vector<double> > lcas(18, std::vector<double>(tn,0));
  //length weigted with probability
  std::vector<std::vector<double> > lcasp(18, std::vector<double>(tn,0));


  //freq switch over
  std::vector<double> fsw(8,0);


  //create sequences
  for(size_t j=0; j < (1<<(tn)); ++j){

    boost::dynamic_bitset<> x(tn, j);
    std::vector<size_t> ct(tn,0.);

    // for (boost::dynamic_bitset<>::size_type k = 0; k < x.size(); ++k)
    //   std::cout << x[k];
    // std::cout << "\n";

    coltab(x, ct);
    lencas(w, ct, lcasp, lcas);
    freqswitchover(w, ct, fsw);

    // for(int i=0; i<4; ++i)
    //   std::cout<<fsw[i]<<" ";
    // std::cout<<std::endl;
  }

  // std::vector<std::vector<double> > freq(9, std::vector<double>(1,0.));
  // outputlcas(tn, "freqcascades.csv", freq);

  //output lcas
  outputlcas(tn, "lencascades.csv", lcas);
  outputlcas(tn, "lencascadesp.csv", lcasp);
  outputsw("switchover.csv", fsw);

  return ;
}






//calculate and output probability of sequence,
// efficiency and public belief
template <typename T>
void seqprob(levelkworld<T> w, size_t tn){


  //probabilities for the 9 cases + 1 chksum + 1 efficiency +1 public belief
  std::vector<std::vector<double> > pt(12, std::vector<double>(tn,0.));

  //create sequences
  for(size_t i=0; i < tn ; ++i) {
    for(size_t j=0; j < (1<<(i+1)); ++j){

      boost::dynamic_bitset<> x(i+1, j);
      std::vector<size_t> ct(i+1,0.);

      // for (boost::dynamic_bitset<>::size_type k = 0; k < x.size(); ++k)
      //   std::cout << x[k];
      // std::cout << "\n";

      coltab(x, ct);

      // for (size_t k = 0; k < x.size(); ++k)
      //   std::cout << w.p[ct[k]]<<", ";
      // std::cout << "\n";



      std::vector<double> tab2res(9,0.);
      tab2(w,ct,i,tab2res);

      for(size_t k=0; k<9; ++k){
        pt[k][i] += tab2res[k];
        pt[9][i] += pt[k][i];
      }

    }
    //chk sum
    pt[9][i] = 0.;
    for(size_t k=0; k<9; ++k)
      pt[9][i] += pt[k][i];
  }


  //calculate efficiency
  efficiency(w, pt);
  //calculate public belief
  publicbelief(w, pt);
  //output results
  output(w, tn, pt);

  return;
}





int main() {

  //read ini file
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini("../config.ini", pt);
  size_t t0 = pt.get<size_t>("levelkworld.t0");
  size_t tn = pt.get<size_t>("levelkworld.tn");

  //create lvl k world
  levelkworld<double> w(pt.get<double>("levelkworld.lvl0ratio"), pt.get<double>("levelkworld.lvl1ratio"),
                        pt.get<double>("levelkworld.lvl2ratio"), pt.get<double>("levelkworld.lvl3ratio"),
                        pt.get<double>("levelkworld.signal") );


  seqprob(w, tn);
  cascadestats(w, tn);




  return EXIT_SUCCESS;
}
