#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>


void Tree_Reader(){



  gROOT->ProcessLine(".L ./Loader.C+");
  TFile *f = new TFile("skim4_Tree1.root");
  TTree *t1 = (TTree*)f->Get("skim4_Tree");

  vector<TLorentzVector> *v_p4=0;


  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;


  vector<TLorentzVector> *v_vertex=0;

  vector<double> *v_beta=0;

  Double_t start_time;
  vector<double> *energy=0;
  vector<double> *charge=0;
  vector<double> *PID=0;
  vector<double> *chi2PID=0;

  Int_t readchargetracks;
  Int_t readprotonno;
  Int_t readpipno;
  Int_t readpimno;
  Int_t readelno;
  Int_t readkaonpno;



  t1->SetBranchAddress("p4",&v_p4);

  t1->SetBranchAddress("vertex",&v_vertex);

  t1->SetBranchAddress("beta",&v_beta);

  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);

  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&energy);
  t1->SetBranchAddress("charge",&charge);
  t1->SetBranchAddress("PID",&PID);
  t1->SetBranchAddress("chi2PID",&chi2PID);
  t1->SetBranchAddress("chargetracks",&readchargetracks);
  t1->SetBranchAddress("protonno",&readprotonno);
  t1->SetBranchAddress("pipno",&readpipno);
  t1->SetBranchAddress("pimno",&readpimno);
  t1->SetBranchAddress("elno",&readelno);
  t1->SetBranchAddress("kaonpno",&readkaonpno);

  TFile fileOutput1("Reduction_Factor.root","recreate");

  //Creating histograms for delta beta of K^+ after different cuts
  auto* h1=new TH1F("1","#Delta#beta of K^{+};P [GeV];Counts",200,0,13);
  auto* h2=new TH2F("2","#Delta#beta of K^{+};P [GeV];Counts",200,0,13,200,-0.2,0.2);
  auto* h3=new TH2F("3","#Delta#beta of K^{+};P [GeV];Counts",200,0,13,200,-0.2,0.2);
  auto* h4=new TH1F("4","#Delta#beta of K^{+};P [GeV];Counts",200,0,13);
  auto* h5=new TH1F("5","#Delta#beta of K^{+};P [GeV];Counts",200,0,13);
  auto* h6=new TH1F("6","#Delta#beta of K^{+};P [GeV];Counts",200,0,13);
  auto* h7=new TH1F("7","#Delta#beta of K^{+};P [GeV];Counts",200,0,13);
  auto* h8=new TH1F("8","#Delta#beta of K^{+};P [GeV];Counts",200,0,13);

  auto* hmiss_1=new TH1F("miss_1","MM(e' K^{+});P [GeV];Counts",200,-2,5);
  auto* hmiss_2=new TH1F("miss_2","MM(e' K^{+} K^{+});P [GeV];Counts",200,-2,5);

  auto* hl0_1=new TH1F("l0_1","M(p #pi^{-});P [GeV];Counts",200,0,5);

  auto* hP=new TH1F("P","momentum of K^{+};P [GeV];Counts",200,0,13);
  auto* hP1=new TH1F("P1","momentum of K^{+};P [GeV];Counts",200,0,13);
  auto* hP2=new TH1F("P2","momentum of K^{+};P [GeV];Counts",200,0,13);
  auto* hP3=new TH1F("P3","momentum of K^{+};P [GeV];Counts",200,0,13);

  vector<TLorentzVector> v_kp;

  TLorentzVector miss;
  TLorentzVector miss2;

  TLorentzVector l0;

  Double_t beta_tof_pip;
  Double_t P_pip;
  Double_t beta_calc_pip;
  Double_t delta_beta_pip;

  Double_t beta_tof_pim;
  Double_t P_pim;
  Double_t beta_calc_pim;
  Double_t delta_beta_pim;

  Double_t beta_tof_pr;
  Double_t P_pr;
  Double_t beta_calc_pr;
  Double_t delta_beta_pr;

  Double_t P_el;

  vector<Double_t> v_beta_tof_kp;
  Double_t beta_tof_kp;
  vector<Double_t> v_P_kp;
  Double_t P_kp;
  vector<Double_t> v_beta_calc_kp;
  Double_t beta_calc_kp;
  vector<Double_t> v_delta_beta_kp;
  Double_t delta_beta_kp;

  TLorentzVector el;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;




  Long64_t nentries = t1->GetEntries();

  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);
    if (i % 1000 == 0){
      fprintf (stderr, "%lld\r", i/1000);
      fflush (stderr);
        }
      Int_t Nparticles = v_p4->size();

      for(Int_t j=0; j<Nparticles; j++){

        if(PID->at(j)==11){
          el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        }
        else if(PID->at(j)==211){
          pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pip = v_beta->at(j);
        }
        else if(PID->at(j)==-211){
          pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pim = v_beta->at(j);
        }
        else if(PID->at(j)==2212){
          pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pr = v_beta->at(j);
        }
        else if(PID->at(j)==321){
          v_kp.clear();
          v_beta_tof_kp.clear();
          kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_kp = v_beta->at(j);
          v_kp.push_back(kp);
          v_beta_tof_kp.push_back(beta_tof_kp);

        }

      }

      P_pip = sqrt((pow(pip.Px(),2))+(pow(pip.Py(),2))+(pow(pip.Pz(),2)));
      beta_calc_pip = P_pip/(sqrt((pow(P_pip,2))+(pow(pip.M(),2))));
      delta_beta_pip = beta_calc_pip-beta_tof_pip;

      P_pim = sqrt((pow(pim.Px(),2))+(pow(pim.Py(),2))+(pow(pim.Pz(),2)));
      beta_calc_pim = P_pim/(sqrt((pow(P_pim,2))+(pow(pim.M(),2))));
      delta_beta_pim = beta_calc_pim-beta_tof_pim;

      P_pr = sqrt((pow(pr.Px(),2))+(pow(pr.Py(),2))+(pow(pr.Pz(),2)));
      beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(pr.M(),2))));
      delta_beta_pr = beta_calc_pr-beta_tof_pr;

      P_el = sqrt((pow(el.Px(),2))+(pow(el.Py(),2))+(pow(el.Pz(),2)));

      P_kp = sqrt((pow(kp.Px(),2))+(pow(kp.Py(),2))+(pow(kp.Pz(),2)));
      beta_calc_kp = P_kp/(sqrt((pow(P_kp,2))+(pow(kp.M(),2))));
      delta_beta_kp = beta_calc_kp-beta_tof_kp;

      for(int k=0;k<readkaonpno;k++){
        v_P_kp.clear();
        v_beta_calc_kp.clear();
        v_delta_beta_kp.clear();
        v_P_kp.push_back(P_kp);
        v_beta_calc_kp.push_back(beta_calc_kp);
        v_delta_beta_kp.push_back(delta_beta_kp);
      }
      //K^+ lambda^0 channel
      miss = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - kp;
      l0 = pr + pim;
      miss2 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - v_kp[0] - v_kp[1];

      hP->Fill(P_kp);
      hP1->Fill(v_P_kp[0]);
      hP2->Fill(v_P_kp[1]);
      hP3->Fill(v_P_kp[2]);

      if(readkaonpno==1){
        hmiss_1->Fill(miss.M());
        hl0_1->Fill(l0.M());
      }

      if(readkaonpno==2){
        hmiss_2->Fill(miss2.M());

      }
      //K^+ K^+ cascade channel
      // miss2 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - kp -kp;


      if(readkaonpno>0){
        for(int k=0;k<readkaonpno;k++){
          h1->Fill(v_P_kp[k]);
          if(abs(v_delta_beta_kp[k])<0.02 && v_P_kp[k]>0.2)h4->Fill(v_P_kp[k]);
        }
      }
      if(readkaonpno>1){
        for(int k=0;k<readkaonpno;k++){
          h2->Fill(v_P_kp[k],v_delta_beta_kp[k]);
          h5->Fill(v_P_kp[k]);
        }
      }
      if(readkaonpno>2){
        for(int k=0;k<readkaonpno;k++){
          h3->Fill(v_P_kp[k],v_delta_beta_kp[k]);

        }
        h6->Fill(v_P_kp[0]);
        h7->Fill(v_P_kp[1]);
        h8->Fill(P_kp);

      }

  }


fileOutput1.Write();


}
