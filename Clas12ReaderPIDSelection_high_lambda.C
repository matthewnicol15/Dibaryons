#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void Clas12ReaderPIDSelection_high_lambda(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();



  //Choosing file to run macro on
  TString inputFile("/work/clas12/rg-a/trains/v16_v2/skim4_inclusive/skim4_5*.hipo");

  cout<<"Analysing hipo file "<<inputFile<<endl;

  TChain fake("hipo");
  fake.Add(inputFile.Data());

  //get the hipo data
  //   reader.open(inputFile.Data());
  auto files=fake.GetListOfFiles();

  //particle identification
  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.6,10.6); //use for skim4_5038
  //TLorentzVector beam(0,0,7.6,7.6);  //use for skim4_5700
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  //TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
  TLorentzVector kp(0,0,0,db->GetParticle(321)->Mass());

  TFile fileOutput1("skim4_all_high_lambda.root","recreate");

  //creating histograms for missing mass, invariant mass and delta beta
  //MM(e' K^{+})
  auto* hmissl1=new TH1F("hmissl1","MM(e' K^{+}) after #Delta#beta selection of K^{+} and p and MM(e' p K^{+}) for 2nd peak",200,0,4);

  //MM(e' p K^{+})
  auto* hmisspim1=new TH1F("hmisspim1","MM(e' p K^{+}) after #Delta#beta selection of K^{+} and p and MM(e' K^{+}) for 2nd peak",200,-0.6,0.6);
  auto* hmisspim2=new TH1F("hmisspim2","MM(e' p K^{+}) after #Delta#beta selection of K^{+} and p  and MM(e' K^{+}) for 3rd peak",200,-0.6,0.6);

  //setting axis titles
  hmissl1->GetXaxis()->SetTitle("MM(e' K^{+}) [GeV/c^{2}]");
  hmissl1->GetYaxis()->SetTitle("Counts");

  hmisspim1->GetXaxis()->SetTitle("MM(e' p K^{+})^2 [GeV/c^{2}]");
  hmisspim1->GetYaxis()->SetTitle("Counts");

  hmisspim2->GetXaxis()->SetTitle("MM(e' p K^{+})^2 [GeV/c^{2}]");
  hmisspim2->GetYaxis()->SetTitle("Counts");

  gBenchmark->Start("timer");
  int counter=0;

  for(Int_t i=0;i<files->GetEntries();i++){
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());
    //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}


    while(c12.next()==true){
      // c12.event()->getStartTime(); //hipo4
      // c12.head()->getStartTime(); //hipo3
      //Loop over all particles to see how to access detector info.

 	    for(auto& p : c12.getDetParticles()){
        //  get predefined selected information
        p->getTime();
	      p->par()->getBeta();
	      p->par()->getP();
	      p->par()->getCharge();
	      p->getDetEnergy();
	      p->getDeltaEnergy();


        // get any detector information (if exists for this particle)
        // there should be a get function for any entry in the bank
        switch(p->getRegion()) {
          case FD :
	        p->cal(PCAL)->getEnergy();
	        p->cal(ECIN)->getEnergy();
	        p->cal(ECOUT)->getEnergy();
	        p->sci(FTOF1A)->getEnergy();
	        p->sci(FTOF1B)->getEnergy();
	        p->sci(FTOF2)->getEnergy();
	        p->trk(DC)->getSector();
	        p->che(HTCC)->getNphe();
	        p->che(LTCC)->getNphe();
	        //trajectories
	        p->traj(LTCC)->getX();
	        // p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
	        break;
	        case FT :
      	  p->ft(FTCAL)->getEnergy();
      	  p->ft(FTHODO)->getEnergy();
      	  break;
      	  case CD:
      	  p->sci(CTOF)->getEnergy();
      	  p->sci(CND)->getEnergy();
      	  break;
        }
        //   covariance matrix (comment in to see!)
        // p->covmat()->print();
        p->cmat();
      }

      // get particles by type
      auto electrons=c12.getByID(11);
      auto kaonps=c12.getByID(321);
      auto protons=c12.getByID(2212);
      auto pims=c12.getByID(-211);


      //select correct events
      if(electrons.size()==1 && kaonps.size()==1 && protons.size()==1){

        // set the particle momentum
	      SetLorentzVector(el,electrons[0]);
	      SetLorentzVector(kp,kaonps[0]);
        SetLorentzVector(pr,protons[0]);

        //Creating TLorentzVectors to fill histograms
	      Double_t DeltaBetaKp=kp.Beta()-kaonps[0]->par()->getBeta();
        Double_t DeltaBetapr=pr.Beta()-protons[0]->par()->getBeta();
        TLorentzVector misspim=beam+target-el-kp-pr;
        TLorentzVector missl=beam+target-el-kp;

        //Filling histograms
        //cuts on delta beta of kaons and kaon momentum
	      if(abs(DeltaBetaKp)<0.02 && kp.P()>0.2 && abs(DeltaBetapr)<0.01 && pr.P()>0.2){

            //Cut on MM(e' K^{+})
            if(missl.M()>1.4 && missl.M()<1.65){
              hmisspim1->Fill(misspim.M2());
            }
            if(missl.M()>1.7 && missl.M()<2){
              hmisspim2->Fill(misspim.M2());
            }
            //Cut on MM(e' K^{+} p)
            if(misspim.M2()>0.2 && misspim.M2()<0.3){
              hmissl1->Fill(missl.M());
            }
          }
        }
      }

      counter++;
     }


   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");

   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
   fileOutput1.Write();
 }
