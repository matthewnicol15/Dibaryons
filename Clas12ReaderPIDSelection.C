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

void Clas12ReaderPIDSelection(){
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

  TFile fileOutput1("skim4.root","recreate");

  //creating histograms for missing mass, invariant mass and delta beta
  //MM(e' K^{+})
  auto* hmissl1=new TH1F("missl1","MM(e' K^{+})",200,0,4);
  auto* hmissl2=new TH1F("hmissl2","MM(e' K^{+}) after #Delta#beta selection of K^{+}",200,0,4);
  auto* hmissl3=new TH1F("hmissl3","MM(e' K^{+}) after #Delta#beta selection of K^{+} and p",200,0,4);
  auto* hmissl4=new TH1F("hmissl4","MM(e' K^{+}) after #Delta#beta selection of K^{+} and p and MM(e' p K^{+})",200,0,4);

  //MM(e' K^{+}) against momentum of individual particles
  auto* hmisslp1=new TH2F("hmisslp1","MM(e' K^{+}) after #Delta#beta selection of K^{+} and p and MM(e' p K^{+}) against K^{+} momentum",200,0,12,200,0,4);
  auto* hmisslp2=new TH2F("hmisslp2","MM(e' K^{+}) after #Delta#beta selection of K^{+} and p and MM(e' p K^{+}) against p momentum",200,0,12,200,0,4);
  auto* hmisslp3=new TH2F("hmisslp3","MM(e' K^{+}) after #Delta#beta selection of K^{+} and p and MM(e' p K^{+}) against e' momentum",200,0,12,200,0,4);

  //MM(e' p K^{+})^2 against MM(e'K^{+})^2
  auto* hmiss=new TH2F("hmiss","MM(e' p K^{+})^2 against MM(e'K^{+})",200,-2,2,200,-1,4);

  //MM(e' p K^{+})
  auto* hmisspim1=new TH1F("hmisspim1","MM(e' p K^{+})",200,-0.6,0.6);
  auto* hmisspim2=new TH1F("hmisspim2","MM(e' p K^{+}) after #Delta#beta selection of K^{+}",200,-0.6,0.6);
  auto* hmisspim3=new TH1F("hmisspim3","MM(e' p K^{+}) after #Delta#beta selection of K^{+} and p",200,-0.6,0.6);
  auto* hmisspim4=new TH1F("hmisspim4","MM(e' p K^{+}) after #Delta#beta selection of K^{+} and p  and MM(e' K^{+})",200,-0.6,0.6);

  //invariant mass p pi^-
  auto* himl=new TH1F("himl","M(p #pi^{-}) after #Delta#beta selection of K^{+} and p",200,0,4);

  //delta beta plots
  auto* hBetaKp=new TH2F("hBetaKp","#Delta#beta of K^{+}(post PID)",100,0,3,200,-0.02,0.02);
  auto* hBetapr=new TH2F("hBetapr","#Delta#beta of p(post PID)",200,0,3,200,-0.02,0.02);
  auto* hBetapr2=new TH2F("hBetapr2","#Delta#beta of p(post PID) after #Delta#beta selection of K^{+}",200,0,3,200,-0.02,0.02);


  //setting axis titles
  hmissl1->GetXaxis()->SetTitle("MM(e' K^{+}) [GeV/c^{2}]");
  hmissl1->GetYaxis()->SetTitle("Counts");

  hmissl2->GetXaxis()->SetTitle("MM(e' K^{+}) [GeV/c^{2}]");
  hmissl2->GetYaxis()->SetTitle("Counts");

  hmissl3->GetXaxis()->SetTitle("p(K^{+}) [GeV/c^{2}]");
  hmissl3->GetYaxis()->SetTitle("MM(e' K^{+}) [GeV/c^{2}]");

  hmisslp1->GetXaxis()->SetTitle("MM(e' p K^{+}) [GeV/c^{2}]");
  hmisslp1->GetYaxis()->SetTitle("MM(e' K^{+}) [GeV/c^{2}]");

  hmisslp2->GetXaxis()->SetTitle("MM(e' p K^{+}) [GeV/c^{2}]");
  hmisslp2->GetYaxis()->SetTitle("MM(e' K^{+}) [GeV/c^{2}]");

  hmisslp3->GetXaxis()->SetTitle("MM(e' p K^{+}) [GeV/c^{2}]");
  hmisslp3->GetYaxis()->SetTitle("MM(e' K^{+}) [GeV/c^{2}]");

  hmisspim1->GetXaxis()->SetTitle("MM(e' p K^{+})^2 [GeV/c^{2}]");
  hmisspim1->GetYaxis()->SetTitle("Counts");

  hmisspim2->GetXaxis()->SetTitle("MM(e' p K^{+})^2 [GeV/c^{2}]");
  hmisspim2->GetYaxis()->SetTitle("Counts");

  hmisspim3->GetXaxis()->SetTitle("MM(e' p K^{+})^2 [GeV/c^{2}]");
  hmisspim3->GetYaxis()->SetTitle("Counts");

  hmisspim4->GetXaxis()->SetTitle("MM(e' p K^{+})^2 [GeV/c^{2}]");
  hmisspim4->GetYaxis()->SetTitle("Counts");

  himl->GetXaxis()->SetTitle("M (p #pi^{-})[GeV/c^{2}]");
  himl->GetYaxis()->SetTitle("counts");

  hBetaKp->GetXaxis()->SetTitle("P [GeV/c]");
  hBetaKp->GetYaxis()->SetTitle("#Delta#beta");

  hBetapr->GetXaxis()->SetTitle("P [GeV/c]");
  hBetapr->GetYaxis()->SetTitle("#Delta#beta");

  hBetapr2->GetXaxis()->SetTitle("P [GeV/c]");
  hBetapr2->GetYaxis()->SetTitle("#Delta#beta");

  gBenchmark->Start("timer");
  int counter=0;


  for(Int_t i=0;i<files->GetEntries();i++){
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());
    //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}


    while(c12.next()==true){


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
        TLorentzVector lambda = pr+pim;

        //Filling histograms
        //PID selection cut
        hmissl1->Fill(missl.M());
        hmisspim1->Fill(misspim.M2());
        hBetaKp->Fill(kp.P(),DeltaBetaKp);
        hBetapr->Fill(pr.P(),DeltaBetapr);
        hmiss->Fill(misspim.M2(),missl.M());

        //cuts on delta beta of kaons and kaon momentum
	      if(abs(DeltaBetaKp)<0.02 && kp.P()>0.2){
          hmissl2->Fill(missl.M());
          hmisspim2->Fill(misspim.M2());
          hBetapr2->Fill(pr.P(),DeltaBetapr);

          //Cuts on delta beta of protons and proton momentum
          if(abs(DeltaBetapr)<0.01 && pr.P()>0.2){
            hmissl3->Fill(missl.M());
            hmisspim3->Fill(misspim.M2());

            //Cut on MM(e' K^{+})
            if(missl.M()>1 && missl.M()<1.3){
              hmisspim4->Fill(misspim.M2());
            }

            //Cut on MM(e' K^{+} p)
            if(misspim.M2()>-0.05 && misspim.M2()<0.1){
              hmissl4->Fill(missl.M());
              hmisslp1->Fill(kp.P(),missl.M());
              hmisslp2->Fill(pr.P(),missl.M());
              hmisslp3->Fill(el.P(),missl.M());
            }
          }
        }
      }

      counter++;
     }
   }


   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");

   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
   fileOutput1.Write();
 }
