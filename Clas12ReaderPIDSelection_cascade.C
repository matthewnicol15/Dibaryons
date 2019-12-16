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
#include <vector>

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void Clas12ReaderPIDSelection_cascade(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //Choose the data files that you would like to run the macro over
   TString inputFile("/work/clas12/rg-a/trains/v16_v2/skim4_inclusive/skim4_50*.hipo");



   cout<<"Analysing hipo file "<<inputFile<<endl;

   TChain fake("hipo");
   fake.Add(inputFile.Data());

   //get the hipo data
   //   reader.open(inputFile.Data());
   auto files=fake.GetListOfFiles();

   //particle identification
   auto db=TDatabasePDG::Instance();
   TLorentzVector beam(0,0,10.6,10.6);
   TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
   TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
   TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
   TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
   TLorentzVector kp1(0,0,0,db->GetParticle(321)->Mass());
   TLorentzVector kp2(0,0,0,db->GetParticle(321)->Mass());
   Int_t eventno1=0;
   Int_t eventno2=0;
   Int_t eventno3=0;
   Int_t eventno4=0;
   Int_t eventno5=0;
   Int_t eventno6=0;
   Int_t eventno7=0;

   TFile fileOutput1("skim4_50_cascade.root","recreate");

   //creating histograms
   //histograms for finding missing mass (cascade)
   auto* hmissc1=new TH1F("missc1","MM(e' K^{+} K^{+}) after PID selection",200,0,3.5);
   auto* hmissc2=new TH1F("missc2","MM(e' K^{+} K^{+}) after PID and #Delta#beta selection of both K^{+}",200,0,5);
   auto* hmissc3=new TH1F("missc3","MM(e' K^{+} K^{+}) after PID and #Delta#beta selection of both K^{+} and proton",200,0,3.5);
   auto* hmissc4=new TH1F("missc4","MM(e' K^{+} K^{+}) after PID, #Delta#beta selection of both K^{+} and proton and MM(e' K^{+} K^{+} #pi^{-}) cut",200,0,3.5);

   //histograms for finding missing mass (lambda^{0})
   auto* hmissl01=new TH1F("hmissl01","MM(e' K^{+} K^{+} #pi^{-}) after PID selection",200,0,5);
   auto* hmissl02=new TH1F("hmissl02","MM(e' K^{+} K^{+} #pi^{-}) after PID and #Delta#beta selection of both K^{+}",200,0,5);
   auto* hmissl03=new TH1F("hmissl03","MM(e' K^{+} K^{+} #pi^{-}) after PID and #Delta#beta selection of both K^{+} and proton",200,0,5);
   auto* hmissl04=new TH1F("hmissl04","MM(e' K^{+} K^{+} #pi^{-}) after PID, #Delta#beta selection of both K^{+} and proton and MM(e' K^{+} K^{+}) cut",200,0,5);

   // histograms for missing mass (cascade) against missing mass (lambda^{0}) after deltabeta selection of both K^{+}
   auto* hmissc1l01=new TH2F("hmissc1l01","MM(e' K^{+} K^{+}) against MM(e' K^{+} K^{+} #pi^{-} after PID selection)",200,0,5,200,0,5);
   auto* hmissc2l02=new TH2F("hmissc2l02","MM(e' K^{+} K^{+}) against MM(e' K^{+} K^{+} #pi^{-}) after PID and #Delta#beta selection of both K^{+} and proton",200,0,5,200,0,5);

   //deltabeta of K^{+}
   auto* hBetakp1=new TH2F("hBetakp1","#Delta#beta of K^{+}(post PID)",100,0,12,200,-0.15,0.15);
   auto* hBetakp2=new TH2F("hBetakp2","#Delta#beta of K^{+}(post PID)",100,0,12,200,-0.15,0.15);

   //deltabeta of p
   auto* hBetapr=new TH2F("hBetapr","#Delta#beta of p(post PID)",100,0,12,200,-0.15,0.15);

   //kaon1 momentum against kaon 2 momentum
   auto* hkaons=new TH2F("hkaons","momentum of Kaon[0] against momentum of Kaon[1]",200,0,12,200,0,12);


   //setting axis titles
   hmissc1->GetXaxis()->SetTitle("MM(e' K^{+} K^{+}) [GeV/c^{2}]");
   hmissc1->GetYaxis()->SetTitle("Counts");

   hmissc2->GetXaxis()->SetTitle("MM(e' K^{+} K^{+}) [GeV/c^{2}]");
   hmissc2->GetYaxis()->SetTitle("Counts");

   hmissc3->GetXaxis()->SetTitle("MM(e' K^{+} K^{+}) [GeV/c^{2}]");
   hmissc3->GetYaxis()->SetTitle("Counts");

   hmissl01->GetXaxis()->SetTitle("MM(e' K^{+} K^{+} #pi^{-}) [GeV/c^{2}]");
   hmissl01->GetYaxis()->SetTitle("Counts");

   hmissl02->GetXaxis()->SetTitle("MM(e' K^{+} K^{+} #pi^{-}) [GeV/c^{2}]");
   hmissl02->GetYaxis()->SetTitle("Counts");

   hmissc2l02->GetXaxis()->SetTitle("MM(e' K^{+} K^{+} #pi^{-}) [GeV/c^{2}]");
   hmissc2l02->GetYaxis()->SetTitle("MM(e' K^{+} K^{+}) [GeV/c^{2}]");

   hBetakp1->GetXaxis()->SetTitle("P [GeV/c]");
   hBetakp1->GetYaxis()->SetTitle("#Delta#beta");

   hBetakp2->GetXaxis()->SetTitle("P [GeV/c]");
   hBetakp2->GetYaxis()->SetTitle("#Delta#beta");


   gBenchmark->Start("timer");
   int counter=0;

   for(Int_t i=0;i<files->GetEntries();i++){
     //create the event reader
     clas12reader c12(files->At(i)->GetTitle());


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

       eventno1 = eventno1 + 1;

       // get particles by type
       auto electrons=c12.getByID(11);
       auto kaonps=c12.getByID(321);
       auto protons=c12.getByID(2212);
       auto pims=c12.getByID(-211);

       //PID selection
       if(electrons.size()==1){
         eventno2 = eventno2 + 1;
         if(kaonps.size()>0){
           eventno3 = eventno3 + 1;
           if(kaonps.size()==2){
             eventno4 = eventno4 + 1;
             if(protons.size()==1){
               eventno5 = eventno5 + 1;

            	 // set the particle momentum
            	 SetLorentzVector(el,electrons[0]);
            	 SetLorentzVector(kp1,kaonps[0]);
               SetLorentzVector(kp2,kaonps[1]);
               SetLorentzVector(pr,protons[0]);


               //calculation of delta beta for kp1 and kp2
	              Double_t DeltaBetakp1=kp1.Beta()-kaonps[0]->par()->getBeta();
                Double_t DeltaBetakp2=kp2.Beta()-kaonps[1]->par()->getBeta();
                Double_t DeltaBetapr=pr.Beta()-protons[0]->par()->getBeta();

                //MM(e' K^{+} K^{+}), looking for cascade^[-]
                TLorentzVector misscasc=beam+target-el-kp1-kp2;
                TLorentzVector missl0=beam+target-el-kp1-kp2-pim;

                //Filling histograms
                //PID selection only
                hkaons->Fill(kp1.P(),kp2.P());
                hmissc1->Fill(misscasc.M());
                hBetakp1->Fill(kp1.P(),DeltaBetakp1);
                hBetakp2->Fill(kp2.P(),DeltaBetakp2);
                hBetapr->Fill(pr.P(),DeltaBetapr);
                hmissl01->Fill(missl0.M());
                hmissc1l01->Fill(missl0.M(),misscasc.M());

                //Performing deltabeta of kaonps cuts
                if(abs(DeltaBetakp1)<0.01 && abs(DeltaBetakp2)<0.01 && kp1.P()>0.2 && kp2.P()>0.2){
                  eventno6 = eventno6 + 1;
                  hmissc2->Fill(misscasc.M());
                  hmissl02->Fill(missl0.M());

                  //Performing deltabeta of proton cuts
                  if(abs(DeltaBetapr)<0.01 && abs(DeltaBetapr)<0.01 && pr.P()>0.2){
                    eventno7 = eventno7 + 1;
                    hmissc2l02->Fill(missl0.M(),misscasc.M());
                    hmissl03->Fill(missl0.M());
                    hmissc3->Fill(misscasc.M());


                    //Performing lambda0 missing mass cuts
                    if(missl0.M()>1 && missl0.M()<1.3){
                      hmissc4->Fill(misscasc.M());
                    }
                    if(misscasc.M()>1.1 && misscasc.M()<1.5){
                      hmissl04->Fill(missl0.M());
                    }
                  }
                }
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
   // TCanvas *can1=new TCanvas("can1","MM(e' K^{+} K^{+})", 600, 600);
   // hmissc1->Draw();

fileOutput1.Write();
}
