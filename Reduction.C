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

void Reduction(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //Choose the data files that you would like to run the macro over
   TString inputFile("/work/clas12/rg-a/trains/v16_v2/skim4_inclusive/skim4_5201.hipo");

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
   // TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
   // TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
   TLorentzVector kp1(0,0,0,db->GetParticle(321)->Mass());
   TLorentzVector kp2(0,0,0,db->GetParticle(321)->Mass());
   TLorentzVector kp3(0,0,0,db->GetParticle(321)->Mass());


   TFile fileOutput1("Reduction_5201.root","recreate");

   //creating histograms
   TH1F *hel1=new TH1F("hel1","momentum of electron;P [GeV/c];Counts",200,0,12);
   TH1F *hel2=new TH1F("hel2","momentum of electron;P [GeV/c];Counts",200,0,12);
   TH1F *hel3=new TH1F("hel3","momentum of electron;P [GeV/c];Counts",200,0,12);
   TH1F *hel4=new TH1F("hel4","momentum of electron;P [GeV/c];Counts",200,0,12);
   TH1F *hel5=new TH1F("hel5","momentum of electron;P [GeV/c];Counts",200,0,12);
   TH1F *hel6=new TH1F("hel6","momentum of electron;P [GeV/c];Counts",200,0,12);
   TH1F *hel7=new TH1F("hel7","momentum of electron;P [GeV/c];Counts",200,0,12);
   TH1F *hel8=new TH1F("hel8","momentum of electron;P [GeV/c];Counts",200,0,12);


   gBenchmark->Start("timer");
   int counter=0;

   for(Int_t i=0;i<files->GetEntries();i++){
     //create the event reader

     clas12reader c12(files->At(i)->GetTitle());

      while(c12.next()==true){
        auto particles = c12.getDetParticles();

        for(auto& p : particles){

	 switch(p->getRegion()) {
	 case FD :
	   break;
	 case FT :
	   break;
	 case CD:
	   break;

     // get particles by type
     auto electrons=c12.getByID(11);
     auto kaonps=c12.getByID(321);
     auto protons=c12.getByID(2212);

     // set the particle momentum
    SetLorentzVector(el,electrons[0]);
    SetLorentzVector(kp1,kaonps[0]);
    SetLorentzVector(kp2,kaonps[1]);
    SetLorentzVector(kp3,kaonps[2]);
    // SetLorentzVector(pr,protons[0]);

    Double_t DeltaBetakp1=kp1.Beta()-kaonps[0]->par()->getBeta();
    Double_t DeltaBetakp2=kp2.Beta()-kaonps[1]->par()->getBeta();
    Double_t DeltaBetakp3=kp3.Beta()-kaonps[2]->par()->getBeta();





     //PID selection
     if(electrons.size()==1){
       hel2->Fill(el.P());

       if(kaonps.size()==1){
         hel3->Fill(el.P());
         if(kp1.P()>0.2 && abs(DeltaBetakp1)<0.02){
           hel4->Fill(el.P());
         }
       }

       if(kaonps.size()==2){
           hel5->Fill(el.P());
           if(kp1.P()>0.2 && abs(DeltaBetakp1)<0.02 && kp2.P()>0.2 && abs(DeltaBetakp2)<0.02){
             hel6->Fill(el.P());
           }
         }
         if(kaonps.size()==3){
           hel7->Fill(el.P());
           if(kp1.P()>0.2 && abs(DeltaBetakp1)<0.02 && kp2.P()>0.2 && abs(DeltaBetakp2)<0.02 && kp3.P()>0.2 && abs(DeltaBetakp3)<0.02){
             hel8->Fill(el.P());
           }
         }
       }

	 }
       }



        TH2F *hratio1 = (TH2F*)hel2->Clone("hratio1");
        hratio1->Divide(hel1);

        TH2F *hratio2 = (TH2F*)hel3->Clone("hratio2");
        hratio2->Divide(hel2);

        TH2F *hratio3 = (TH2F*)hel4->Clone("hratio3");
        hratio3->Divide(hel3);

        TH2F *hratio4 = (TH2F*)hel5->Clone("hratio4");
        hratio4->Divide(hel4);

        TH2F *hratio5 = (TH2F*)hel6->Clone("hratio5");
        hratio5->Divide(hel5);

        TH2F *hratio6 = (TH2F*)hel7->Clone("hratio6");
        hratio6->Divide(hel6);

        TH2F *hratio7 = (TH2F*)hel8->Clone("hratio7");
        hratio7->Divide(hel7);


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
