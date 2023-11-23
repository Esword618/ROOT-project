#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TString.h>

using namespace std;

TString outputFile = "a.root";
TString inputFile = "i.root";

int evtID;
int MuMult;
int PDG[25];
double MuInitPosx[25];
double MuInitPosy[25];
double MuInitPosz[25];
double MuInitPx[25];
double MuInitPy[25];
double MuInitPz[25];
double MuInitKine[25];
double TrackLengthInRock[25];
double TrackLengthInVetoWater[25];
double TrackLengthInCDWater[25];
double TrackLengthInAcrylic[25];
double TrackLengthInSteel[25];
double TrackLengthInScint[25];
double ELossInRock[25];
double ELossInVetoWater[25];
double ELossInCDWater[25];
double ELossInAcrylic[25];
double LossInSteel[25];
double ELossInScint[25];
double MuExitPosx[25];
double MuExitPosy[25];
double MuExitPosz[25];
int mustpMaterial[25];


int evtID_n;
int stepNumber;
int fnPDG[10000];
double fnPositionx[10000];
double fnPositiony[10000];
double fnPositionz[10000];
double fnGlobalTime[10000];
double fnTotalEnergy[10000];
double fnKineticEnergy[10000];
double fnEnergyDeposit[10000];
double fnQEnergyDeposit[10000];
double fnCaptime[10000];



void ana(TString fileName, TString outName); 

int main(int argc, char** argv) {

  bool noOutput = true;
  bool noInputPath = true;
  for(int i=1; i<argc; i++) {
				if(strcmp(argv[i],"-o") == 0) {
						outputFile = argv[i+1];
						noOutput = false;
				}
				if(strcmp(argv[i],"-i") == 0) {
						inputFile = argv[i+1];
						noInputPath = false;
				}
  }

  if(noOutput) cout << "Info: Use default OutputFile!" << endl;
  if(noInputPath) cout << "Info: Use default InputDataPath!" << endl;

  cout << "**********   Arguments   ***********" << endl ;
  cout << "OutputFile       | "  << outputFile << endl ;
  cout << "InputDataPath    | "  << inputFile  << endl ;
  cout << "**********      End      ***********\n" << endl ;

  //TFile *f_hists = new TFile(outputFile, "recreate");
  ana(inputFile, outputFile);
}

void ana(TString fileName, TString outName) {
  
  TFile* f = new TFile(fileName);
  int nentries;
  
  TTree* muT = (TTree*)f->Get("mu");
  muT->SetBranchAddress("evtID", &evtID);
  muT->SetBranchAddress("MuMult", &MuMult);
  muT->SetBranchAddress("PDG",PDG);
  muT->SetBranchAddress("MuInitPosx",MuInitPosx);
  muT->SetBranchAddress("MuInitPosy",MuInitPosy);
  muT->SetBranchAddress("MuInitPosz",MuInitPosz);
  muT->SetBranchAddress("MuInitPx",MuInitPx);
  muT->SetBranchAddress("MuInitPy",MuInitPy);
  muT->SetBranchAddress("MuInitPz",MuInitPz);
  muT->SetBranchAddress("MuInitKine",MuInitKine);
  muT->SetBranchAddress("TrackLengthInRock",TrackLengthInRock);
  muT->SetBranchAddress("TrackLengthInVetoWater",TrackLengthInVetoWater);
  muT->SetBranchAddress("TrackLengthInCDWater",TrackLengthInCDWater);
  muT->SetBranchAddress("TrackLengthInAcrylic",TrackLengthInAcrylic);
  muT->SetBranchAddress("TrackLengthInSteel",TrackLengthInSteel);
  muT->SetBranchAddress("TrackLengthInScint",TrackLengthInScint);
  muT->SetBranchAddress("ELossInRock",ELossInRock);
  muT->SetBranchAddress("ELossInVetoWater",ELossInVetoWater);
  muT->SetBranchAddress("ELossInCDWater",ELossInCDWater);
  muT->SetBranchAddress("ELossInAcrylic",ELossInAcrylic);
  //t->SetBranchAddress("ELossInSteel",LossInSteel);
  muT->SetBranchAddress("ELossInScint",ELossInScint);
  muT->SetBranchAddress("MuExitPosx",MuExitPosx);
  muT->SetBranchAddress("MuExitPosy",MuExitPosy);
  muT->SetBranchAddress("MuExitPosz",MuExitPosz);
  muT->SetBranchAddress("mustpMaterial",mustpMaterial);

  TTree* fastT = (TTree*)f->Get("mufastn");
  fastT->SetBranchAddress("evtID", &evtID_n);
  fastT->SetBranchAddress("stepNumber", &stepNumber);
  fastT->SetBranchAddress("fnPDG", fnPDG);
  fastT->SetBranchAddress("fnPositionx", fnPositionx);
  fastT->SetBranchAddress("fnPositiony", fnPositiony);
  fastT->SetBranchAddress("fnPositionz", fnPositionz);
  fastT->SetBranchAddress("fnGlobalTime", fnGlobalTime);
  fastT->SetBranchAddress("fnTotalEnergy", fnTotalEnergy);
  fastT->SetBranchAddress("fnEnergyDeposit", fnEnergyDeposit);
  fastT->SetBranchAddress("fnQEnergyDeposit", fnQEnergyDeposit);
  fastT->SetBranchAddress("fnKineticEnergy", fnKineticEnergy);
  fastT->SetBranchAddress("fnCaptime", fnCaptime);
 
  TH1F*muTL = new TH1F("muTL", "muTL", 40000, 0, 40); //--------------> mu track length --> m
  TH1F*muKE = new TH1F("muKE", "muKE", 10000, 0, 10000); //--------------> mu kinetic energy --> GeV
  TH1F*mudepE = new TH1F("mudepE", "mudepE", 4000, 0, 4000); //--------------> mu kinetic energy --> GeV

  TFile*f_merged = new TFile(outName, "recreate");
  TTree*t_merged = new TTree("Merged", "production tree");
  //int fileId;
  //int eventId;

  const int Nmax = 1000;
  //------------ neutron info
  int Nsignals;
  int pdg[Nmax];
  double x[Nmax];
  double y[Nmax];
  double z[Nmax];
  double t[Nmax];
  double Edep[Nmax];
  double QE[Nmax];
  double tcap[Nmax];
  double KE[Nmax];
  //t_merged->Branch("fileId", &fileId, "fileId/I");
  //t_merged->Branch("eventId", &eventId, "eventId/I");
  t_merged->Branch("Nsignals", &Nsignals, "Nsignals/I");
  t_merged->Branch("pdg", pdg, "pdg[Nsignals]/I");
  t_merged->Branch("x", x, "x[Nsignals]/D");
  t_merged->Branch("y", y, "y[Nsignals]/D");
  t_merged->Branch("z", z, "z[Nsignals]/D");
  t_merged->Branch("t", t, "t[Nsignals]/D");
  t_merged->Branch("Edep", Edep, "Edep[Nsignals]/D");
  t_merged->Branch("QE", QE, "QE[Nsignals]/D");
  t_merged->Branch("tcap", tcap, "tcap[Nsignals]/D");
  t_merged->Branch("KE", KE, "KE[Nsignals]/D");
  


  //-----------> ana Mu 
  int ne = (int)muT->GetEntries();
  for(int jj=0; jj<ne; jj++) {
    muT->GetEntry(jj);
    for(int mm=0; mm<MuMult; mm++) {
      if(TrackLengthInScint[mm]>0) {
        muTL->Fill(TrackLengthInScint[mm]*1.e-3);
        muKE->Fill(MuInitKine[mm]*1.e-3);
        mudepE->Fill(-1.*ELossInScint[mm]*1.e-3);
      }
    }
  }

  nentries = (int)fastT->GetEntries();
  for(int ii=0; ii<nentries; ii++) {
    fastT->GetEntry(ii);
    //cout << "-------> EventID : " << evtID_n << endl;
    if(stepNumber>=1) {
      bool tagMuFlag = false;
      bool neuFlag = false;
      muT->GetEntry(evtID_n);
      for(int mm=0; mm<MuMult; mm++) {
				double trackInScint = TrackLengthInScint[mm];
				if(trackInScint>0) {
            tagMuFlag = true; 
				} else {
            tagMuFlag = false;
        }        
      }
//      if(rockMuFlag) {
//        cout << "-------> EventID : " << evtID_n << endl;
//      }
      for(int stepN=0; stepN<stepNumber; stepN++) {
        
        if(fnCaptime[stepN] > 0 && fnPDG[stepN] == 2112) {
          neuFlag = true;
//          cout << "-------> EventID : " << evtID_n << "------> fastn ID" << ii <<endl;
        }
      }
      
      if(tagMuFlag && neuFlag) {
        cout << "-------> EventID : " << evtID_n << endl;
        int t_PDG = 0;
        double t_x=0.,t_y=0.,t_z=0.;
        double t_t=0.;
        double t_edep=0.;
        double t_qe=0.;
        double t_tcap =0.;
        double t_ke =0;
        double t_win = 300.;
        Nsignals = 0;
        int steps=0;
        for(int j=0; j<Nmax; j++) {
          pdg[j] =0;
          x[j]=0;
          y[j]=0;
          z[j]=0;
          t[j]=0;
          Edep[j]=0;
          QE[j]=0;
          tcap[j]=0;
          KE[j]=0;
        } 
        for(int stepN=0; stepN<stepNumber; stepN++) {
          if(stepN==0) {
            t_t = fnGlobalTime[stepN];
            t_x = fnPositionx[stepN]*fnEnergyDeposit[stepN];
            t_y = fnPositiony[stepN]*fnEnergyDeposit[stepN];
            t_z = fnPositionz[stepN]*fnEnergyDeposit[stepN];
            t_edep = fnEnergyDeposit[stepN];
            t_qe = fnQEnergyDeposit[stepN];
            t_ke = fnKineticEnergy[stepN];
            steps = steps+1;
            t_PDG = 0;
            t_tcap = 0;
            if(fnCaptime[stepN] > 0 && fnPDG[stepN] == 2112) {
              t_PDG = 2112;
              t_tcap = fnCaptime[stepN];
            } 
          } else {
            double dt = fnGlobalTime[stepN] - t_t;
            if(dt<0) cout << " Warning  ::: Call Jie : please check input root file"<< endl;
            if(dt>0 && dt < t_win) {
              t_x += fnPositionx[stepN]*fnEnergyDeposit[stepN];
              t_y += fnPositiony[stepN]*fnEnergyDeposit[stepN];
              t_z += fnPositionz[stepN]*fnEnergyDeposit[stepN];
              t_edep += fnEnergyDeposit[stepN];
              t_qe += fnQEnergyDeposit[stepN];
              t_ke += fnKineticEnergy[stepN];
              steps = steps+1;
              if(fnCaptime[stepN] > 0 && fnPDG[stepN] == 2112) {
                if(t_PDG == 2112) cout << "Warning  ::: Call Jie : please check input root file for double neutrons"<< endl;
                t_PDG = 2112;
                t_tcap = fnCaptime[stepN];
              } 
            }else if(dt >= t_win) {
              pdg[Nsignals] = t_PDG;
              x[Nsignals] = t_x/t_edep;
              y[Nsignals] = t_y/t_edep;
              z[Nsignals] = t_z/t_edep;
              Edep[Nsignals] = t_edep;
              QE[Nsignals] = t_qe;
              t[Nsignals] = t_t*1.e-3;
              tcap[Nsignals] = t_tcap;
              KE[Nsignals] = t_ke/steps;
              steps=0;
              Nsignals ++ ;
              
              t_PDG = 0;
              t_tcap = 0;
              t_t = fnGlobalTime[stepN];
              t_x = fnPositionx[stepN]*fnEnergyDeposit[stepN];
              t_y = fnPositiony[stepN]*fnEnergyDeposit[stepN];
              t_z = fnPositionz[stepN]*fnEnergyDeposit[stepN];
              t_edep = fnEnergyDeposit[stepN];
              t_qe = fnQEnergyDeposit[stepN];
              t_ke = fnKineticEnergy[stepN];
              steps=steps+1;
              if(fnCaptime[stepN] > 0 && fnPDG[stepN] == 2112) {
                t_PDG = 2112;
                t_tcap = fnCaptime[stepN];
              } 
            }
          }
          if(stepN == stepNumber-1) {
            pdg[Nsignals] = t_PDG;
            x[Nsignals] = t_x/t_edep;
            y[Nsignals] = t_y/t_edep;
            z[Nsignals] = t_z/t_edep;
            Edep[Nsignals] = t_edep;
            QE[Nsignals] = t_qe;
            t[Nsignals] = t_t*1.e-3;
            tcap[Nsignals] = t_tcap;
            KE[Nsignals] = t_ke/steps;
            steps=0;
            Nsignals ++ ;
          }
        }
        t_merged->Fill();
      }
    }
  }
  f_merged->cd();
  muTL->Write();
  muKE->Write();
  mudepE->Write();
  t_merged->Write();
  f_merged->Close();
}

