#include <iostream>
#include <fstream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>

void decorrelation(char* direct, double order, int eventStep, double etaRefMin, double etaRefMax, double etaTestMin, double etaTestMax, int Npart_min, int Npart_max, int isSym)
{
  cout << endl << endl << "Calculating decorrelation..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Cuts
  const double ptMinCut = 0.4;
  const double ptMaxCut = 4.0;
  const int nBins = 5;
  const double deta = (etaTestMax-etaTestMin)/nBins;

  // Adding root files from directory into a tree
  TString dirname = direct;
  TChain *tree = new TChain("treefin");
  TSystemDirectory dir ("rootfiles",dirname);
  TList *files = dir.GetListOfFiles();
  TSystemFile *file;
  TIter next(files);
  while ((file=(TSystemFile*)next()))
  {
    TString fname = file->GetName();
    if(strstr(fname,".root"))
      tree->Add(dirname+fname);
  }

  // Maximum number of particles - length of arrays
  const int NP = 50000;

  Float_t px [NP], py[NP], pz[NP];
  Short_t ele[NP];
  Int_t npart, Nparticipants;
  int nevents = tree->GetEntries();
  cout << "Total number of events in directory: " << nevents << endl;
  int centrality_events = 0; // variable that counts number of events in given centrality bin
  tree->SetBranchAddress("px",&px[0]);
  tree->SetBranchAddress("py",&py[0]);
  tree->SetBranchAddress("pz",&pz[0]);
  tree->SetBranchAddress("ele",&ele[0]);
  tree->SetBranchAddress("npart",&npart);
  tree->SetBranchAddress("Nparticipants",&Nparticipants);

  double rn[nBins]={0.0}, rn_[nBins]={0.0};
  double rnErr[nBins]={0.0}, rnErr_[nBins]={0.0};
  double rnNom[nBins]={0.0}, rnDenom[nBins]={0.0};
  double rnNom_[nBins]={0.0}, rnDenom_[nBins]={0.0};
  double rnNomErr[nBins]={0.0}, rnDenomErr[nBins]={0.0};
  double rnNomErr_[nBins]={0.0}, rnDenomErr_[nBins]={0.0};
  double rnNomSD[nBins]={0.0}, rnDenomSD[nBins]={0.0};
  double rnNomSD_[nBins]={0.0}, rnDenomSD_[nBins]={0.0};

  for(int iev=0; iev<nevents; iev+= eventStep)
  {
    double QxRef = 0.0, QyRef = 0.0;
    double QxRef_ = 0.0, QyRef_ = 0.0;
    double Qx[nBins] = {0.0}, Qy[nBins] = {0.0};
    double Qx_[nBins] = {0.0}, Qy_[nBins] = {0.0};
    int NRef = 0, NRef_ = 0, NTest[nBins] = {0}, NTest_[nBins] = {0};

    for(int k = 0; k < eventStep; k++)
    {
      tree->GetEntry(iev+k);
      if ((iev+k)%((int)nevents/20) == 0) cout << (int)100*(iev+k+1)/nevents << "%" << endl;
      // centrality cut
      if (Nparticipants > Npart_min && Nparticipants < Npart_max)
      {
        for(int i=0; i<npart; i++)
        {
          const float pabs = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
          const float pt = sqrt(px[i]*px[i]+py[i]*py[i]);
          const float phi = atan2(py[i],px[i]);
          double eta = 0.5*log((pabs+pz[i])/(pabs-pz[i]));
          if(pt > ptMinCut && pt < ptMaxCut && ele[i] != 0)
          {
            if(eta < etaRefMax && eta > etaRefMin)
            {
              QxRef += cos(order*phi);  
              QyRef += sin(order*phi);
              NRef++;
            }
            else if (-eta < etaRefMax && -eta > etaRefMin)
            {
              QxRef_ += cos(order*phi);  
              QyRef_ += sin(order*phi);
              NRef_++;
            }
            else if (eta < etaTestMax && eta > etaTestMin)
            {
              int etaBin = (eta-etaTestMin)/deta;
              Qx[etaBin] += cos(order*phi);  
              Qy[etaBin] += sin(order*phi);
              NTest[etaBin]++;
            }
            else if (-eta < etaTestMax && -eta > etaTestMin)
            {
              int etaBin = (-eta-etaTestMin)/deta;
              Qx_[etaBin] += cos(order*phi);  
              Qy_[etaBin] += sin(order*phi);
              NTest_[etaBin]++;
            }
          }
        }
      }
    }
    if (NRef > 0)
    {
      centrality_events++;
      QxRef = (double)QxRef/NRef;
      QyRef = (double)QyRef/NRef;
    }
    if (NRef_ > 0)
    {
      QxRef_ = (double)QxRef_/NRef_;
      QyRef_ = (double)QyRef_/NRef_;
    }
    for (int i = 0; i < nBins; i++)
    {
      if (NTest[i] > 0)
      {
        Qx[i] = Qx[i]/NTest[i];
        Qy[i] = Qy[i]/NTest[i];
      }
      if (NTest_[i] > 0)
      {
        Qx_[i] = Qx_[i]/NTest_[i];
        Qy_[i] = Qy_[i]/NTest_[i];
      }
    }

    for (int i = 0; i < nBins; i++)
    {
      rnNom[i] += (Qx_[i]*QxRef+Qy_[i]*QyRef);
      rnDenom[i] += (Qx[i]*QxRef+Qy[i]*QyRef);
      rnNom_[i] += (Qx[i]*QxRef_+Qy[i]*QyRef_);
      rnDenom_[i] += (Qx_[i]*QxRef_+Qy_[i]*QyRef_);
      rnNomSD[i] += (Qx_[i]*QxRef+Qy_[i]*QyRef)*(Qx_[i]*QxRef+Qy_[i]*QyRef);
      rnDenomSD[i] += (Qx[i]*QxRef+Qy[i]*QyRef)*(Qx[i]*QxRef+Qy[i]*QyRef);
      rnNomSD_[i] += (Qx[i]*QxRef_+Qy[i]*QyRef_)*(Qx[i]*QxRef_+Qy[i]*QyRef_);
      rnDenomSD_[i] += (Qx_[i]*QxRef_+Qy_[i]*QyRef_)*(Qx_[i]*QxRef_+Qy_[i]*QyRef_);
    }
  }

  cout << "Number of oversampled events in chosen centrality interval: " << centrality_events << endl;

  for(int i=0; i<nBins; i++)
  {
    rn[i] = (rnNom[i]) / (rnDenom[i]);
    rn_[i] = (rnNom_[i]) / (rnDenom_[i]);
    rnNomErr[i] = sqrt(rnNomSD[i]/centrality_events - rnNom[i]*rnNom[i]/pow(centrality_events,2))/sqrt(centrality_events);
    rnNomErr_[i] = sqrt(rnNomSD_[i]/centrality_events - rnNom_[i]*rnNom_[i]/pow(centrality_events,2))/sqrt(centrality_events);
    rnDenomErr[i] = sqrt(rnDenomSD[i]/centrality_events - rnDenom[i]*rnDenom[i]/pow(centrality_events,2))/sqrt(centrality_events);
    rnDenomErr_[i] = sqrt(rnDenomSD_[i]/centrality_events - rnDenom_[i]*rnDenom_[i]/pow(centrality_events,2))/sqrt(centrality_events);
    rnErr[i] = rn[i] * centrality_events * sqrt(pow(rnNomErr[i]/rnNom[i],2) + pow(rnDenomErr[i]/rnDenom[i],2));
    rnErr_[i] = rn_[i] * centrality_events * sqrt(pow(rnNomErr_[i]/rnNom_[i],2) + pow(rnDenomErr_[i]/rnDenom_[i],2));
  }

  double *rapBin = new double [nBins];
  double *rnFin = new double [nBins];
  double *rnFinErr = new double [nBins];

  // Write results into the text file (append)
  ofstream fout;
  fout.open("decorrelation.dat", ofstream::app);
  fout << direct << " r_" << order << endl;
  for(int irap=0; irap<nBins; irap++)
  {
    rapBin[irap] = etaTestMin + (irap+0.5)*deta;
    if (isSym == 1)
    {
      rnFin[irap] = rn[irap];
      rnFinErr[irap] = rnErr[irap];
    }
    else if (isSym == 0)
    {
      rnFin[irap] = sqrt(rn[irap]*rn_[irap]);
      rnFinErr[irap] = 0.5 * rnFin[irap] * sqrt(pow(rnErr[irap]/rn[irap],2)+pow(rnErr_[irap]/rn_[irap],2));
    }
    fout << rapBin[irap] << "\t" << rnFin[irap] << "\t" << rnFinErr[irap] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'decorrelation.dat'" << endl;
}

