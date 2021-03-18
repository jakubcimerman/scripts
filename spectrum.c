#include <iostream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>

void spectrum(const char* direct, int Npart_min, int Npart_max)
{
  cout << endl << endl << "Calculating pT spectrum of pi+, K+ and protons..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Borders of eta
  const double ptMin = 0.2;
  const double ptMax = 2.0;
  const double yCut = 0.1;
  const int nBins = 36;
  const double pi = 3.1415926535897932;
  const double dpt = (ptMax-ptMin)/nBins;

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
  const int NP = 60000;

  Int_t id[NP];
  Float_t px [NP], py[NP], pz[NP], E[NP];
  Int_t npart, Nparticipants;

  int nevents = tree->GetEntries();
  cout << "Total number of events in directory: " << nevents << endl;
  int centrality_events = 0;
  tree->SetBranchAddress("px",&px[0]);
  tree->SetBranchAddress("py",&py[0]);
  tree->SetBranchAddress("pz",&pz[0]);
  tree->SetBranchAddress("E",&E[0]);
  tree->SetBranchAddress("id",&id[0]);
  tree->SetBranchAddress("npart",&npart);
  tree->SetBranchAddress("Nparticipants",&Nparticipants);

  double dNpi[nBins]={0.0}, dNK[nBins]={0.0}, dNp[nBins]={0.0};
  for(int iev=0; iev<nevents; iev++)
  {
    if ((iev)%((int)nevents/20) == 0) cout << (int)100*iev/nevents << "%" << endl;
    tree->GetEntry(iev);
    if (Nparticipants > Npart_min && Nparticipants < Npart_max)
    {
      centrality_events++;
      for(int i=0; i<npart; i++)
      {
        float rap = 0.5*log((E[i]+pz[i])/(E[i]-pz[i]));
        float pt = sqrt(px[i]*px[i]+py[i]*py[i]);
        if (abs(rap) < yCut && pt > ptMin && pt < ptMax) {
          int ptBin = (pt-ptMin)/dpt ;
          if (id[i] == 211) dNpi[ptBin]++;
          if (id[i] == 321) dNK[ptBin]++;
          if (id[i] == 2212) dNp[ptBin]++;
        }
      }
    }
  }

  cout << "Number of events in chosen centrality interval: " << centrality_events << endl;
  nevents = centrality_events;

  // Write results into the text file (append)
  ofstream fout;
  fout.open("spectrum.dat", ofstream::app);
  fout << direct << endl;
  for (int i = 0; i < nBins; i++)
  {
    dNpi[i] /= (nevents*2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNK[i] /= (nevents*2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNp[i] /= (nevents*2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    fout << ptMin + (i+0.5)*dpt  << "\t" << dNpi[i] << "\t" << dNK[i] << "\t" << dNp[i] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'spectrum.dat'" << endl;

}
