#include <iostream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>

void dndeta(const char* direct, int Npart_min, int Npart_max);

void dndeta(const char* direct, int Npart_min, int Npart_max)
{
  // Borders of eta
  const double etaMin = -5.0;
  const double etaMax =  5.0;
  const int nBins = 50;
  const double deta = (etaMax-etaMin)/nBins;

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
    cout << "[" << fname.Data() << "]";
    if(strstr(fname,".root"))
      tree->Add(dirname+fname);
  }

  // Maximum number of particles - length of arrays
  const int NP = 60000;

  Float_t px [NP], py[NP], pz[NP];
  Short_t ele[NP];
  Int_t npart, Nparticipants;

  int nevents = tree->GetEntries();
  int centrality_events = 0;
  tree->SetBranchAddress("px",&px[0]);
  tree->SetBranchAddress("py",&py[0]);
  tree->SetBranchAddress("pz",&pz[0]);
  tree->SetBranchAddress("ele",&ele[0]);
  tree->SetBranchAddress("npart",&npart);
  tree->SetBranchAddress("Nparticipants",&Nparticipants);

  double dN[nBins]={0.0};
  for(int iev=0; iev<nevents; iev++)
  {
    if (iev%((int)nevents/100) == 0) cout << (int)100*iev/nevents << "%" << endl;
    tree->GetEntry(iev);
    if (Nparticipants > Npart_min && Nparticipants < Npart_max)
    {  
      centrality_events++;
      for(int i=0; i<npart; i++)
      {
        const float eta = 0.5*log((sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])+pz[i])/(sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-pz[i]));
        if (ele[i]!=0 && eta < etaMax && eta > etaMin)
        {
          int etaBin = (eta-etaMin)/deta;
          dN[etaBin]++;
        }
      }
    }
  }

  cout << centrality_events << endl;

  // Write results into the text file (append)
  ofstream fout;
  fout.open("dndeta.dat", ofstream::app);
  fout << direct << endl;
  for (int i = 0; i < nBins; i++)
  {
    dN[i] /= (centrality_events*deta);
    fout << etaMin + (i+0.5)*deta  << "\t" << dN[i] << endl;
  }
  fout << endl;
  fout.close();

}
