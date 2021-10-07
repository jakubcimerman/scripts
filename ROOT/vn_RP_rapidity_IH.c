#include <iostream>
#include <fstream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>

using namespace std;

//------ main block ----------
void vn_RP_rapidity_IH(const char* direct, int Npart_min, int Npart_max, double order, int pid)
{
  cout << endl << endl << "Calculating v_" << (int)order << "(rapidity) of identified hadrons (" << pid << ")..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Cuts
  const double etaCut = 1.0;
  const double ptMinCut = 0.2;
  const double ptMaxCut = 1.6;
  const double yMin = -1.0;
  const double yMax =  1.0;
  const int nBins = 10;
  const double dy = (yMax-yMin)/nBins;

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
  Short_t ele[NP];
  Int_t npart, Nparticipants;

  int nevents = tree->GetEntries();
  cout << "Total number of events in directory: " << nevents << endl;
  int centrality_events[nBins] = {0};
  tree->SetBranchAddress("px",&px[0]);
  tree->SetBranchAddress("py",&py[0]);
  tree->SetBranchAddress("pz",&pz[0]);
  tree->SetBranchAddress("E",&E[0]);
  tree->SetBranchAddress("id",&id[0]);
  tree->SetBranchAddress("ele",&ele[0]);
  tree->SetBranchAddress("npart",&npart);
  tree->SetBranchAddress("Nparticipants",&Nparticipants);

  double vn[nBins]={0.0};
  double sd1[nBins]={0.0}, vnerr[nBins]={0.0};

  for(int iev=0; iev<tree->GetEntries(); iev++)
  {
    if ((iev)%((int)nevents/20) == 0) cout << (int)100*iev/nevents << "%" << endl;
    double Qx[nBins] = {0.0};
    int Np[nBins] = {0};

    tree->GetEntry(iev);
    if (Nparticipants > Npart_min && Nparticipants < Npart_max)
    {
      for(int i=0; i<npart; i++)
      {
        const float pabs = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
        const float pt = sqrt(px[i]*px[i]+py[i]*py[i]);
        const float phi = atan2(py[i],px[i]);
        const float eta = 0.5*log((pabs+pz[i])/(pabs-pz[i]));
        const float rap = 0.5*log((E[i]+pz[i])/(E[i]-pz[i]));

        if(fabs(eta)<etaCut && pt>ptMinCut && pt<ptMaxCut && id[i]==pid && rap<yMax && rap>yMin)
        {
          int yBin = (rap-yMin)/dy;
          Np[yBin]++;
          Qx[yBin] += cos(order*phi);
        }
      }

      for(int i=0; i<nBins; i++)
      {
        if (Np[i] > 0)
        {
          centrality_events[i]++;
          Qx[i] /= Np[i];
          vn[i] += Qx[i];
          sd1[i] += Qx[i] * Qx[i];
        }
      }
    }
  }

  cout << "Number of events in chosen centrality interval: " << centrality_events[(int)(nBins/2)] << endl;

  for(int i=0; i<nBins; i++)
  {
    vn[i] /= centrality_events[i];
    vnerr[i] = sqrt(sd1[i]/centrality_events[i] - vn[i]*vn[i]);
  }

  double *rapBin = new double [nBins];

  // Write results into the text file (append)
  ofstream fout;
  fout.open("vn_RP_rapidity_IH.dat", ofstream::app);
  fout << direct << "\t" << (int)order << "\t" << pid << endl;
  for(int irap = 0; irap < nBins; irap++)
  {
    vnerr[irap] = vnerr[irap] / sqrt(centrality_events[irap]);
    rapBin[irap] = yMin + (irap+0.5)*dy;
    cout << rapBin[irap] << "\t" << vn[irap] << "\t" << vnerr[irap] << endl;
    fout << rapBin[irap] << "\t" << vn[irap] << "\t" << vnerr[irap] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'vn_RP_rapidity_IH.dat'" << endl;

  delete tree;
}
