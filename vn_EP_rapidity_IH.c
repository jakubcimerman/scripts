#include <iostream>
#include <fstream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>

using namespace std;

//------ main block ----------
void vn_EP_rapidity_IH(const char* direct, int Npart_min, int Npart_max, double order, int pid)
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
  int centrality_events = 0;
  tree->SetBranchAddress("px",&px[0]);
  tree->SetBranchAddress("py",&py[0]);
  tree->SetBranchAddress("pz",&pz[0]);
  tree->SetBranchAddress("E",&E[0]);
  tree->SetBranchAddress("id",&id[0]);
  tree->SetBranchAddress("ele",&ele[0]);
  tree->SetBranchAddress("npart",&npart);
  tree->SetBranchAddress("Nparticipants",&Nparticipants);

  double Etotal, Etotaver = 0.0;
  double vn_obs[nBins]={0.0}, Rn = 0.0;
  double sd1[nBins]={0.0}, vnerr[nBins]={0.0};

  for(int iev=0; iev<tree->GetEntries(); iev++)
  {
    if ((iev)%((int)nevents/20) == 0) cout << (int)100*iev/nevents << "%" << endl;
    double Qx = 0.0, Qy = 0.0;
    tree->GetEntry(iev);
    if (Nparticipants > Npart_min && Nparticipants < Npart_max)
    {
      centrality_events++;
      for(int i=0; i<npart; i++)
      {
        const float pabs = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
        const float pt = sqrt(px[i]*px[i]+py[i]*py[i]);
        const float phi = atan2(py[i],px[i]);
        const float eta = 0.5*log((pabs+pz[i])/(pabs-pz[i]));

        if(fabs(eta)<etaCut && pt>ptMinCut && pt<ptMaxCut)
        {
          Qx += pt*cos(order*phi);
          Qy += pt*sin(order*phi);
        }
      }

      Etotal = 0.0;
      double _vn_obs[nBins] = {0.0};
      double QxA=0.0, QxB=0.0, QyA=0.0, QyB=0.0;
      int _nvn[nBins] = {0};

      for(int i=0; i<npart; i++)
      {
        const float pabs = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
        const float pt = sqrt(px[i]*px[i]+py[i]*py[i]);
        const float phi = atan2(py[i],px[i]);

        if(pt>ptMinCut && pt<ptMaxCut && ele[i]!=0 && id[i]==pid)
        {
          const double _Qx = Qx - pt*cos(order*phi);  // corrected flow vector
          const double _Qy = Qy - pt*sin(order*phi);  // corrected flow vector
          const double cosn = cos(order*phi);
          const double sinn = sin(order*phi);
          const double psin = atan2(_Qy,_Qx)/order;
          double rap = 0.5*log((E[i]+pz[i])/(E[i]-pz[i]));
          int yBin = floor((rap-yMin)/dy);
          if(yBin>=0 && yBin<nBins)
          {
            _vn_obs[yBin] += (cosn*cos(order*psin) + sinn*sin(order*psin));
            _nvn[yBin]++;
          }
        }
        if(fabs(0.5*log((pabs+pz[i])/(pabs-pz[i])))<etaCut && pt>ptMinCut && pt<ptMaxCut)
        {
          if(i%2==0)
          { // subevent A
            QxA += pt*cos(order*phi);
            QyA += pt*sin(order*phi);
          }
          else
          { // subevent B
            QxB += pt*cos(order*phi);
            QyB += pt*sin(order*phi);
          }
        }
        Etotal += E[i];
      }

      const double psi2A = atan2(QyA,QxA)/order;
      const double psi2B = atan2(QyB,QxB)/order;

      for(int i=0; i<nBins; i++)
      {
        if(_nvn[i]>0) _vn_obs[i] /= _nvn[i];
        vn_obs[i] += _vn_obs[i];
        sd1[i] += _vn_obs[i] * _vn_obs[i];
      }

      Etotaver += Etotal;
      Rn += cos(order*(psi2A-psi2B));
    }
  }

  nevents = centrality_events;
  cout << "Number of events in chosen centrality interval: " << centrality_events << endl;

  for(int i=0; i<nBins; i++)
  {
    vn_obs[i] /= nevents;
    vnerr[i] = sqrt(sd1[i]/nevents - vn_obs[i]*vn_obs[i]);
  }

  Etotaver /= nevents;
  cout << "vn_obs = " << vn_obs << endl;

  Rn = sqrt(Rn/nevents); // now it is R_n^sub
  cout << "R_n^sub = " << Rn << endl;
  const double pf = sqrt(TMath::Pi())/(2.0*sqrt(2.0));
  double ksiMin=0., ksiMax = 20.;
  while(ksiMax-ksiMin>0.01)
  {
    double ksi = 0.5*(ksiMin+ksiMax);
    double R = pf*ksi*exp(-0.25*ksi*ksi)*(TMath::BesselI0(0.25*ksi*ksi)+TMath::BesselI1(0.25*ksi*ksi));
    if(R>Rn) ksiMax = ksi;
    else ksiMin = ksi;
  }

  const double ksi = sqrt(2)*0.5*(ksiMin+ksiMax);
  cout << "ksi = " << ksi << endl;
  Rn = pf*ksi*exp(-0.25*ksi*ksi)*(TMath::BesselI0(0.25*ksi*ksi)+TMath::BesselI1(0.25*ksi*ksi));
  double *rapBin = new double [nBins];
  double *vn = new double [nBins];

  // Write results into the text file (append)
  ofstream fout;
  fout.open("vn_EP_rapidity_IH.dat", ofstream::app);
  fout << direct << "\t" << (int)order << "\t" << pid << endl;
  for(int irap = 0; irap < nBins; irap++)
  {
    vn[irap] = vn_obs[irap]/Rn;
    vnerr[irap] = vnerr[irap] / Rn / sqrt(nevents);
    rapBin[irap] = yMin + (irap+0.5)*dy;
    cout << rapBin[irap] << "\t" << vn[irap] << "\t" << vnerr[irap] << endl;
    fout << rapBin[irap] << "\t" << vn[irap] << "\t" << vnerr[irap] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'vn_EP_rapidity_IH.dat'" << endl;

  delete tree;
}
