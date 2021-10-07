#include <iostream>
#include <fstream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TMath.h>

using namespace std;

//------ main block ----------
void vn_EP_mTm0_strange(const char* direct, int Npart_min, int Npart_max, double order)
{
  cout << endl << endl << "Calculating v_" << (int)order << "(mT-m0)..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Cuts
  const double etaCut = 1.0;
  const double ptMin = 0.;
  const double ptMax = 3.0;
  const int nBins = 15;
  const double dpt = (ptMax-ptMin)/nBins;

  // List of particles that will be calculated
  const int NoP = 8;
  int pids[NoP] = {
  	-2212, //pbar
  	-211, //pi-
  	-3122, //lambdabar
  	-321, //K-
  	333, //phi
  	310, //K^0_S
  	3312, //xi-
  	3334 //omega-
  };

  // Adding root files from directory into a tree
  TString dirname = direct;
  TChain *tree = new TChain("treeini");
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
  double vn_obs[nBins][NoP]={0.0}, Rn = 0.0;
  double sd1[nBins][NoP]={0.0}, vnerr[nBins][NoP]={0.0};
  int NoE[nBins][NoP]={0}; // number of non-empty events in each bin 

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
        float pabs = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
        float pt = sqrt(px[i]*px[i]+py[i]*py[i]);
        float phi = atan2(py[i],px[i]);
        float eta = 0.5*log((pabs+pz[i])/(pabs-pz[i]));
        //float mtm0 = sqrt(E[i]*E[i] - pz[i]*pz[i]) - sqrt(E[i]*E[i] - pabs*pabs);

        if(fabs(eta)<etaCut && pt>ptMin && pt<ptMax)
        {
          Qx += pt*cos(order*phi);
          Qy += pt*sin(order*phi);
        }
      }

      Etotal = 0.0;
      double _vn_obs[nBins][NoP] = {0.0};
      double QxA=0.0, QxB=0.0, QyA=0.0, QyB=0.0;
      int _nvn[nBins][NoP] = {0};

      for(int i=0; i<npart; i++)
      {
        float pabs = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
        float pt = sqrt(px[i]*px[i]+py[i]*py[i]);
        float phi = atan2(py[i],px[i]);
        float eta = 0.5*log((pabs+pz[i])/(pabs-pz[i]));
        //float mtm0 = sqrt(E[i]*E[i] - pz[i]*pz[i]) - sqrt(E[i]*E[i] - pabs*pabs);

        if(pt>ptMin && pt<ptMax)
        {
          for (int j = 0; j < NoP; j++)
          {
            if (id[i] == pids[j])
            {
              const double _Qx = Qx - pt*cos(order*phi);  // corrected flow vector
              const double _Qy = Qy - pt*sin(order*phi);  // corrected flow vector
              const double cosn = cos(order*phi);
              const double sinn = sin(order*phi);
              const double psin = atan2(_Qy,_Qx)/order;
              int ptBin = (pt-ptMin)/dpt;
              if(ptBin>=0 && ptBin<nBins)
              {
                _vn_obs[ptBin][j] += (cosn*cos(order*psin) + sinn*sin(order*psin));
                _nvn[ptBin][j]++;
              }
            }
          }
        }
        if(fabs(eta)<etaCut && pt>ptMin && pt<ptMax)
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
        for (int j = 0; j < NoP; j++)
        {
          if(_nvn[i][j]>0) 
          {
            NoE[i][j]++;
            _vn_obs[i][j] /= _nvn[i][j];
            vn_obs[i][j] += _vn_obs[i][j];
            sd1[i][j] += _vn_obs[i][j] * _vn_obs[i][j];
          }
        }
      }

      Etotaver += Etotal;
      Rn += cos(order*(psi2A-psi2B));
    }
  }

  nevents = centrality_events;
  cout << "Number of events in chosen centrality interval: " << centrality_events << endl;

  for(int i=0; i<nBins; i++)
  {
    for (int j = 0; j < NoP; j++)
    {
      vn_obs[i][j] /= NoE[i][j];
      vnerr[i][j] = sqrt(sd1[i][j]/NoE[i][j] - vn_obs[i][j]*vn_obs[i][j]);
    }
  }

  Etotaver /= nevents;

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
  double *ptBin = new double [nBins];
  //double *vn = new double [nBins];

  // Write results into the text file (append)
  ofstream fout;
  fout.open("vn_EP_mTm0.dat", ofstream::app);
  fout << direct << "\t" << (int)order << endl;
  for(int i = 0; i < nBins; i++)
  {
    ptBin[i] = ptMin + (i+0.5)*dpt;
    cout << ptBin[i] << "\t";
    fout << ptBin[i] << "\t";
    for (int j = 0; j < NoP; j++)
    {
      vn_obs[i][j] = vn_obs[i][j]/Rn;
      vnerr[i][j] = vnerr[i][j] / Rn / sqrt(NoE[i][j]);
      cout << vn_obs[i][j] << "\t" << vnerr[i][j] << "\t";
      fout << vn_obs[i][j] << "\t" << vnerr[i][j] << "\t";
    }
    cout << endl;
    fout << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'vn_EP_mTm0.dat'" << endl;

  delete tree;
}
