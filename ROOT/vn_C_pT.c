#include <iostream>
#include <fstream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>

using namespace std;

//------ main block ----------
void vn_C_pT(const char* direct, int Npart_min, int Npart_max, double order)
{
  cout << endl << endl << "Calculating v_" << (int)order << "{2} (pT)..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Cuts
  const double etaCut = 1.0;
  const double ptMinCut = 0.2;
  const double ptMaxCut = 3.0;
  const int ptBins = 14;
  const int eventStep = 50; // number of events in one super-event (to increase statistics)
  const int NS = 10; // number of subsamples to estimate the error
  const double dpt = (ptMaxCut-ptMinCut)/ptBins;

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
  tree->SetBranchAddress("px",&px[0]);
  tree->SetBranchAddress("py",&py[0]);
  tree->SetBranchAddress("pz",&pz[0]);
  tree->SetBranchAddress("E",&E[0]);
  tree->SetBranchAddress("id",&id[0]);
  tree->SetBranchAddress("ele",&ele[0]);
  tree->SetBranchAddress("npart",&npart);
  tree->SetBranchAddress("Nparticipants",&Nparticipants);

  // arrays for calculation v_n
  double vn2[ptBins] = {0.0}, vn4[ptBins] = {0.0};
  double dn2[ptBins] = {0.0}, dn2_nom[ptBins] = {0.0}, dn2_denom[ptBins] = {0.0};
  double dn4[ptBins] = {0.0}, dn4_nom[ptBins] = {0.0}, dn4_denom[ptBins] = {0.0};
  double cn2 = 0.0, cn2_nom = 0.0, cn2_denom = 0.0;
  double cn4 = 0.0, cn4_nom = 0.0, cn4_denom = 0.0;
  double vn2_err[ptBins] = {0.0}, vn4_err[ptBins] = {0.0};

  // arrays for calculation error of v_n
  double cn2_sub[NS] = {0.0}, cn2_nom_sub[NS] = {0.0}, cn2_denom_sub[NS] = {0.0};
  double cn4_sub[NS] = {0.0}, cn4_nom_sub[NS] = {0.0}, cn4_denom_sub[NS] = {0.0};
  double ** dn2_sub = new double*[NS];
  double ** dn2_nom_sub = new double*[NS];
  double ** dn2_denom_sub = new double*[NS];
  double ** dn4_sub = new double*[NS];
  double ** dn4_nom_sub = new double*[NS];
  double ** dn4_denom_sub = new double*[NS];
  double ** vn2_sub = new double*[NS];
  double ** vn4_sub = new double*[NS];
  for (int i = 0; i < NS; i++) {
    dn2_sub[i] = new double[ptBins];
    dn2_nom_sub[i] = new double[ptBins];
    dn2_denom_sub[i] = new double[ptBins];
    dn4_sub[i] = new double[ptBins];
    dn4_nom_sub[i] = new double[ptBins];
    dn4_denom_sub[i] = new double[ptBins];
    vn2_sub[i] = new double[ptBins];
    vn4_sub[i] = new double[ptBins];
    for (int j = 0; j < ptBins; j++) {
      dn2_sub[i][j] = 0.0;
      dn2_nom_sub[i][j] = 0.0;
      dn2_denom_sub[i][j] = 0.0;
      dn4_sub[i][j] = 0.0;
      dn4_nom_sub[i][j] = 0.0;
      dn4_denom_sub[i][j] = 0.0;
      vn2_sub[i][j] = 0.0;
      vn4_sub[i][j] = 0.0;
    }
  }

  // loop through events
  for(int iev=0; iev<tree->GetEntries(); iev+=eventStep)
  {
    int RFP = 0, POI[ptBins] = {0}; // RFP = M (all particles), POI = m (particles in given ptBin)
    double cumulant2 = 0.0, diff_cumulant2[ptBins] = {0.0};
    double cumulant4 = 0.0, diff_cumulant4[ptBins] = {0.0};
    double Qx = 0.0, Qy = 0.0;
    double Qx2 = 0.0, Qy2 = 0.0;
    double qx[ptBins] = {0.0}, qy[ptBins] = {0.0};
    double qx2[ptBins] = {0.0}, qy2[ptBins] = {0.0};
    complex<double> Q, Q2, q[ptBins], q2[ptBins];


    // loop through super-event
    for(int k = 0; k < eventStep; k++)
    {
      if ((iev+k)%((int)nevents/20) == 0) cout << (int)100*(iev+k)/nevents << "%" << endl;
      tree->GetEntry(iev+k); 
      // centrality cut 
      if (Nparticipants > Npart_min && Nparticipants < Npart_max)
      {
        // loop through particles
        for(int i=0; i<npart; i++)
        {
          double pabs = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
          double pt = sqrt(px[i]*px[i]+py[i]*py[i]);
          double eta = 0.5*log((pabs+pz[i])/(pabs-pz[i]));
          double phi = atan2(py[i], px[i]);

          if(fabs(eta)<etaCut && pt>ptMinCut && pt<ptMaxCut)
          {
            RFP++; 
            int ptBin = (pt-ptMinCut)/dpt;
            POI[ptBin]++;
            Qx += cos(order*phi);
            Qy += sin(order*phi);
            Qx2 += cos(2*order*phi);
            Qy2 += sin(2*order*phi);
            qx[ptBin] += cos(order*phi);
            qy[ptBin] += sin(order*phi);
            qx2[ptBin] += cos(2*order*phi);
            qy2[ptBin] += sin(2*order*phi);
          }
        }
      }
    }

    if (RFP > 0)
    {
      Q = complex<double>(Qx, Qy);
      Q2 = complex<double>(Qx2, Qy2);
      // calculation of cumulants
      cumulant2 = (double)(abs(Q)*abs(Q)-RFP)/(RFP*(RFP-1));
      cumulant4 = (pow(abs(Q),4) + pow(abs(Q2),2) - 2*real(Q2*conj(Q)*conj(Q)) - 4*(RFP-2)*abs(Q)*abs(Q) 
        + 2*RFP*(RFP-3)) / (RFP*(RFP-1)*(RFP-2)*(RFP-3));

      cn2_nom += RFP*(RFP-1)*cumulant2;
      cn2_denom += RFP*(RFP-1);
      cn4_nom += (RFP*(RFP-1)*(RFP-2)*(RFP-3))*cumulant4;
      cn4_denom += RFP*(RFP-1)*(RFP-2)*(RFP-3);
      
      // distribute events to NS samples for error calculation
      int r = rand() % NS;
      cn2_nom_sub[r] += RFP*(RFP-1)*cumulant2;
      cn2_denom_sub[r] += RFP*(RFP-1);
      cn4_nom_sub[r] += (RFP*(RFP-1)*(RFP-2)*(RFP-3))*cumulant4;
      cn4_denom_sub[r] += RFP*(RFP-1)*(RFP-2)*(RFP-3);

      for (int ipt = 0; ipt < ptBins; ipt++)
      {
        if (POI[ipt] > 0)
        {
          q[ipt] = complex<double>(qx[ipt], qy[ipt]);
          q2[ipt] = complex<double>(qx2[ipt], qy2[ipt]);
          diff_cumulant2[ipt] = real(q[ipt]*conj(Q)-complex<double>(POI[ipt],0))/(POI[ipt]*RFP-POI[ipt]);
          complex<double> term1 = q[ipt]*Q*conj(Q)*conj(Q);
          complex<double> term2 = q2[ipt]*conj(Q)*conj(Q);
          complex<double> term3 = q[ipt]*Q*conj(Q2);
          complex<double> term4 = complex<double>((9-2*RFP),0)*q[ipt]*conj(Q);
          complex<double> term5 = 2*POI[ipt]*abs(Q)*abs(Q);
          complex<double> term6 = Q*conj(q[ipt]);
          complex<double> term7 = q2[ipt]*conj(Q2);
          diff_cumulant4[ipt] = real(term1 - term2 - term3 + term4 - term5 - term6 + term7
          + complex<double>(2*POI[ipt]*RFP - 6*POI[ipt],0)) / ((POI[ipt]*RFP-3*POI[ipt])*(RFP-1)*(RFP-2));

          dn2_nom[ipt] += (POI[ipt]*RFP - POI[ipt]) * diff_cumulant2[ipt];
          dn2_denom[ipt] += POI[ipt]*RFP - POI[ipt];

          dn4_nom[ipt] += (POI[ipt]*RFP-3*POI[ipt])*(RFP-1)*(RFP-2) * diff_cumulant4[ipt];
          dn4_denom[ipt] += (POI[ipt]*RFP-3*POI[ipt])*(RFP-1)*(RFP-2);

          dn2_nom_sub[r][ipt] += (POI[ipt]*RFP - POI[ipt]) * diff_cumulant2[ipt];
          dn2_denom_sub[r][ipt] += POI[ipt]*RFP - POI[ipt];

          dn4_nom_sub[r][ipt] += (POI[ipt]*RFP-3*POI[ipt])*(RFP-1)*(RFP-2) * diff_cumulant4[ipt];
          dn4_denom_sub[r][ipt] += (POI[ipt]*RFP-3*POI[ipt])*(RFP-1)*(RFP-2);
        }
      }
    }
  }

  // final calculation of v_n
  cn2 = cn2_nom / cn2_denom;
  cn4 = cn4_nom / cn4_denom - 2 * cn2 * cn2;

  for (int j = 0; j < NS; j++)
  {
    cn2_sub[j] = cn2_nom_sub[j] / cn2_denom_sub[j];
    cn4_sub[j] = cn4_nom_sub[j] / cn4_denom_sub[j] - 2 * cn2_sub[j] * cn2_sub[j];
  }

  for (int ipt = 0; ipt < ptBins; ipt++)
  {
    dn2[ipt] = dn2_nom[ipt] / dn2_denom[ipt];
    vn2[ipt] = dn2[ipt] / sqrt(cn2);

    dn4[ipt] = dn4_nom[ipt] / dn4_denom[ipt] - 2 * dn2[ipt] * cn2;
    vn4[ipt] = -dn4[ipt] / pow(-cn4, 0.75);

    double vn2_mean = 0.0, vn2_sd = 0.0, vn4_mean = 0.0, vn4_sd = 0.0;
    for (int j = 0; j < NS; j++)
    {
      dn2_sub[j][ipt] = dn2_nom_sub[j][ipt] / dn2_denom_sub[j][ipt];
      vn2_sub[j][ipt] = dn2_sub[j][ipt] / sqrt(cn2_sub[j]);

      dn4_sub[j][ipt] = dn4_nom_sub[j][ipt] / dn4_denom_sub[j][ipt] - 2 * dn2_sub[j][ipt] * cn2_sub[j];
      vn4_sub[j][ipt] = -dn4_sub[j][ipt] / pow(-cn4_sub[j], 0.75);

      vn2_mean += vn2_sub[j][ipt];
      vn2_sd += vn2_sub[j][ipt] * vn2_sub[j][ipt];
      vn4_mean += vn4_sub[j][ipt];
      vn4_sd += vn4_sub[j][ipt] * vn4_sub[j][ipt];
    }
    vn2_mean /= NS;
    vn4_mean /= NS;

    vn2_err[ipt] = sqrt(vn2_sd / NS - vn2_mean * vn2_mean) / sqrt(NS);
    vn4_err[ipt] = sqrt(vn4_sd / NS - vn4_mean * vn4_mean) / sqrt(NS);
  }
  

  // Write results into the text file (append)
  ofstream fout;
  fout.open("vn_cumulants_pT.dat", ofstream::app);
  fout << direct << "\t" << (int)order << endl;
  for(int ipt = 0; ipt < ptBins; ipt++)
  {
    double ptBin = (ipt+0.5)*dpt + ptMinCut;
    cout << ptBin << "\t" << vn2[ipt] << "\t" << vn2_err[ipt] << "\t" << vn4[ipt] << "\t" << vn4_err[ipt] << endl;
    fout << ptBin << "\t" << vn2[ipt] << "\t" << vn2_err[ipt] << "\t" << vn4[ipt] << "\t" << vn4_err[ipt] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'vn_cumulants_pT.dat'" << endl;

  for (int i = 0; i < NS; i++)
  {
    delete[] dn2_sub[i];
    delete[] dn2_nom_sub[i];
    delete[] dn2_denom_sub[i];
    delete[] dn4_sub[i];
    delete[] dn4_nom_sub[i];
    delete[] dn4_denom_sub[i];
    delete[] vn2_sub[i];
    delete[] vn4_sub[i];
  }
  delete[] dn2_sub;
  delete[] dn2_nom_sub;
  delete[] dn2_denom_sub;
  delete[] dn4_sub;
  delete[] dn4_nom_sub;
  delete[] dn4_denom_sub;
  delete[] vn2_sub;
  delete[] vn4_sub;

  delete tree;
}
