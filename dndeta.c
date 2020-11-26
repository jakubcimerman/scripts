#include <iostream>
#include <iomanip>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
const int SOL = 1 ;
const int DASH = 7 ;
const int DOT = 3 ;
TLegend *leg ;
void plotEmission(const char* filename, int Npart_min, int Npart_max);

void configureSizes(TH1 *g)
{
 g->SetTitle("");
 g->GetXaxis()->SetTitle("#eta");
 g->GetXaxis()->CenterTitle();
 g->GetXaxis()->SetTitleOffset(1.1);
 g->GetXaxis()->SetTitleSize(0.045);
 g->GetXaxis()->SetLabelSize(0.045);
 g->GetYaxis()->SetTitle("dN_{ch}/d#eta");
 g->GetYaxis()->CenterTitle();
 g->GetYaxis()->SetTitleOffset(1.1);
 g->GetYaxis()->SetTitleSize(0.045);
 g->GetYaxis()->SetLabelSize(0.045);
}

void dndeta()
{
// ---- visc ----
  plotEmission("../hybrid/escape5.out/bigrun-27-550-glis-R1/", 232, 321) ;
  plotEmission("../hybrid/escape5.out/bigrun-27-550-glis-R1/", 165, 232) ;
  plotEmission("../hybrid/escape5.out/bigrun-27-550-glis-R1/", 114, 165) ;
  plotEmission("../hybrid/escape5.out/bigrun-27-550-glis-R1/", 76, 114) ;
 // plotEmission("../hybrid/escape5.out/bigrun-27-550-glis-R1/", 60, 93) ;
  plotEmission("../hybrid/escape5.out/bigrun-27-550-trento-R1/", 232, 321) ;
  plotEmission("../hybrid/escape5.out/bigrun-27-550-trento-R1/", 165, 232) ;
  plotEmission("../hybrid/escape5.out/bigrun-27-550-trento-R1/", 114, 165) ;
  plotEmission("../hybrid/escape5.out/bigrun-27-550-trento-R1/", 76, 114) ;
 // plotEmission("../hybrid/escape5.out/bigrun-27-550-trento-R1/", 60, 93) ;
}


void plotEmission(const char* direct, int Npart_min, int Npart_max)
// autocorrelations removed
{
 //const double etaCut = 7.0 ;
 //const double ptMinCut = 0.2 ;
 //const double ptMaxCut = 2.0 ;
 const double etaMin = -5.0 ;
 const double etaMax =  5.0 ;
 const int nBins = 50 ;
 const int eventStep = 1;
 const double deta = (etaMax-etaMin)/nBins ;
 TString dirname = direct;
 TChain *tree = new TChain("treefin");
 TSystemDirectory dir ("rootfiles",dirname);
 TList *files = dir.GetListOfFiles();
 TSystemFile *file;
 TIter next(files);
 while ((file=(TSystemFile*)next())) {
  TString fname = file->GetName();
  cout << "[" << fname.Data() << "]";
  if(strstr(fname,".root"))
   tree->Add(dirname+fname);
 }
 //TChain* tree = new TChain("treefin;1") ;
 //tree->Add(filename) ;
 static int iplot=0 ;
 const int NP = 60000 ;
 Int_t id[NP] ;
 Float_t px [NP], py[NP], pz[NP], E[NP] ;
 Short_t ele[NP] ;
 Int_t npart, Nparticipants ;
 int nevents = tree->GetEntries() ;
 int centrality_events = 0;
 tree->SetBranchAddress("px",&px[0]) ;
 tree->SetBranchAddress("py",&py[0]) ;
 tree->SetBranchAddress("pz",&pz[0]) ;
 tree->SetBranchAddress("E",&E[0]) ;
 tree->SetBranchAddress("id",&id[0]) ;
 tree->SetBranchAddress("ele",&ele[0]) ;
 tree->SetBranchAddress("npart",&npart) ;
 tree->SetBranchAddress("Nparticipants",&Nparticipants) ;
 char hname [255] ;
 sprintf(hname,"hv2ch_%i",iplot) ;

 double dN[nBins]={0.0};
 for(int iev=0; iev<nevents; iev++){
    // loop #1, total flow vector
    tree->GetEntry(iev) ;
    if (Nparticipants > Npart_min && Nparticipants < Npart_max) {  
      centrality_events++;
      for(int i=0; i<npart; i++){
  //  if(!(E[i]<1000. && E[i]>0.)) continue ;
      const float eta = 0.5*log((sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])+pz[i])/(sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-pz[i])) ;
      const float pt = sqrt(px[i]*px[i]+py[i]*py[i]) ;
   
      if (ele[i]!=0 && eta < etaMax && eta > etaMin) {
        int etaBin = (eta-etaMin)/deta ;
        dN[etaBin]++;
      }
    } // end particle loop #1
  }
}

cout << centrality_events << endl;

ofstream fout;
 fout.open("dndeta.dat", ofstream::app);
fout << direct << endl;
for (int i = 0; i < nBins; i++) {
  dN[i] /= (centrality_events*deta);
  fout << etaMin + (i+0.5)*deta  << "\t" << dN[i] << endl;
}
fout << endl;
fout.close();

}

