#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void spectrum(const char* direct, int NoF, int NoE)
{
  cout << endl << endl << "Calculating dN/dp_t..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Borders of pt
  const double ptMin = 0.2;
  const double ptMax = 2.0;
  const double yCut = 0.1;
  const int nBins = 36;
  const double pi = 3.1415926535897932;
  const double dpt = (ptMax-ptMin)/nBins;

  // Maximum number of particles - length of arrays
  //const int NP = 60000;

  double px, py, pz, m, E;
  double ele;
  int id;
  int npart, nevents = 0;
  int npar, maxpar;

  double dNpi[nBins]={0.0}, dNK[nBins]={0.0}, dNp[nBins]={0.0};
  double dNpi_ev[nBins]={0.0}, dNK_ev[nBins]={0.0}, dNp_ev[nBins]={0.0};
  double dNpi_err[nBins]={0.0}, dNK_err[nBins]={0.0}, dNp_err[nBins]={0.0};

  // Loop over files
  for (int ifls = 1; ifls < NoF + 1; ifls++) {
    FILE *infile;
    char line[500];
    char delims[] = " ,\n\t";
    char *strtokresult = NULL;
    char *pars[12];

    string filename;
    string buffer;
    filename.append(direct);
    filename.append(to_string(ifls));
    filename.append("_fin.oscar");

    infile = fopen(filename.c_str(), "r");
    if (infile == NULL)
    {
      cout << "Warning: Missing file #" << ifls << endl;
    }
    else
    {
      fgets(line,500,infile);
      fgets(line,500,infile);
      fgets(line,500,infile);
      // Loop over events in file
      for (int iev = 0; iev < NoE; iev++)
      {
        fgets(line,500,infile);
        strtokresult = strtok(line, delims);
        npar = 0;
        maxpar = 5;
        while( (strtokresult != NULL) && (npar < maxpar) )
        {
          pars[npar]= strtokresult;
          strtokresult = strtok( NULL, delims );
          npar += 1;
        }
        int npart = atoi(pars[4]);
        if (npart > 0)
        {
          nevents++;
          for (int ibin = 0; ibin < nBins; ibin++)
          {
            dNpi_ev[ibin] = 0.0;
            dNK_ev[ibin] = 0.0;
            dNp_ev[ibin] = 0.0;
          }
          // Loop over particles
          for (int i = 0; i < npart; i++)
          {
            fgets(line,500,infile);
            strtokresult = strtok(line, delims);
            npar = 0;
            maxpar = 12;
            while( (strtokresult != NULL) && (npar < maxpar) )
            {
              pars[npar]= strtokresult;
              strtokresult = strtok( NULL, delims );
              npar += 1;
            }
            m  = atof(pars[4]);
            E  = atof(pars[5]);
            px = atof(pars[6]);
            py = atof(pars[7]);
            pz = atof(pars[8]);
            id = atoi(pars[9]);

            float pt = sqrt(px*px+py*py);
            //float p = sqrt(px*px+py*py+pz*pz);
            //float E = sqrt(m*m+p*p);
            float rap = 0.5*log((E+pz)/(E-pz));
            if (abs(rap) < yCut && pt > ptMin && pt < ptMax) {
              int ptBin = (pt-ptMin)/dpt ;
              if (id == 211) dNpi_ev[ptBin]++;
              if (id == 321) dNK_ev[ptBin]++;
              if (id == 2212) dNp_ev[ptBin]++;
            }

          }

          for (int ibin = 0; ibin < nBins; ibin++)
          {
            dNpi[ibin] += dNpi_ev[ibin];
            dNK[ibin] += dNK_ev[ibin];
            dNp[ibin] += dNp_ev[ibin];
            dNpi_err[ibin] += dNpi_ev[ibin] * dNpi_ev[ibin];
            dNK_err[ibin] += dNK_ev[ibin] * dNK_ev[ibin];
            dNp_err[ibin] += dNp_ev[ibin] * dNp_ev[ibin];
          }
          fgets(line,500,infile);
        }
      }
    }
    if ((ifls)%((int)NoF/20) == 0) cout << (int)100*ifls/NoF << "%" << endl;
  }

  cout << "Total number of events: " << nevents << endl;

  // Write results into the text file (append)
  ofstream fout;
  fout.open("spectrum.dat", ofstream::app);
  fout << direct << endl;
  for (int i = 0; i < nBins; i++)
  {
    dNpi[i] /= nevents;
    dNK[i] /= nevents;
    dNp[i] /= nevents;
    dNpi_err[i] = sqrt(dNpi_err[i] / nevents - dNpi[i] * dNpi[i]) / sqrt(nevents);
    dNK_err[i] = sqrt(dNK_err[i] / nevents - dNK[i] * dNK[i]) / sqrt(nevents);
    dNp_err[i] = sqrt(dNp_err[i] / nevents - dNp[i] * dNp[i]) / sqrt(nevents);
    dNpi[i] /= (2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNK[i] /= (2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNp[i] /= (2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNpi_err[i] /= (2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNK_err[i] /= (2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNp_err[i] /= (2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    fout << ptMin + (i+0.5)*dpt  << "\t" << dNpi[i] << "\t" << dNpi_err[i] << "\t"
         << dNK[i] << "\t" << dNK_err[i] << "\t" << dNp[i] << "\t" << dNp_err[i] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'spectrum.dat'" << endl;
}
