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

  double px, py, pz, m, E, id;
  double ele;
  int npart, nevents = 0;
  int npar, maxpar;

  double dNpi[nBins]={0.0}, dNK[nBins]={0.0}, dNp[nBins]={0.0};

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
            id = atof(pars[9]);

            float pt = sqrt(px*px+py*py);
            //float p = sqrt(px*px+py*py+pz*pz);
            //float E = sqrt(m*m+p*p);
            float rap = 0.5*log((E+pz)/(E-pz));
            if (abs(rap) < yCut && pt > ptMin && pt < ptMax) {
              int ptBin = (pt-ptMin)/dpt ;
              if (id == 211) dNpi[ptBin]++;
              if (id == 321) dNK[ptBin]++;
              if (id == 2212) dNp[ptBin]++;
            }

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
    dNpi[i] /= (nevents*2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNK[i] /= (nevents*2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    dNp[i] /= (nevents*2*pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt));
    fout << ptMin + (i+0.5)*dpt  << "\t" << dNpi[i] << "\t" << dNK[i] << "\t" << dNp[i] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'spectrum.dat'" << endl;
}
