#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void dndeta(const char* direct, int NoF, int NoE)
{
  cout << endl << endl << "Calculating dN/deta..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Borders of eta
  const double etaMin = -5.0;
  const double etaMax =  5.0;
  const int nBins = 50;
  const double deta = (etaMax-etaMin)/nBins;
  double dN[nBins] = {0.0};

  // Maximum number of particles - length of arrays
  //const int NP = 60000;

  double px, py, pz;
  double ele;
  int npart, nevents = 0;

  int npar, maxpar;

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
            px = atof(pars[6]);
            py = atof(pars[7]);
            pz = atof(pars[8]);
            ele = atoi(pars[11]);
            float eta = 0.5*log((sqrt(px*px+py*py+pz*pz)+pz)/(sqrt(px*px+py*py+pz*pz)-pz));
            if (ele != 0 && eta < etaMax && eta > etaMin)
            {
              int etaBin = (eta-etaMin)/deta;
              dN[etaBin]++;
            }
          }

          fgets(line,500,infile);
        }
      }
    }
    if ((ifls)%((int)NoF/20) == 0) cout << (int)100*ifls/NoF << "%" << endl;
  }

  cout << "Number of events in chosen centrality interval: " << nevents << endl;

  // Write results into the text file (append)
  ofstream fout;
  fout.open("dndeta.dat", ofstream::app);
  fout << direct << endl;
  for (int i = 0; i < nBins; i++)
  {
    dN[i] /= (nevents*deta);
    fout << etaMin + (i+0.5)*deta  << "\t" << dN[i] << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to 'dndeta.dat'" << endl;
}
