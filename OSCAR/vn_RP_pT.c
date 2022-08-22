#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// NoF - number of files
// NoE - number of events

void vn_RP_pT(const char* direct, int NoF, int NoE, double order){

    cout << "Processing events from directory: " << direct << endl;

    // Cuts
    const double etaCut = 0.8 ;
    const double ptMinCut = 0.2 ;
    const double ptMaxCut = 3.0 ;
    const int nBins = 14 ;
    const double dpt = (ptMaxCut-ptMinCut)/nBins ;

    double vn[nBins]={0.0};
    double sd1[nBins]={0.0}, vnerr[nBins]={0.0};

    int nevents[nBins] = {0};

    int npar, maxpar;
    Int_t npart, Nparticipants ;

    double Rn = 0.0;

    double E, px, py, pz;
    int ele, id;


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

        if (infile == NULL){

            cout << "Warning: Missing file #" << ifls << endl;
        }
        else{

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
                    double Qx[nBins] = {0.0};
                    int Np[nBins] = {0};

                    // First loop over particles - calculate flow vector
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

                        id = atof(pars[9]);
                        E = atof(pars[5]);
                        px = atof(pars[6]);
                        py = atof(pars[7]);
                        pz = atof(pars[8]);
                        ele = atoi(pars[11]);

                        double pabs = sqrt(px*px+py*py+pz*pz);
                        double pt = sqrt(px*px+py*py);
                        double eta = 0.5*log((pabs+pz)/(pabs-pz));
                        double phi = atan2(py,px);
                        double rap =  0.5*log((E+pz)/(E-pz));

                        if (abs(eta)<etaCut && pt>ptMinCut && pt<ptMaxCut)
                        {
                            int ptBin = (pt - ptMinCut) / dpt;
                            Np[ptBin]++;
                            Qx[ptBin] += cos(order*phi);
                        }
                    }

                    for (int i=0; i<nBins; i++)
                    {
                        if (Np[i] > 0)
                        {
                            nevents[i]++;
                            Qx[i] /= Np[i];
                            vn[i] += Qx[i];
                            sd1[i] += Qx[i] * Qx[i];
                        }
                    }
                }
                fgets(line,500,infile);
            }
        }
    }

    cout << "Number of events: " << nevents[0] << endl;

    for (int i=0; i<nBins; i++)
    {
        vn[i] /= nevents[i];
        vnerr[i] = sqrt(sd1[i]/nevents[i] - vn[i]*vn[i]);
    }

    // Write results into the text file (append)
    ofstream fout;
    fout.open("vn_RP_pT.dat", ofstream::app);
    fout << direct << "\t" << (int)order << endl;

    double ptBin[nBins] = {0.0};

    for (int ipt=0; ipt<nBins; ipt++)
    {
        vnerr[ipt] = vnerr[ipt] / sqrt(nevents[ipt]);
        ptBin[ipt] = ptMinCut + (ipt + 0.5) * dpt;
        cout << ptBin[ipt] << "\t" << vn[ipt] << "\t" << vnerr[ipt] << endl;
        fout << ptBin[ipt] << "\t" << vn[ipt] << "\t" << vnerr[ipt] << endl;
    }

    fout << endl;
    fout.close();

    cout << "Results have been written to 'vn_RP_pT.dat'" << endl;
}
