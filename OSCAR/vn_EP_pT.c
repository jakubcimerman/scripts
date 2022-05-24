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

void vn_EP_pT(const char* direct, int NoF, int NoE, double order){

    cout << "Processing events from directory: " << direct << endl;

    // Cuts
    const double etaCut = 1.0 ;
    const double ptMinCut = 0.2 ;
    const double ptMaxCut = 3.0 ;
    const int nBins = 14 ;
    const double dpt = (ptMaxCut-ptMinCut)/nBins ;

    double vn_obs[nBins]={0.0};
    double sd1[nBins]={0.0}, vnerr[nBins]={0.0};

    int nevents[nBins] = {0};

    int npar, maxpar;
    Int_t npart, Nparticipants ;
    int centrality_events = 0;

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
                    double Qx = 0.0, Qy = 0.0;
                    vector<double> E_ev, px_ev, py_ev, pz_ev;
                    vector<int> ele_ev, id_ev;

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

                        if (abs(eta)<etaCut && pt>ptMinCut && pt<ptMaxCut)
                        {
                            Qx += pt*cos(order*phi);
                            Qy += pt*sin(order*phi);
                        }

                        id_ev.push_back(id);
                        E_ev.push_back(E);
                        px_ev.push_back(px);
                        py_ev.push_back(py);
                        pz_ev.push_back(pz);
                        ele_ev.push_back(ele);
                    }

                    double _vn_obs[nBins] = {0.0};
                    double QxA = 0.0, QxB = 0.0, QyA = 0.0, QyB = 0.0;
                    int _nvn[nBins] = {0};

                    // Second loop over particles
                    for (int i = 0; i < npart; i++)
                    {
                        double pabs = sqrt(px_ev[i]*px_ev[i]+py_ev[i]*py_ev[i]+pz_ev[i]*pz_ev[i]);
                        double pt = sqrt(px_ev[i]*px_ev[i]+py_ev[i]*py_ev[i]);
                        double eta = 0.5*log((pabs+pz_ev[i])/(pabs-pz_ev[i]));
                        double phi = atan2(py_ev[i],px_ev[i]);

                        if (pt>ptMinCut && pt<ptMaxCut && abs(eta)<etaCut && ele_ev[i]!=0)
                        {
                            double _Qx = Qx - pt*cos(order*phi);  // corrected flow vector
                            double _Qy = Qy - pt*sin(order*phi);  // corrected flow vector
                            double cosn = cos(order*phi);
                            double sinn = sin(order*phi);
                            double psin = atan2(_Qy,_Qx)/order;
                            int ptBin = floor((pt-ptMinCut)/dpt);
                            _vn_obs[ptBin] += (cosn*cos(order*psin) + sinn*sin(order*psin));
                            _nvn[ptBin]++;
                        }

                        if (abs(eta)<etaCut && pt>ptMinCut && pt<ptMaxCut)
                        {
                            if (i%2 == 0)
                            {
                                QxA += pt*cos(order*phi);
                                QyA += pt*sin(order*phi);
                            }
                            else
                            {
                                QxB += pt*cos(order*phi);
                                QyB += pt*sin(order*phi);
                            }
                        }
                    }

                    double psi2A = atan2(QyA,QxA)/order;
                    double psi2B = atan2(QyB,QxB)/order;
                    Rn += cos(order*(psi2A-psi2B));

                    for (int i = 0; i < nBins; i++)
                    {
                        if (_nvn[i]>0)
                        {
                             nevents[i]++;
                             _vn_obs[i] /= _nvn[i];
                        }
                        vn_obs[i] += _vn_obs[i];
                        sd1[i] += _vn_obs[i] * _vn_obs[i];
                    }
                }
                fgets(line,500,infile);
            }
        }
    }

    cout << "Number of events: " << nevents[0] << endl;

    for (int i=0; i<nBins; i++)
    {
        vn_obs[i] /= nevents[i];
        vnerr[i] = sqrt(sd1[i]/nevents[i] - vn_obs[i]*vn_obs[i]);
    }

    Rn = sqrt(Rn/nevents[0]) ;
    cout << "Rn^sub = " << Rn << endl;

    const double pf = sqrt(M_PI)/(2.0*sqrt(2.0)) ;
    double ksiMin=0., ksiMax = 20. ;

    while(ksiMax-ksiMin>0.01)
    {
        double ksi = 0.5*(ksiMin+ksiMax) ;
        double R = pf*ksi*exp(-0.25*ksi*ksi)*(TMath::BesselI0(0.25*ksi*ksi)+TMath::BesselI1(0.25*ksi*ksi)) ;
        if(R>Rn) ksiMax = ksi ;
        else ksiMin = ksi ;
    }

    const double ksi = sqrt(2)*0.5*(ksiMin+ksiMax) ;

    cout << "ksi = " << ksi << endl;

    Rn = pf*ksi*exp(-0.25*ksi*ksi)*(TMath::BesselI0(0.25*ksi*ksi)+TMath::BesselI1(0.25*ksi*ksi));

    // Write results into the text file (append)
    ofstream fout;
    fout.open("vn_EP_pT.dat", ofstream::app);
    fout << direct << "\t" << (int)order << endl;

    double vn[nBins] = {0.0};
    double ptBin[nBins] = {0.0};

    for (int ipt=0; ipt<nBins; ipt++)
    {
        vn[ipt] = vn_obs[ipt] / Rn;
        vnerr[ipt] = vnerr[ipt] / Rn / sqrt(nevents[ipt]);
        ptBin[ipt] = ptMinCut + (ipt + 0.5) * dpt;
        cout << ptBin[ipt] << "\t" << vn[ipt] << "\t" << vnerr[ipt] << endl;
        fout << ptBin[ipt] << "\t" << vn[ipt] << "\t" << vnerr[ipt] << endl;
    }

    fout << endl;
    fout.close();

    cout << "Results have been written to 'vn_EP_pT.dat'" << endl;
}
