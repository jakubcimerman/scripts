#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// NoF - number of files
// NoE - number of events

void vn_EP_pT(const char* direct, int NoF, int NoE){

    cout << "Processing events from directory: " << direct << endl;

    // Cuts
    const double etaCut = 1.0 ;
    const double ptMinCut = 0.2 ;
    const double ptMaxCut = 3.0 ;
    const int nBins = 14 ;
    const int eventStep = 1;
    const double dpt = (ptMaxCut-ptMinCut)/nBins ;

    static int iplot=0 ;
    char hname [255] ;
    sprintf(hname,"hv2ch_%i",iplot) ;

    TH2D *havcos2ch = new TH2D(hname,hname, nBins, ptMinCut, ptMaxCut, 40, -etaCut, etaCut) ;
    float *ptBin = new float [100] ;
    float *v2 = new float [100] ;

    double vn[nBins]={0.0};
    double sd1[nBins]={0.0}, vnerr[nBins]={0.0};

    int nevents = 0;

    int npar, maxpar;
    Int_t npart, Nparticipants ;
    int centrality_events = 0;

    double Rn = 0.0;

    vector<vector<double>> E, px, py, pz;
    vector<vector<double>> pabs, pt, eta;
    vector<vector<int>> ele, id;


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
            for (int iev = 0; iev < NoE; iev++){

                vector<double> E_ev, px_ev, py_ev, pz_ev;
                vector<double> pabs_ev, pt_ev, eta_ev;
                vector<int> ele_ev, id_ev;

                fgets(line,500,infile);
                strtokresult = strtok(line, delims);
                npar = 0;
                maxpar = 5;

                while( (strtokresult != NULL) && (npar < maxpar) ){

                    pars[npar]= strtokresult;
                    strtokresult = strtok( NULL, delims );
                    npar += 1;
                }

                int npart = atoi(pars[4]);

                nevents++;

                // Loop over particles 
                for (int i = 0; i < npart; i++){

                    fgets(line,500,infile);
                    strtokresult = strtok(line, delims);
                    npar = 0;
                    maxpar = 12;

                    while( (strtokresult != NULL) && (npar < maxpar) ){

                        pars[npar]= strtokresult;
                        strtokresult = strtok( NULL, delims );
                        npar += 1;
                    }

                    int id_ = atof(pars[9]);
                    double E_ = atof(pars[5]);
                    double px_ = atof(pars[6]);
                    double py_ = atof(pars[7]);
                    double pz_ = atof(pars[8]);
                    double ele_ = atoi(pars[11]);

                    const float pabs_ = sqrt(px_*px_+py_*py_+pz_*pz_);
                    const float pt_ = sqrt(px_*px_+py_*py_);
                    const float eta_ = 0.5*log((pabs_+pz_)/(pabs_-pz_));

                    id_ev.push_back(id_);
                    pabs_ev.push_back(pabs_);
                    pt_ev.push_back(pt_);
                    eta_ev.push_back(eta_);
                    E_ev.push_back(E_);
                    px_ev.push_back(px_);
                    py_ev.push_back(py_);
                    pz_ev.push_back(pz_);
                    ele_ev.push_back(ele_);

                }

                id.push_back(id_ev);
                pabs.push_back(pabs_ev);
                pt.push_back(pt_ev);
                eta.push_back(eta_ev);
                E.push_back(E_ev);
                px.push_back(px_ev);
                py.push_back(py_ev);
                pz.push_back(pz_ev);
                ele.push_back(ele_ev);


                fgets(line,500,infile);
            }
        }
    }

    for (int k=0; k<px.size(); k++){

        double Qx = 0.0, Qy = 0.0;


        // Loop over particles #1
        for (int i = 0; i < px[k].size(); i++){


            if(fabs(0.5*log((pabs.at(k).at(i)+pz.at(k).at(i))/(pabs.at(k).at(i)-pz.at(k).at(i))))<etaCut && pt.at(k).at(i)>ptMinCut && pt.at(k).at(i)<ptMaxCut){

                Qx += (px.at(k).at(i)*px.at(k).at(i)-py.at(k).at(i)*py.at(k).at(i))/pt.at(k).at(i) ;  // which is pt*cos(2*phi)
                Qy += 2.*px.at(k).at(i)*py.at(k).at(i)/pt.at(k).at(i) ;       //  which is pt*cos(2*phi)

            }
        }

        double QxA=0.0, QxB=0.0, QyA=0.0, QyB=0.0 ;
        int _nv2 = 0 ;

        // Loop over particles #2
        for (int i = 0; i < px[k].size(); i++){


            if(pt.at(k).at(i)>ptMinCut && pt.at(k).at(i)<ptMaxCut && ele.at(k).at(i)!=0 && fabs(eta.at(k).at(i))<etaCut){

                const double _Qx = Qx - (px.at(k).at(i)*px.at(k).at(i)-py.at(k).at(i)*py.at(k).at(i))/pt.at(k).at(i) ;  //corrected flow vector
                const double _Qy = Qy - 2.*px.at(k).at(i)*py.at(k).at(i)/pt.at(k).at(i) ;       //corrected flow vector
    	        const double cos2 = (px.at(k).at(i)*px.at(k).at(i)-py.at(k).at(i)*py.at(k).at(i))/(pt.at(k).at(i)*pt.at(k).at(i)) ;
    	        const double sin2 = 2.*px.at(k).at(i)*py.at(k).at(i)/(pt.at(k).at(i)*pt.at(k).at(i)) ;
    	        const double psi2 = 0.5*atan2(_Qy,_Qx) ;
                _nv2++ ;

	            havcos2ch->Fill(sqrt(px.at(k).at(i)*px.at(k).at(i)+py.at(k).at(i)*py.at(k).at(i)), cos2*cos(2.0*psi2) + sin2*sin(2.0*psi2)) ;

            }

	        if(fabs(0.5*log((pabs.at(k).at(i)+pz.at(k).at(i))/(pabs.at(k).at(i)-pz.at(k).at(i))))<etaCut && pt.at(k).at(i)>ptMinCut && pt.at(k).at(i)<ptMaxCut){

                if(i%2==0){ //subevent A
      		        QxA += (px.at(k).at(i)*px.at(k).at(i)-py.at(k).at(i)*py.at(k).at(i))/pt.at(k).at(i) ;  // which is pt*cos(2*phi)
      		        QyA += 2.*px.at(k).at(i)*py.at(k).at(i)/pt.at(k).at(i) ;       // which is pt*cos(2*phi)
    		    }else{ // subevent B
      		        QxB += (px.at(k).at(i)*px.at(k).at(i)-py.at(k).at(i)*py.at(k).at(i))/pt.at(k).at(i) ;  //which is pt*cos(2*phi)
      		        QyB += 2.*px.at(k).at(i)*py.at(k).at(i)/pt.at(k).at(i) ;       // which is pt*cos(2*phi)
    		    }
  	        } // end if

        }

        if(_nv2==0) continue ;
  	    const double psi2A = 0.5*atan2(QyA,QxA) ;
  	    const double psi2B = 0.5*atan2(QyB,QxB) ;
	
	    Rn += cos(2.0*(psi2A-psi2B)) ;

    }//event loop

    Rn = sqrt(Rn/nevents) ;

    cout << "Rn^sub = " << Rn << endl;

 	const double pf = sqrt(TMath::Pi())/(2.0*sqrt(2.0)) ;
 	double ksiMin=0., ksiMax = 20. ;

 	while(ksiMax-ksiMin>0.01){

   	    double ksi = 0.5*(ksiMin+ksiMax) ;
   	    double R = pf*ksi*exp(-0.25*ksi*ksi)*(TMath::BesselI0(0.25*ksi*ksi)+TMath::BesselI1(0.25*ksi*ksi)) ;
   	    if(R>Rn) ksiMax = ksi ;
   	    else ksiMin = ksi ;
 	}

	const double ksi = sqrt(2)*0.5*(ksiMin+ksiMax) ;

    cout << "ksi = " << ksi << endl;

 	Rn = pf*ksi*exp(-0.25*ksi*ksi)*(TMath::BesselI0(0.25*ksi*ksi)+TMath::BesselI1(0.25*ksi*ksi)) ;

    // plotting
    ofstream fout;
    fout.open("vn_EP_pT.dat", ofstream::app);
    fout << direct << endl;

    for(int ipt=1; ipt<havcos2ch->GetNbinsX()+1; ipt++){

        sprintf(hname,"v2pi_pt_avcos%i_%i",iplot,ipt) ;
        TH1D* hvch = havcos2ch->ProjectionY(hname,ipt,ipt) ;
        ptBin[ipt]= havcos2ch->GetXaxis()->GetBinCenter(ipt) ;
        v2[ipt] = hvch->GetMean()/Rn ;
        cout << setw(14) << ptBin[ipt] << setw(14) << v2[ipt] << endl ;
        fout << ptBin[ipt] << "\t" << v2[ipt] << endl ;
        delete hvch ;
    }

    fout << endl;
    fout.close();

}

/*int main(){

    vn_EP_pT("/Users/jakub/Desktop/", 100, 100);

}*/
