{

//open Files
ifstream in, inbg;
in.open("scinti_data/07_Na22_01.TKA");
//inbg.open("../Testmessung/run010_time.TKA");


//Define histos
//
TH1F *hRaw = new TH1F("Time Spectrum","Time Spectrum;channel;counts",4096,1,4096);
TH1F *hBg = new TH1F("Background","Background",4096,1,4096);
TH1F *hBeta = new TH1F("beta spectrum","beta spectrum;channel;counts",4096,1,4096);

//Fill histos
//
double raw,bg;
double t,tbg;
int nlines=0;
while(!(in.eof() /*|| inbg.eof()*/)){
    if(in>>raw /*&& inbg>>bg*/){
    nlines++;
    cout <<nlines << endl;
    if (nlines==1) {
        t=raw;
        //tbg=bg;
        printf("Meas Duration Spectrum: %.1f s and Backgroung: %.1f s\n",t,tbg);
    }
    if (nlines<3) continue; //Skip lines with time information
    //Put Data into histos
    hRaw->SetBinContent(nlines-2,raw);
    //hBg->SetBinContent(nlines-2,bg);
    hBeta->SetBinContent(nlines-2,1.0*raw/t/*-1.0*bg/tbg*/);
    //Do error calculation 
    //
    double sRaw=TMath::Sqrt(1.0*raw)/t;
    double sBg=TMath::Sqrt(1.0*bg)/tbg;
    double sBeta=TMath::Sqrt(/*1.0*sBg*sBg+*/1.0*sRaw*sRaw);
    
    //Set Errors to histos
    //
    hRaw->SetBinError(nlines-2,sRaw);
    //hBg->SetBinError(nlines-2,sBg);
    hBeta->SetBinError(nlines-2,sBeta);
    }
}
//Define Fitfunc
TF1 *fitfunc=new TF1("fitfunc","[0]*exp(-x/[1])+[2]");
fitfunc->SetParameters(30,1000,0);

hRaw->Fit("fitfunc","","",40,200);

//time calibration
//
double T[4]={2.44,4.57,6.32,8.6};
double Ch[4]={895.7,1772.3,2508.97,3506.0};
double sT[4]={0.1,0.1,0.1,0.1};
double sCh[4]={1,1,1,1};


TGraphErrors *ECalibGraph=new TGraphErrors(4,Ch,T,sCh,sT);
ECalibGraph->SetTitle("Time Calibration");
TF1* TCalibFunction=new TF1("TCalib","[0]*x+[1]");
TCanvas c2;
c2.cd();
ECalibGraph->Fit("TCalib");
ECalibGraph->Draw();

//conversion to µs
//
double TCh=fitfunc->GetParameter(1);
double sTCh=fitfunc->GetParError(1);

double TMu=TCalibFunction->GetParameter(0)*TCh;

double slopeSquare=TCalibFunction->GetParameter(0)*TCalibFunction->GetParameter(0);
double sslopeSquare=TCalibFunction->GetParError(0)*TCalibFunction->GetParError(0);
double TChSquare=TCh*TCh;


double sTMu=TMath::Sqrt((sslopeSquare/slopeSquare +sTCh*sTCh/TChSquare))*TCalibFunction->GetParameter(0)*TCh;

printf("tau is (%f +- %f)Ch \n",TCh,sTCh);
printf("tau is (%f +- %f)mus \n",TMu,sTMu);



}
