{
//create maps to store histograms and fitstograms 	
map<string,TF1*> fits; //fits for peak 1	 
map<string,TF1*> fits_;	//fits for peak 2
map<int,TF1*> eufits;	
map<int,TF1*> thfits;
map<string, TH1F*> h;

//arrays for timing info
int time_1[5]={};
int time_2[5]={};

//actual values of energy and empty channel one to fill for tgraph
float energy[10]={511.0,1173.2,344.0,1277.0,1332.5,779.0,224.0,964.0,1112.0, 1408.0};//kev
float energy_er[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};//kev
float channel[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};//kev
float channel_er[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};//kev

//arrays for auto-calculating fit parameters from 511
float zeroes[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
float ones[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
float twos[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
float eur[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

//set 511 values by just looking lol---------------------------
float z=355.0;
float o=8.0;
float t=0.3;

//automise naming into array
string directory="scinti_data/";
string names[5]={"na", "co", "eu", "ug", "th"};
string labels[5]={"Sodium-22", "Cobalt-60", "Europium-152", "Background", "Thorium"};
string file_type = ".TKA";
string file_names[5]={};
for(int i =0; i<5; i++){
file_names[i]=directory+names[i]+file_type;
}

//open fstream
ifstream in;

//find number of bins
string s; int sTotal; int bins;
in.open(file_names[0]);
while(!in.eof()) {	getline(in, s);	sTotal ++;	}
in.close();	bins = sTotal -1;

//estimate all other fit parameters from na first peak
zeroes[0]=z; ones[0]=o; twos[0]=t; //first na peak
zeroes[5]=z*2.4;ones[5]=o*1.2;twos[5]=t*0.14; //second na peak 
zeroes[1]=z*2.23; ones[1]=o*1.5;twos[1]=t*0.24;//first cobalt peak
zeroes[6]=z*2.52; ones[6]=o*1.3;twos[6]=t*0.2;//second co peak
zeroes[2]=z*0.685; ones[2]=o*0.8;twos[2]=t*4;//first eu peak
zeroes[7]=z*1.5; ones[7]=o*0.9;twos[7]=t;//second eu peak
//all extra europium peaks
eur[0]=z*0.49; eur[1]=z*1.84; eur[2]=z*2.08; eur[3]=z*2.66; eur[4]=o*0.5; eur[5]=o; eur[6]=o; eur[7] =o;
eur[8]=t*2.6;eur[9]=t*1.2; eur[10]=t*1.3; eur[11]=t;



//make the canvas
TCanvas c1("spectra","spectra");
c1.Divide(2,4); // divides the canvas into three rows and three columns

for(int i =0; i<3; i++){
//intf("current i: %i \n", i);
	
c1.cd(i+1); //use different pad for each spectrum
cout << names[i] <<endl;
h[names[i]]=new TH1F(labels[i].c_str(),labels[i].c_str(),bins,1,bins);
//h[names[i]]->Sumw2();
in.open(file_names[i]);

//Fill histos
double raw;
int t_1,t_2;
int nlines=0;
while(!(in.eof() /*|| inbg.eof()*/)){
    if(in>>raw /*&& inbg>>bg*/){
    nlines++;
    //cout <<nlines << endl;
    if (nlines==1) {
        t_1=raw;
        time_1[i]=t_1;
	continue;
    }
    if (nlines==2) {
	    t_2=raw;
        time_2[i]=t_2;
	continue;}
    //Put Data into histos
    h[names[i]]->SetBinContent(nlines-2,raw/t_2);
    //Do error calculation 
    double sRaw=TMath::Sqrt(1.0*raw)/t_2;
    //Set Errors to histos
    h[names[i]]->SetBinError(nlines-2,sRaw);
    }

    //printf("nlines: %i \n", nlines);
}
in.close();
//Define Fitfunc
//TF1 *fitfunc=new TF1("fitfunc","[2]*TMath::Gaus(x,[0],[1])");
fits[names[i]]=new TF1(names[i].c_str(),"[2]*TMath::Gaus(x,[0],[1])");
fits_[names[i]]=new TF1(names[i].c_str(),"[2]*TMath::Gaus(x,[0],[1])");
//fitfunc->SetParameters(350,10,2000);
fits[names[i]]->SetParameters(zeroes[i], ones[i], twos[i]);
fits_[names[i]]->SetParameters(zeroes[i+5], ones[i+5], twos[i+5]);

channel[i]=fits[names[i]]->GetParameter(0);
channel_er[i]=fits[names[i]]->GetParError(0);
channel[i+3]=fits_[names[i]]->GetParameter(0);
channel_er[i+3]=fits_[names[i]]->GetParError(0);

printf("Trying to get parameter for %i gives %f \n", i,fits[names[i]]->GetParameter(0));
printf("Trying to get 2nd  parameter for %i gives %f \n", i,fits_[names[i]]->GetParameter(0));
h[names[i]]->Fit(fits[names[i]],"","",zeroes[i]-2.5*ones[i],zeroes[i]+2.5*ones[i]);
h[names[i]]->Fit(fits_[names[i]],"Q+","",zeroes[i+5]-2.5*ones[i],zeroes[i+5]+2.5*ones[i+5]);
h[names[i]]->Draw("E");
}

//Extra fits on eu
for(int eu_i=0; eu_i<4; eu_i++){
	eufits[eu_i]=new TF1(Form("eu_%i",eu_i),"[2]*TMath::Gaus(x,[0],[1])");
	eufits[eu_i]->SetParameters(eur[eu_i], eur[eu_i+4], eur[eu_i+8]);
	h["eu"]->Draw();
	h["eu"]->Fit(Form("eu_%i",eu_i),"Q+","",eur[eu_i]-3*eur[eu_i+4],eur[eu_i]+3*eur[eu_i+4]);
	channel[eu_i+6]=eufits[eu_i]->GetParameter(0);
	channel_er[eu_i+6]=eufits[eu_i]->GetParError(0);

}
//c1.cd(4); //use different pad for each spectrum
//cout << names[i] <<endl;
for(int ith=3;ith<5;ith++){
h[names[ith]]=new TH1F(labels[ith].c_str(),labels[ith].c_str(),bins,1,bins);
in.open(file_names[ith]);
c1.cd(ith+1);

//Fill histos
double raw;
int t_1,t_2;
int nlines=0;
while(!(in.eof() /*|| inbg.eof()*/)){
    if(in>>raw /*&& inbg>>bg*/){
    nlines++;
    //cout <<nlines << endl;
    if (nlines==1) {
        t_1=raw;
        time_1[ith]=t_1;
	continue;
    }
    if (nlines==2) {
	    t_2=raw;
        time_2[ith]=t_2;
	continue;}
    //Put Data into histos
  h[names[ith]]->SetBinContent(nlines-2,raw/t_2);
    //Do error calculation 
    double sRaw=TMath::Sqrt(1.0*raw)/t_2;
    //Set Errors to histos
h[names[ith]]->SetBinError(nlines-2,sRaw);
    }

    //printf("nlines: %i \n", nlines);
}
in.close();
//h[names[ith]]->Draw();
}
h["th_minus"]=new TH1F("Thorium_no_back","Thorium_no_back",bins,1,bins);
//float scale = -time_1[4]/time_1[3];
for(int bin=0;bin<bins;bin++){
float err_th = h["th"]->GetBinError(bin);
float err_ug = h["ug"]->GetBinError(bin);//*scale;
float new_bin_err = TMath::Sqrt(TMath::Power(err_th,2)+TMath::Power(err_ug,2));	
float new_bin_c = h["th"]->GetBinContent(bin)-h["ug"]->GetBinContent(bin);
h["th_minus"]->SetBinError(bin, new_bin_err);
h["th_minus"]->SetBinContent(bin, new_bin_c);

}

//h["th_minus"]->Sumw2();
h["th_minus_log"]=h["th_minus"];
c1.cd(4);
h["th_minus"]->Draw();
c1.cd(5);
gPad->SetLogy(1);
h["th_minus_log"]->Draw();





// Fit Thorium boissss
float th_z[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
float th_o[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
float th_t[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


th_z[0]=z*0.175;
th_z[1]=z*0.315;
th_z[2]=z*0.68;
th_z[3]=z*0.8;
th_z[4]=z*1.15;
th_z[5]=z*1.6;
th_z[6]=z*0.19;
th_z[7]=z*4.89;

th_o[0]=o*0.4;
th_o[1]=o*0.4;
th_o[2]=o*0.5;
th_o[3]=o*0.5;
th_o[4]=o;
th_o[5]=o;
th_o[6]=o*0.4;
th_o[7]=o*2; 

th_t[0]=t*10;
th_t[1]=t*40;
th_t[2]=t*38;
th_t[3]=t*20;
th_t[4]=t*0.23;
th_t[5]=t*0.3;
th_t[6]=t*8; 
th_t[7]=t*0.15; 

for (int thi =1; thi<8; thi++){
	thfits[thi]=new TF1(Form("th_%i",thi),"[2]*TMath::Gaus(x,[0],[1])");
	thfits[thi]->SetParameters(th_z[thi], th_o[thi], th_t[thi]);
	//if(thi== 0){h["th_minus"]->Fit(Form("th_%i",thi),"Q+","",th_z[thi]-2*th_o[thi],th_z[thi]+2*th_o[thi]);}
	std::cout << Form("th_%i",thi) << std::endl;
	if(thi < 5){h["th_minus"]->Fit(Form("th_%i",thi),"Q+","",th_z[thi]-3*th_o[thi],th_z[thi]+3*th_o[thi]);}
	if(thi>4 && thi != 6){h["th_minus_log"]->Fit(Form("th_%i",thi),"Q+","",th_z[thi]-3*th_o[thi],th_z[thi]+3*th_o[thi]);}

}
th_z[8]=z*0.479; 
th_o[8]=o*0.5; 
th_t[8]=t*8;

th_z[9]=z*0.535; 
th_o[9]=o*0.5; 
th_t[9]=t*8;


thfits[8]=new TF1("mega","[2]*TMath::Gaus(x,[0],[1])+[5]*TMath::Gaus(x,[3],[4])");
thfits[8]->SetParameters(th_z[8], th_o[8], th_t[8], th_z[9], th_o[9], th_t[9]);
h["th_minus"]->Fit("mega","Q+","",th_z[8]-2.5*th_o[8],th_z[9]+2.5*th_o[9]);

thfits[0]=new TF1("mega2","[2]*TMath::Gaus(x,[0],[1])+[5]*TMath::Gaus(x,[3],[4])");
thfits[0]->SetParameters(th_z[0], th_o[0], th_t[0], th_z[6], th_o[6], th_t[6]);
h["th_minus"]->Fit("mega2","Q+","",th_z[0]-2.8*th_o[0],th_z[6]+2.5*th_o[6]);

auto gr = new TGraphErrors(10);//nergy[n_points],channel[n_points],energy_er[n_points],channel_er[n_points]);

for(int n=0;n<10;n++){
energy_er[n]=100;
channel_er[n]=100;
gr->SetPoint(n,energy[n], channel[n]);
gr->SetPointError(n,energy_er[n], channel_er[n]);
//printf("Setting point %i with e %f and channel %f \n", n, energy[n], channel[n]);
}

gr->SetTitle("TGraphErrors");
gr->SetMarkerColor(38);
gr->SetMarkerStyle(1);
gr->SetMarkerSize(5);
gr->GetXaxis()->SetTitle("Energy (kev)");
gr->GetYaxis()->SetTitle("Channel");
gr->Fit("pol1", "Q");

c1.cd(6);
gr->Draw("ap");

//////////////////////////////try a tgraph on the peaks man to have a lil check
//int range = o*2;
//auto grna = new TGraphErrors(range);//nergy[n_points],channel[n_points],energy_er[n_points],channel_er[n_points]);
//for(int n=0;n<range;n++){
//grna->SetPoint(n,z-o+n,h["na"]->GetBinContent(z-o+n));
//grna->SetPointError(n,5, h["na"]->GetBinError(z-o+n));
//printf("Setting point %i with e %f and channel %f \n", n, energy[n], channel[n]);
//}
//grna->SetTitle("Na peak as TGraphErrors");
//grna->SetMarkerColor(38);
//grna->SetMarkerStyle(1);
//grna->SetMarkerSize(5);
//grna->GetXaxis()->SetTitle("Channel");
//grna->GetYaxis()->SetTitle("Counts");
//grna->Fit("gaus");
//c1.cd(7);
//grna->Draw("ap");

}
