void fit3()
{

   double B_M;
   double q2;
   double q2_NB;
   double ctl_lhcb;
   double ctk_lhcb;
   double phi_lhcb;

   TChain* chain = new TChain("DecayTree");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2011/Flagged/B02Kst0Jpsi2eeSS-Data-2011-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2011/Flagged/B02Kst0Jpsi2eeSS-Data-2011-MagUp-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2012/Flagged/B02Kst0Jpsi2eeSS-Data-2012-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2012/Flagged/B02Kst0Jpsi2eeSS-Data-2012-MagUp-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2015/Flagged/B02Kst0Jpsi2eeSS-Data-2015-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2015/Flagged/B02Kst0Jpsi2eeSS-Data-2015-MagUp-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2016/Flagged/B02Kst0Jpsi2eeSS-Data-2016-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2016/Flagged/B02Kst0Jpsi2eeSS-Data-2016-MagUp-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2017/Flagged/B02Kst0Jpsi2eeSS-Data-2017-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2017/Flagged/B02Kst0Jpsi2eeSS-Data-2017-MagUp-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2018/Flagged/B02Kst0Jpsi2eeSS-Data-2018-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2018/Flagged/B02Kst0Jpsi2eeSS-Data-2018-MagUp-Translated-Flagged-XGBOutput-q2BDTOutput.root");

   Long64_t GenericPresel , GenericPresel_Additional, MeerkatPresel_Tight, PIDPresel, TighterKst0Presel, TriggerPresel, VetoesPresel, VetoesPresel_Additional, CloneVeto;
   Float_t XGBOutput;
   chain->SetBranchStatus("*", 0);
   chain->SetBranchStatus("B_M"     ,1);
   chain->SetBranchStatus("q2"     ,1);
   chain->SetBranchStatus("q2_NB"     ,1);
   chain->SetBranchStatus("ctl_lhcb"     ,1);
   chain->SetBranchStatus("ctk_lhcb"     ,1);
   chain->SetBranchStatus("phi_lhcb"     ,1);
   chain->SetBranchAddress("B_M"     ,&B_M);
   chain->SetBranchAddress("q2"     ,&q2);
   chain->SetBranchAddress("q2_NB"     ,&q2_NB);
   chain->SetBranchAddress("ctl_lhcb"     ,&ctl_lhcb);
   chain->SetBranchAddress("ctk_lhcb"     ,&ctk_lhcb);
   chain->SetBranchAddress("phi_lhcb"     ,&phi_lhcb);
   chain->SetBranchStatus("GenericPresel",1);
   chain->SetBranchStatus("GenericPresel_Additional",1);
   chain->SetBranchStatus("MeerkatPresel_Tight",1);
   chain->SetBranchStatus("PIDPresel",1);
   chain->SetBranchStatus("TighterKst0Presel",1);
   chain->SetBranchStatus("TriggerPresel",1);
   chain->SetBranchStatus("VetoesPresel",1);
   chain->SetBranchStatus("VetoesPresel_Additional",1);
   chain->SetBranchStatus("CloneVeto",1);
   chain->SetBranchStatus("XGBOutput",1);
   chain->SetBranchAddress("GenericPresel"     ,&GenericPresel);
   chain->SetBranchAddress("GenericPresel_Additional"     ,&GenericPresel_Additional);
   chain->SetBranchAddress("MeerkatPresel_Tight"     ,&MeerkatPresel_Tight);
   chain->SetBranchAddress("PIDPresel"     ,&PIDPresel);
   chain->SetBranchAddress("TighterKst0Presel"     ,&TighterKst0Presel);
   chain->SetBranchAddress("TriggerPresel"     ,&TriggerPresel);
   chain->SetBranchAddress("VetoesPresel"     ,&VetoesPresel);
   chain->SetBranchAddress("VetoesPresel_Additional"     ,&VetoesPresel_Additional);
   chain->SetBranchAddress("CloneVeto"     ,&CloneVeto);
   chain->SetBranchAddress("XGBOutput"     ,&XGBOutput);
   int n_bins = 100;
   TH1D* h_B_M = new TH1D("h_B_M", ";m(K^{*0}e^{-}e^{+}) [MeV];Candidates",n_bins,4500,7000);

   for (unsigned int i=0; i<chain->GetEntries(); ++i)
   {
        chain->GetEntry(i);
        if((q2_NB>14.0)&& (GenericPresel) && (GenericPresel_Additional) && (MeerkatPresel_Tight) && (PIDPresel) && (TighterKst0Presel) && (TriggerPresel) && (VetoesPresel) && (VetoesPresel_Additional) && (CloneVeto))
        {
          if(XGBOutput>0.01)
          {

            h_B_M->Fill(B_M);
          }
        }

  }

    TCanvas *canvas = new TCanvas();
    TF1* f = new TF1("f","gauss");
    f->SetLineColor(kRed);
    gStyle->SetOptFit(11);
    h_B_M->Fit(f);
    h_B_M->Draw("AP");
    h_B_M->SetMinimum(0);
    h_B_M->SetLineColor(kBlue);
    h_B_M->SetStats(0);
    canvas->Draw();
    TString filename = (TString)h_B_M->GetName() + ".pdf";
    canvas->Print(filename);

   }
