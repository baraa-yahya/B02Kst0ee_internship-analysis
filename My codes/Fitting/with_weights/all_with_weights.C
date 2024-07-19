void background()
{
   double B_M;
   double q2;
   double q2_NB;
   double ctl_lhcb;
   double ctk_lhcb;
   double phi_lhcb;
   Int_t B0_BKGCAT;
   double w_Track_Trig_Kin;
   Long64_t GenericPresel , GenericPresel_Additional, MeerkatPresel_Tight, PIDPresel, TighterKst0Presel, TriggerPresel, VetoesPresel, VetoesPresel_Additional, CloneVeto;
   Float_t XGBOutput;

   TChain* chain = new TChain("DecayTree");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2011/Flagged/Bu2Kpipiee-MC-2011-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2011/Flagged/Bu2Kpipiee-MC-2011-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2012/Flagged/Bu2Kpipiee-MC-2012-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2012/Flagged/Bu2Kpipiee-MC-2012-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2015/Flagged/Bu2Kpipiee-MC-2015-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2015/Flagged/Bu2Kpipiee-MC-2015-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2016/Flagged/Bu2Kpipiee-MC-2016-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2016/Flagged/Bu2Kpipiee-MC-2016-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2017/Flagged/Bu2Kpipiee-MC-2017-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2017/Flagged/Bu2Kpipiee-MC-2017-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2018/Flagged/Bu2Kpipiee-MC-2018-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/Bu2Kpipiee/2018/Flagged/Bu2Kpipiee-MC-2018-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2011/Flagged/B02Kst0eeSignal-MC-2011-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2011/Flagged/B02Kst0eeSignal-MC-2011-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2012/Flagged/B02Kst0eeSignal-MC-2012-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2012/Flagged/B02Kst0eeSignal-MC-2012-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2015/Flagged/B02Kst0eeSignal-MC-2015-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2015/Flagged/B02Kst0eeSignal-MC-2015-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2016/Flagged/B02Kst0eeSignal-MC-2016-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2016/Flagged/B02Kst0eeSignal-MC-2016-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2017/Flagged/B02Kst0eeSignal-MC-2017-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2017/Flagged/B02Kst0eeSignal-MC-2017-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2018/Flagged/B02Kst0eeSignal-MC-2018-MagDown-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");
   chain2->Add("/sps/lhcb/pietrzyk/analysis/samples/MC/B02Kst0eeSignal/2018/Flagged/B02Kst0eeSignal-MC-2018-MagUp-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root");

   chain->SetBranchStatus("*", 0);
   chain->SetBranchStatus("B_M"     ,1);
   chain->SetBranchStatus("q2"     ,1);
   chain->SetBranchStatus("ctl_lhcb"     ,1);
   chain->SetBranchStatus("ctk_lhcb"     ,1);
   chain->SetBranchStatus("phi_lhcb"     ,1);
   chain->SetBranchStatus("B0_BKGCAT",1);

   // Set branch status for w_Track_Trig_Kin
  chain->SetBranchStatus("w_Track_Trig_Kin", 1);
  // Declare variable to hold w_Track_Trig_Kin
  chain->SetBranchAddress("w_Track_Trig_Kin", &w_Track_Trig_Kin);

   chain->SetBranchAddress("B0_BKGCAT"      , &B0_BKGCAT);
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
   TH1D* h_B_M_part_rec = new TH1D("h_B_M", ";m(K^{*0}e^{-}e^{+}) [MeV];Candidates",n_bins,4500,6000);
   TH1D* h_B_M_part_rec_with_weight = new TH1D("h_B_M_with_weight", ";m(K^{*0}e^{-}e^{+}) [MeV];Candidates",n_bins,4500,6000);
   
   TH1D* h_B_M0 = new TH1D("h_B_M0", ";m(K^{*0}e^{-}e^{+}) [MeV];Candidates",n_bins,4500,6000);



   vector<TH1D*> vec_hist = {
     h_B_M,
 
  };
  vector<TH1D*> vec_hist_with_weight = {
     h_B_M_with_weight,

  };


   for (unsigned int i=0; i<chain->GetEntries(); ++i)
   {
        chain->GetEntry(i);
        if(q2_NB>14.0)
        {
          if(((GenericPresel) && (GenericPresel_Additional) && (MeerkatPresel_Tight) && (PIDPresel) && (TighterKst0Presel) && (TriggerPresel) && (VetoesPresel) && (VetoesPresel_Additional) && (CloneVeto) && (XGBOutput>0.1)) &&((B0_BKGCAT)==0 || (B0_BKGCAT)==10 || (B0_BKGCAT)==50))
          {
            h_B_M->Fill(B_M);
            h_B_M_with_weight->Fill(B_M, w_Track_Trig_Kin);
          }
        }
  }

  for (unsigned int i=0; i<vec_hist.size(); ++i)
  {
    TCanvas *canvas = new TCanvas();
    vec_hist[i]->SetMinimum(0);
    vec_hist_with_weight[i]->SetMinimum(0);
    vec_hist[i]->SetLineColor(kRed);
    vec_hist_with_weight[i]->SetLineColor(kBlue);
    vec_hist[i]->SetStats(0);
    vec_hist_with_weight[i]->SetStats(0);
    vec_hist_with_weight[i]->Scale(1.0*vec_hist[i]->Integral()/vec_hist_with_weight[i]->Integral());
    vec_hist[i]->Draw();
    vec_hist_with_weight[i]->Draw("same");
    cout<< vec_hist[i]->Integral() << endl;
    cout<< vec_hist_with_weight[i]->Integral()<<endl;
    canvas->Draw();
    TString filename = (TString)vec_hist[i]->GetName() + ".pdf";
    canvas->Print(filename);
  }
// Delete histograms
    for (auto hist : vec_hist)
        delete hist;
    for (auto hist : vec_hist_with_weight)
        delete hist;    
}












