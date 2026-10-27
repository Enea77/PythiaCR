// jetshape.cpp -- correction chain and jet shapes for the jet-V0 correlations of correlations/jetV0_analyzer.cpp.
// Port of JetAnalyzer/jetshape_V0_rebinned.cpp (the nominal version, driven by jetshape_loop.cpp) with the SAME
// algorithms (jet_shape_V0.h copied unchanged). Changes only where the new analyzer output differs:
//   - ptw = true reads the pT-weighted correlations (histogram names with suffix "_ptw"), output prefix "ptw_"
//   - jet sparses have multiplicity on axis 3 (was 4; correlations still axis 4), see CorrCfg::Ax
//   - no dependence on input_variables.h (none of its globals was used)
// Stages: raw -> signal_over_mixing -> seagull_corrected / sideband_bkg_corrected -> fully_corrected
//
//   root -l -b -q 'jetshape.cpp("reco","reco","MC_Lam","input/",true,60,2,-1,-1,false)'
//   ptEdges (last argument), e.g. "3,10": one V0 pT bin 3-10 GeV (output jetshapes_rebinned_ptbins3to10_...)
#include "jet_shape_V0.h"
#include "TObjString.h"
#include "injet_eff.h"

// UE (sideband) region used in the background subtraction: gUEMin < |Delta eta| < gUEMax (nominal 1.5-2.5).
// Systematic variations (run_jetshapes_data.C argument ue = "1.5,2.0"): output tag "ue<min>to<max>_".
float gUEMin = 1.5, gUEMax = 2.5;

void jetshape(string jet_type = "reco", string V0_type = "gen",
                 TString filename = "K0s_Jptmin80_nojettrig_VzReweighted_pthat50andabove_WTA_GenV0",
                 TString filepath = "files_output/MC_Unembedded/",
                 bool remove_mass_bkg = true,
                 float jetptmin = 80,
                 int bkg_method = 2, // 1 = Method A (subtract before mix), 2 = Method B (subtract after mix)
                 int multmin = -1,
                 int multmax = -1,
                 bool ptw = false,         // true: pT-weighted correlations (names with suffix _ptw)
                 TString ptEdges = "",     // V0 pT bin edges, e.g. "3,10" (default: 1.5,2,3,5,7,10); tag "ptbins3to10_" in the output name
                 TString injetEff = ""     // reco V0 only: in-jet efficiency file (make_injet_eff.C, injet_eff.h); the corrected
                                           // correlation (before the UE subtraction) is divided by c(pT, dr); tag "ijc_" in the output name
                 )
{
    InJet::Maps ijMaps;
    // optional systematic variation of the in-jet factors: injetEff = "<file>:up" / ":down" / ":half" (tag ijcup_ / ijcdown_ / ijchalf_)
    TString ijVar = "";
    InJet::gVar = 0;
    if (injetEff.Contains(":")) {
        ijVar = injetEff(injetEff.Last(':') + 1, injetEff.Length()); injetEff = injetEff(0, injetEff.Last(':'));
        InJet::gVar = ijVar == "up" ? 1 : ijVar == "down" ? -1 : ijVar == "half" ? 2 : ijVar == "fit" ? 3 : ijVar == "fitup" ? 4 : ijVar == "fitdown" ? 5 : 0;
        if (InJet::gVar == 0) { cout << "----Error: unknown in-jet variation " << ijVar << "----" << endl; return; }
    }
    const bool doInJet = (injetEff != "" && V0_type == "reco");
    if (doInJet) { ijMaps = InJet::Load(injetEff); if (!ijMaps.gen) { cout << "----Error: in-jet efficiency file not usable----" << endl; return; } }
    std::vector<double> pt_edges = {1.5,2, 3.0,5,7,10.0};
    TString pt_tag = "";
    if (ptEdges != "") {
        pt_edges.clear();
        TObjArray* tok = ptEdges.Tokenize(",");
        for (int k = 0; k < tok->GetEntries(); ++k) pt_edges.push_back(((TObjString*)tok->At(k))->GetString().Atof());
        delete tok;
        pt_tag = Form("ptbins%gto%g_", pt_edges.front(), pt_edges.back());
        if (pt_edges.size() > 2) pt_tag = Form("ptbins%gto%gn%zu_", pt_edges.front(), pt_edges.back(), pt_edges.size() - 1);
    }

    TFile *f1 = TFile::Open(filepath + filename + ".root");
    if (!f1 || f1->IsZombie()) { cout << "----Error: Input file not found!----" << endl; return; }

    TString method_str = "";
    if (remove_mass_bkg) method_str = (bkg_method == 1) ? "MethodA_" : "";
    //rebinned_
    const TString sfx = ptw ? "_ptw" : "";
    const TString ue_tag = (std::fabs(gUEMin - 1.5) > 1e-4 || std::fabs(gUEMax - 2.5) > 1e-4) ? TString(Form("ue%gto%g_", gUEMin, gUEMax)) : TString("");
    TString output_filename = Form(filepath+"jetshapes_rebinned_"+ue_tag+pt_tag+(doInJet ? TString("ijc") + ijVar + "_" : TString(""))+(ptw ? TString("ptw_") : TString(""))+(multmin>-1? Form("mult%.1ito%.1i_",multmin,multmax):"")+
                        (remove_mass_bkg ? method_str.Data() : "no_bkg_removal_")+
                        "%s_%s_",jet_type.c_str(),V0_type.c_str())+filename+".root";
    TFile *f2 = new TFile(output_filename, "RECREATE");
    if (!f2 || f2->IsZombie()) {
        cout << "----Error: Could not create output file: " << output_filename.Data() << "----" << endl;
        if (f2) {
            f2->Close();
            delete f2;
        }
        return;
    }

    int mult = 0;

    //common selections
    pair<float, float> signal_deta = make_pair(-1.0, 1.0);   
    pair<float, float> fullranges = make_pair(-10,10);   

    // Jussi selections
    pair<float, float> jussi_nearside = make_pair(-1.0, 1.0);
    pair<float, float> jussi_awayside = make_pair(PI-1, PI+1);
    pair<float, float> jussi_betweenpeaks = make_pair(1.4, 1.7);
    
    // Olga selections
    pair<float, float> olga_nearside = make_pair(-0.1, 0.1);
    pair<float, float> olga_awayside = make_pair(PI-0.1, PI+0.1);
    pair<float, float> olga_betweenpeaks = make_pair(PI/2.-0.1, PI/2+0.1); 
       
    bool dooneregion = false;
    bool dobinnorm = false;

    //input vector of histograms -> 2D
    std::vector<TH2 *> rawsignal;
    std::vector<TH2 *> mixing;
    std::vector<TH2 *> rawbkg; // NEW: mass sideband signal
    std::vector<TH2 *> mixingbkg; // NEW: mass sideband mixing
    
    std::vector<TH2 *> normsignal;
    std::vector<TH2 *> normmixing;
    std::vector<TH2 *> rawsignal_after_mixing_correction;
    std::vector<TH2 *> signal_after_seagull_correction;
    std::vector<TH2 *> signal_after_sideband_bkg_correction;
    std::vector<TH2 *> signal_after_seagull_and_sideband_bkg_correction;
    std::vector<TH2 *> raw_seagull;
    
    std::vector<TH1 *> jet_shape_histo_raw; 
    std::vector<TH1 *> jet_shape_histo_signal_over_mixing; 
    std::vector<TH1 *> jet_shape_histo_seagull_corrected; 
    std::vector<TH1 *> jet_shape_histo_sideband_bkg_corrected; 
    std::vector<TH1 *> jet_shape_histo_fully_corrected; 

    // method B components for the process plots (plot_process.C), reco V0 only; one entry per pT bin
    std::vector<TH2 *> proc_sig_over_mix, proc_norm_bkg, proc_norm_mix_bkg, proc_bkg_over_mix;

    // V0 pT bin edges: pt_edges (set at the top from ptEdges)
    const int New_ptbins = pt_edges.size() - 1;

    TH1D* Integral_vs_pt_bin_widths_scale = new TH1D("Integral_vs_pt_bin_widths_scale","Integral_vs_pt_bin_widths_scale",New_ptbins,pt_edges.data());
    Integral_vs_pt_bin_widths_scale->Sumw2();
    for (int i=1; i<New_ptbins+1; i++) Integral_vs_pt_bin_widths_scale->SetBinContent(i,Integral_vs_pt_bin_widths_scale->GetBinWidth(i)/Integral_vs_pt_bin_widths_scale->GetBinWidth(1));

    float dr_start = 0;
    float dr_end = 0.99;
    Double_t Integral_error = 0;
    TH1D* Integral_vs_pt_raw = new TH1D("Integral_vs_pt_raw","Integral_vs_pt_raw",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_signal_over_mixing = new TH1D("Integral_vs_pt_signal_over_mixing","Integral_vs_pt_signal_over_mixing",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_seagull_corrected = new TH1D("Integral_vs_pt_seagull_corrected","Integral_vs_pt_seagull_corrected",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_sideband_bkg_corrected = new TH1D("Integral_vs_pt_sideband_bkg_corrected","Integral_vs_pt_sideband_bkg_corrected",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_fully_corrected = new TH1D("Integral_vs_pt_fully_corrected","Integral_vs_pt_fully_corrected",New_ptbins,pt_edges.data());

    TH1D* Integral_vs_pt_innercone_raw = new TH1D("Integral_vs_pt_innercone_raw","Integral_vs_pt_innercone_raw",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_innercone_signal_over_mixing = new TH1D("Integral_vs_pt_innercone_signal_over_mixing","Integral_vs_pt_innercone_signal_over_mixing",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_innercone_seagull_corrected = new TH1D("Integral_vs_pt_innercone_seagull_corrected","Integral_vs_pt_innercone_seagull_corrected",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_innercone_sideband_bkg_corrected = new TH1D("Integral_vs_pt_innercone_sideband_bkg_corrected","Integral_vs_pt_innercone_sideband_bkg_corrected",New_ptbins,pt_edges.data());
    TH1D* Integral_vs_pt_innercone_fully_corrected = new TH1D("Integral_vs_pt_innercone_fully_corrected","Integral_vs_pt_innercone_fully_corrected",New_ptbins,pt_edges.data());

    // 2. Update string formatting to look for _V0_ instead of _track_
    TString dir_name = Form("correlation_%s_%s_histograms", jet_type.c_str(), V0_type.c_str());
    THnSparse *signal_sparse = (THnSparse *)f1->Get(Form("%s/hist_correlation_signal_jet_%s_V0_%s%s", dir_name.Data(), jet_type.c_str(), V0_type.c_str(), sfx.Data()));
    THnSparse *mixing_sparse = (THnSparse *)f1->Get(Form("%s/hist_correlation_mixing_jet_%s_V0_%s%s", dir_name.Data(), jet_type.c_str(), V0_type.c_str(), sfx.Data()));
    
    if (!signal_sparse || !mixing_sparse) { cout<<"----Error: Correlation sparses not found!----"<<endl; return; }
    if (multmin > -1 && signal_sparse->GetAxis(4)->GetNbins() == 1) {
        cout<<"----Error: correlation sparses have one multiplicity bin (MC); use multmin = multmax = -1----"<<endl; return;
    }
    if (multmin > -1){
        signal_sparse->GetAxis(4)->SetRangeUser(multmin+0.1,multmax-0.1);
        mixing_sparse->GetAxis(4)->SetRangeUser(multmin+0.1,multmax-0.1);
    }

    THnSparse *jetpt_sparse = (THnSparse *)f1->Get(Form("QA_histograms/hist_gen_jet_weighted"));
    if (jet_type == "reco" || jet_type == "ref") jetpt_sparse = (THnSparse *)f1->Get(Form("correlation_%s_%s_histograms/hist_jet_from_%s_%s_sig",jet_type.c_str(), V0_type.c_str(),jet_type.c_str(), V0_type.c_str()));
    if (!jetpt_sparse || jetpt_sparse->GetEntries() < 10) {cout<<"----Error: No jets found in QA histograms!----"<<endl; return;}

    // Apply Jet pT cut to QA projection
    jetpt_sparse->GetAxis(0)->SetRangeUser(jetptmin+0.0001,8000);
    if (multmin > -1) jetpt_sparse->GetAxis(3)->SetRangeUser(multmin+0.1,multmax-0.1);   // jet sparse: mult is axis 3
    TH1D *ptbins = (TH1D*)jetpt_sparse->Projection(0);

    THnSparse *bkg_sparse = nullptr;
    THnSparse *mixing_bkg_sparse = nullptr;

    // 4. apply background 
    if (V0_type == "reco" ) {
        bkg_sparse = (THnSparse *)f1->Get(Form("%s/hist_correlation_bkg_signal_jet_%s_V0_%s%s", dir_name.Data(), jet_type.c_str(), V0_type.c_str(), sfx.Data()));
        mixing_bkg_sparse = (THnSparse *)f1->Get(Form("%s/hist_correlation_bkg_mixing_jet_%s_V0_%s%s", dir_name.Data(), jet_type.c_str(), V0_type.c_str(), sfx.Data()));
        if (!bkg_sparse || !mixing_bkg_sparse) { cout<<"----Error: Sideband correlation sparses not found!----"<<endl; return; }           
        // Sparses are NO LONGER subtracted here. Handled in 2D loop below.
         if (multmin > -1){
            bkg_sparse->GetAxis(4)->SetRangeUser(multmin+0.1,multmax-0.1);
            mixing_bkg_sparse->GetAxis(4)->SetRangeUser(multmin+0.1,multmax-0.1);
        }
    }

    if (filename.Contains("MC")) {
        TH2D* pthat_2D_weighted = (TH2D* )f1->Get("QA_histograms/pthathist_weighted");
        if (pthat_2D_weighted) pthat_2D_weighted->ProjectionX()->Draw();
    }
   

    
// Loop over the number of custom bins (number of edges - 1)
    for (size_t i = 0; i < pt_edges.size() - 1; i++) {
        double pt_min = pt_edges[i];
        double pt_max = pt_edges[i+1];
        
        int bin_min = signal_sparse->GetAxis(2)->FindBin(pt_min);
        int bin_max = signal_sparse->GetAxis(2)->FindBin(pt_max - 1e-6);

        // 1. Extract Signal
        signal_sparse->GetAxis(2)->SetRange(bin_min, bin_max);
        TH2* fromsparse_2D_sig = (TH2*)signal_sparse->Projection(1,0);
        fromsparse_2D_sig->SetName(Form("rawsig_ptbin_%zu", i)); // <-- ADD THIS
        rawsignal.push_back(fromsparse_2D_sig);
        
        // 2. Extract Mixing
        mixing_sparse->GetAxis(2)->SetRange(bin_min, bin_max);
        TH2* fromsparse_2D_mix = (TH2*)mixing_sparse->Projection(1,0);
        fromsparse_2D_mix->SetName(Form("mixing_ptbin_%zu", i)); // <-- ADD THIS
        mixing.push_back(fromsparse_2D_mix);

        // 3. Extract Background (if applicable)
        if (V0_type == "reco" && remove_mass_bkg && bkg_sparse && mixing_bkg_sparse) {
            bkg_sparse->GetAxis(2)->SetRange(bin_min, bin_max);
            TH2* fromsparse_2D_bkg = (TH2*)bkg_sparse->Projection(1,0);
            fromsparse_2D_bkg->SetName(Form("rawbkg_ptbin_%zu", i)); // <-- ADD THIS
            rawbkg.push_back(fromsparse_2D_bkg);
            
            mixing_bkg_sparse->GetAxis(2)->SetRange(bin_min, bin_max);
            TH2* fromsparse_2D_mixbkg = (TH2*)mixing_bkg_sparse->Projection(1,0);
            fromsparse_2D_mixbkg->SetName(Form("mixingbkg_ptbin_%zu", i)); // <-- ADD THIS
            mixingbkg.push_back(fromsparse_2D_mixbkg);
        }
    }

    cout << "---------------------------------------- Final Histograms Code ----------------------------------------" << endl;

    for (int i = 0; i < rawsignal.size(); i++)
    {
        if(i!=0)cout << "Processing 2D histograms " << (float)100*i/rawsignal.size() << " %" << endl;
        
        TH2* raw_signal = (TH2*)rawsignal[i]->Clone(Form("raw_sig_%d", i));
        TH2* normalized_signal = (TH2*)raw_signal_normalization(raw_signal, ptbins, 3.0);
        normsignal.push_back(normalized_signal); 
        
        TH2* normalized_mixing = (TH2*)mixing_normalization(mixing[i], true, true, 0.2);
        normmixing.push_back(normalized_mixing); 
        
        TH2* normalized_signal_temp = (TH2*)normalized_signal->Clone(Form("norm_sig_temp_%d", i));
        TH2* normalized_mixing_temp = (TH2*)normalized_mixing->Clone(Form("norm_mix_temp_%d", i));
        
        TH2* corrected_via_mixing = nullptr;

        // Apply dynamic Background Subtraction Logic
        if (V0_type == "reco" && remove_mass_bkg && rawbkg.size() > i && mixingbkg.size() > i) {
            TH2* raw_bkg = (TH2*)rawbkg[i]->Clone(Form("raw_bkg_%d", i));
            TH2* normalized_bkg = (TH2*)raw_signal_normalization(raw_bkg, ptbins, 3.0);
            TH2* normalized_mixing_bkg = (TH2*)mixing_normalization(mixingbkg[i], true, true, 0.2);
            
            if (bkg_method == 1) {
                // Method A: (Signal_peak - Signal_bkg) / Mix_peak
                normalized_signal_temp->Add(normalized_bkg, -1.0);
                corrected_via_mixing = (TH2*)corrected_via_mixing_events(normalized_signal_temp, normalized_mixing_temp);
            } else if (bkg_method == 2) {
                // Method B: (Signal_peak / Mix_peak) - (Signal_bkg / Mix_bkg)
                // (corrected_via_mixing_events divides its first argument in place: clones kept before)
                proc_norm_bkg.push_back((TH2*)normalized_bkg->Clone(Form("proc_norm_bkg_%d", i)));
                proc_norm_mix_bkg.push_back((TH2*)normalized_mixing_bkg->Clone(Form("proc_norm_mix_bkg_%d", i)));
                corrected_via_mixing = (TH2*)corrected_via_mixing_events(normalized_signal_temp, normalized_mixing_temp);
                TH2* corrected_bkg = (TH2*)corrected_via_mixing_events(normalized_bkg, normalized_mixing_bkg);
                proc_sig_over_mix.push_back((TH2*)corrected_via_mixing->Clone(Form("proc_sig_over_mix_%d", i)));
                proc_bkg_over_mix.push_back((TH2*)corrected_bkg->Clone(Form("proc_bkg_over_mix_%d", i)));
                corrected_via_mixing->Add(corrected_bkg, -1.0);
                if (doInJet) {   // in-jet efficiency: after the mass-sideband subtraction, before the UE subtraction
                    std::vector<double> ce;
                    const std::vector<double> cf = InJet::Factors(ijMaps, pt_edges[i], pt_edges[i + 1], &ce);
                    InJet::Apply(corrected_via_mixing, cf);
                    cout << Form("in-jet efficiency, pT %.1f-%.1f:", pt_edges[i], pt_edges[i + 1]);
                    for (size_t k = 0; k < cf.size(); ++k) cout << Form("  dr %.1f-%.1f c = %.3f", InJet::kInJetDr[k], InJet::kInJetDr[k + 1], cf[k]);
                    cout << endl;
                }
            }
        } else {
            // No background removal (Gen, Matched, or flag disabled)
            corrected_via_mixing = (TH2*)corrected_via_mixing_events(normalized_signal_temp, normalized_mixing_temp);
        }

        rawsignal_after_mixing_correction.push_back(corrected_via_mixing);
        
        TH2 *seagull = (TH2*)do_seagull_correction(corrected_via_mixing, 0, 0, 0.5,i);
        signal_after_seagull_correction.push_back(seagull); 

        TH2 *raw_seagull_hist = (TH2*)seagull->Clone(Form("raw_seagull_%d", i));
        TH2 *signal_hist = (TH2*)corrected_via_mixing->Clone(Form("sig_hist_%d", i));
        raw_seagull_hist->Divide(signal_hist,raw_seagull_hist,1,1);
        raw_seagull.push_back(raw_seagull_hist);

        TH2 *backgroundSubtractedSeagull = (TH2*)background_subtraction(seagull, seagull, 4, true, gUEMin, gUEMax);
        signal_after_seagull_and_sideband_bkg_correction.push_back(backgroundSubtractedSeagull); 
        
        TH2 *backgroundSubtracted = background_subtraction(corrected_via_mixing, corrected_via_mixing, 4, true, gUEMin, gUEMax);
        signal_after_sideband_bkg_correction.push_back(backgroundSubtracted); 

        //Jet Shapes        
        TH1 *jets_shape_raw = (TH1*)get_jet_shape(normalized_signal,"rawsignal");
        jet_shape_histo_raw.push_back(jets_shape_raw);
        
        TH1 *jet_shape_signal_over_mixing = (TH1*)get_jet_shape(corrected_via_mixing,"mixingcorrected");
        jet_shape_histo_signal_over_mixing.push_back(jet_shape_signal_over_mixing);

        TH1 *jet_shape_seagull_corrected = (TH1*)get_jet_shape(seagull,"seagullcorrected");
        jet_shape_histo_seagull_corrected.push_back(jet_shape_seagull_corrected);

        TH1 *jet_shape_sideband_bkg_corrected = (TH1*)get_jet_shape(backgroundSubtracted,"bkgcorrected");
        jet_shape_histo_sideband_bkg_corrected.push_back(jet_shape_sideband_bkg_corrected);

        TH1 *jets_shape_fully_corrected = (TH1*)get_jet_shape(backgroundSubtractedSeagull,"fullycorrected");
        jet_shape_histo_fully_corrected.push_back(jets_shape_fully_corrected);

        //Integral vs pT
        int first_bin = jets_shape_raw->FindBin(dr_start);
        int last_bin =  jets_shape_raw->FindBin(dr_end);
        Integral_vs_pt_raw->SetBinContent(i+1,jets_shape_raw->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_raw->SetBinError(i+1,Integral_error);
        
        Integral_vs_pt_signal_over_mixing->SetBinContent(i+1,jet_shape_signal_over_mixing->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_signal_over_mixing->SetBinError(i+1,Integral_error);

        Integral_vs_pt_seagull_corrected->SetBinContent(i+1,jet_shape_seagull_corrected->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_seagull_corrected->SetBinError(i+1,Integral_error);

        Integral_vs_pt_sideband_bkg_corrected->SetBinContent(i+1,jet_shape_sideband_bkg_corrected->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_sideband_bkg_corrected->SetBinError(i+1,Integral_error);

        Integral_vs_pt_fully_corrected->SetBinContent(i+1,jets_shape_fully_corrected->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_fully_corrected->SetBinError(i+1,Integral_error);

        //Inner cone Integrals
        last_bin =  jets_shape_raw->FindBin(0.399);
        Integral_vs_pt_innercone_raw->SetBinContent(i+1,jets_shape_raw->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_innercone_raw->SetBinError(i+1,Integral_error);

        Integral_vs_pt_innercone_signal_over_mixing->SetBinContent(i+1,jet_shape_signal_over_mixing->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_innercone_signal_over_mixing->SetBinError(i+1,Integral_error);

        Integral_vs_pt_innercone_seagull_corrected->SetBinContent(i+1,jet_shape_seagull_corrected->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_innercone_seagull_corrected->SetBinError(i+1,Integral_error);

        Integral_vs_pt_innercone_sideband_bkg_corrected->SetBinContent(i+1,jet_shape_sideband_bkg_corrected->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_innercone_sideband_bkg_corrected->SetBinError(i+1,Integral_error);

        Integral_vs_pt_innercone_fully_corrected->SetBinContent(i+1,jets_shape_fully_corrected->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Integral_vs_pt_innercone_fully_corrected->SetBinError(i+1,Integral_error);

    }

    Integral_vs_pt_signal_over_mixing->Divide(Integral_vs_pt_bin_widths_scale);
    Integral_vs_pt_seagull_corrected->Divide(Integral_vs_pt_bin_widths_scale);
    Integral_vs_pt_sideband_bkg_corrected->Divide(Integral_vs_pt_bin_widths_scale);
    Integral_vs_pt_fully_corrected->Divide(Integral_vs_pt_bin_widths_scale);
    
    Integral_vs_pt_innercone_signal_over_mixing->Divide(Integral_vs_pt_bin_widths_scale);
    Integral_vs_pt_innercone_seagull_corrected->Divide(Integral_vs_pt_bin_widths_scale);
    Integral_vs_pt_innercone_sideband_bkg_corrected->Divide(Integral_vs_pt_bin_widths_scale);
    Integral_vs_pt_innercone_fully_corrected->Divide(Integral_vs_pt_bin_widths_scale);
    
    cout << "=========================================================" << endl;
    f2->cd();

    f2->mkdir("general_2D_plots");
    f2->cd("general_2D_plots");
    
    for (int i = 0; i < rawsignal.size(); i++){
        if(i!=0)cout << "Writting 2D histograms: " << (float)100*i/rawsignal.size() << " %" << endl;
        rawsignal[i]->Write(Form("rawsignal_2D_Deta_Dphi_c%i",i));
        mixing[i]->Write(Form("mixing_2D_Deta_Dphi_c%i",i));
        normsignal[i]->Write(Form("normrawsignal_2D_Deta_Dphi_c%i",i));
        normmixing[i]->Write(Form("normmixing_2D_Deta_Dphi_c%i",i));
        rawsignal_after_mixing_correction[i]->Write(Form("rawsignal_over_mixing_2D_Deta_Dphi_c%i",i));    
        signal_after_seagull_correction[i]->Write(Form("signal_seagull_corrected_2D_Deta_Dphi_c%i",i));
        raw_seagull[i]->Write(Form("raw_seagull_2D_Deta_Dphi_c%i",i));
        signal_after_sideband_bkg_correction[i]->Write(Form("signal_bkg_corrected_2D_Deta_Dphi_c%i",i));    
        signal_after_seagull_and_sideband_bkg_correction[i]->Write(Form("signal_seagull_and_bkg_corrected_2D_Deta_Dphi_c%i",i));    
    }

    cout << "=========================================================" << endl;

    f2->mkdir("general_1D_plots_integrated");
    f2->cd("general_1D_plots_integrated");

    for (int i = 0; i < rawsignal.size(); i++){
        
        if(i!=0)cout << "Writting 1D integrated projections: " << (float)100*i/rawsignal.size() << " %" << endl;
        
        //Delta Eta projections     
        TH1 *rawsignal_Deta_hist_int = ProjectRegionDeltaEta(rawsignal[i],fullranges.first,fullranges.second,"",dobinnorm); 
        rawsignal_Deta_hist_int->Write(Form("rawsignal_integrated_Deta_c%i",i));    
        TH1 *mixing_Deta_hist_int = ProjectRegionDeltaEta(mixing[i],fullranges.first,fullranges.second,"",dobinnorm); 
        mixing_Deta_hist_int->Write(Form("mixing_integrated_Deta_c%i",i));
        TH1 *normsignal_Deta_hist_int = ProjectRegionDeltaEta(normsignal[i],fullranges.first,fullranges.second,"",dobinnorm); 
        normsignal_Deta_hist_int->Write(Form("normrawsignal_integrated_Deta_c%i",i));
        TH1 *normmixing_Deta_hist_int = ProjectRegionDeltaEta(normmixing[i],fullranges.first,fullranges.second,"",dobinnorm); 
        normmixing_Deta_hist_int->Write(Form("normmixing_integrated_Deta_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Deta_hist_int = ProjectRegionDeltaEta(rawsignal_after_mixing_correction[i],fullranges.first,fullranges.second,"",dobinnorm); 
        rawsignal_after_mixing_correction_Deta_hist_int->Write(Form("rawsignal_over_mixing_integrated_Deta_c%i",i));
        TH1 *signal_seagull_corrected_correction_Deta_hist_int = ProjectRegionDeltaEta(signal_after_seagull_correction[i],fullranges.first,fullranges.second,"",dobinnorm); 
        signal_seagull_corrected_correction_Deta_hist_int->Write(Form("signal_seagull_corrected_integrated_Deta_c%i",i));

        TH1 *raw_seagull_Deta_hist_int = ProjectRegionDeltaEta(raw_seagull[i],fullranges.first,fullranges.second,"",dobinnorm); 
        raw_seagull_Deta_hist_int->Write(Form("raw_seagull_integrated_Deta_c%i",i));
        
        TH1 *signal_after_sideband_bkg_correction_Deta_hist_int = ProjectRegionDeltaEta(signal_after_sideband_bkg_correction[i],fullranges.first,fullranges.second,"",dobinnorm); 
        signal_after_sideband_bkg_correction_Deta_hist_int->Write(Form("signal_bkg_corrected_integrated_Deta_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Deta_hist_int = ProjectRegionDeltaEta(signal_after_seagull_and_sideband_bkg_correction[i],fullranges.first,fullranges.second,"",dobinnorm); 
        signal_after_seagull_and_sideband_bkg_correction_Deta_hist_int->Write(Form("signal_seagull_and_bkg_corrected_integrated_Deta_c%i",i));
        
        //Delta Phi projections
        TH1 *rawsignal_Dphi_hist_int = ProjectRegionDeltaPhi(rawsignal[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        rawsignal_Dphi_hist_int->Write(Form("rawsignal_integrated_Dphi_c%i",i));
        TH1 *mixing_Dphi_hist_int = ProjectRegionDeltaPhi(mixing[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        mixing_Dphi_hist_int->Write(Form("mixing_integrated_Dphi_c%i",i));
        TH1 *normsignal_Dphi_hist_int = ProjectRegionDeltaPhi(normsignal[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        normsignal_Dphi_hist_int->Write(Form("normrawsignal_integrated_Dphi_c%i",i));
        TH1 *normmixing_Dphi_hist_int = ProjectRegionDeltaPhi(normmixing[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        normmixing_Dphi_hist_int->Write(Form("normmixing_integrated_Dphi_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Dphi_hist_int = ProjectRegionDeltaPhi(rawsignal_after_mixing_correction[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        rawsignal_after_mixing_correction_Dphi_hist_int->Write(Form("rawsignal_over_mixing_integrated_Dphi_c%i",i));
        TH1 *signal_seagull_corrected_correction_Dphi_hist_int = ProjectRegionDeltaPhi(signal_after_seagull_correction[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        signal_seagull_corrected_correction_Dphi_hist_int->Write(Form("signal_seagull_corrected_integrated_Dphi_c%i",i));

        TH1 *raw_seagull_Dphi_hist_int = ProjectRegionDeltaPhi(raw_seagull[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        raw_seagull_Dphi_hist_int->Write(Form("raw_seagull_integrated_Dphi_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Dphi_hist_int = ProjectRegionDeltaPhi(signal_after_sideband_bkg_correction[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        signal_after_sideband_bkg_correction_Dphi_hist_int->Write(Form("signal_bkg_corrected_integrated_Dphi_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Dphi_hist_int = ProjectRegionDeltaPhi(signal_after_seagull_and_sideband_bkg_correction[i],fullranges.first,fullranges.second,"",dooneregion,dobinnorm);  
        signal_after_seagull_and_sideband_bkg_correction_Dphi_hist_int->Write(Form("signal_seagull_and_bkg_corrected_integrated_Dphi_c%i",i));

    }
    
cout << "=========================================================" << endl;

    f2->mkdir("general_1D_plots_diff_Deta");
    f2->cd("general_1D_plots_diff_Deta");

    for (int i = 0; i < rawsignal.size(); i++){
        
        if(i!=0)cout << "Writting 1D Delta Phi projections: " << (float)100*i/rawsignal.size() << " %" << endl;
        //Delta phi projections   
        TH1 *rawsignal_Dphi_hist_sig= ProjectRegionDeltaPhi(rawsignal[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        rawsignal_Dphi_hist_sig->Write(Form("rawsignal_signal_Dphi_c%i",i));    
        TH1 *mixing_Dphi_hist_sig= ProjectRegionDeltaPhi(mixing[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        mixing_Dphi_hist_sig->Write(Form("mixing_signal_Dphi_c%i",i));
        TH1 *normsignal_Dphi_hist_sig= ProjectRegionDeltaPhi(normsignal[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        normsignal_Dphi_hist_sig->Write(Form("normrawsignal_signal_Dphi_c%i",i));
        TH1 *normmixing_Dphi_hist_sig= ProjectRegionDeltaPhi(normmixing[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        normmixing_Dphi_hist_sig->Write(Form("normmixing_signal_Dphi_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Dphi_hist_sig= ProjectRegionDeltaPhi(rawsignal_after_mixing_correction[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        rawsignal_after_mixing_correction_Dphi_hist_sig->Write(Form("rawsignal_over_mixing_signal_Dphi_c%i",i));
        TH1 *signal_seagull_corrected_correction_Dphi_hist_sig= ProjectRegionDeltaPhi(signal_after_seagull_correction[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        signal_seagull_corrected_correction_Dphi_hist_sig->Write(Form("signal_seagull_corrected_signal_Dphi_c%i",i));

        TH1 *raw_seagull_Dphi_hist_sig= ProjectRegionDeltaPhi(raw_seagull[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        raw_seagull_Dphi_hist_sig->Write(Form("raw_seagull_signal_Dphi_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Dphi_hist_sig= ProjectRegionDeltaPhi(signal_after_sideband_bkg_correction[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        signal_after_sideband_bkg_correction_Dphi_hist_sig->Write(Form("signal_bkg_corrected_signal_Dphi_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Dphi_hist_sig= ProjectRegionDeltaPhi(signal_after_seagull_and_sideband_bkg_correction[i],signal_deta.first,signal_deta.second,"",dooneregion,dobinnorm);  
        signal_after_seagull_and_sideband_bkg_correction_Dphi_hist_sig->Write(Form("signal_seagull_and_bkg_corrected_signal_Dphi_c%i",i));  

        TH1 *rawsignal_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(rawsignal[i]); 
        rawsignal_Dphi_hist_bkg->Write(Form("rawsignal_background_Dphi_c%i",i));    
        TH1 *mixing_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(mixing[i]); 
        mixing_Dphi_hist_bkg->Write(Form("mixing_background_Dphi_c%i",i));
        TH1 *normsignal_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(normsignal[i]); 
        normsignal_Dphi_hist_bkg->Write(Form("normrawsignal_background_Dphi_c%i",i));
        TH1 *normmixing_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(normmixing[i]); 
        normmixing_Dphi_hist_bkg->Write(Form("normmixing_background_Dphi_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(rawsignal_after_mixing_correction[i]); 
        rawsignal_after_mixing_correction_Dphi_hist_bkg->Write(Form("rawsignal_over_mixing_background_Dphi_c%i",i));
        TH1 *signal_seagull_corrected_correction_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(signal_after_seagull_correction[i]); 
        signal_seagull_corrected_correction_Dphi_hist_bkg->Write(Form("signal_seagull_corrected_background_Dphi_c%i",i));

        TH1 *raw_seagull_Dphi_hist_bkg = ProjectBackgroundDeltaPhi(raw_seagull[i]); 
        raw_seagull_Dphi_hist_bkg->Write(Form("raw_seagull_background_Dphi_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(signal_after_sideband_bkg_correction[i]); 
        signal_after_sideband_bkg_correction_Dphi_hist_bkg->Write(Form("signal_bkg_corrected_background_Dphi_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Dphi_hist_bkg= ProjectBackgroundDeltaPhi(signal_after_seagull_and_sideband_bkg_correction[i]); 
        signal_after_seagull_and_sideband_bkg_correction_Dphi_hist_bkg->Write(Form("signal_seagull_and_bkg_corrected_background_Dphi_c%i",i));  
        
    }

    cout << "=========================================================" << endl;
    
    f2->mkdir("Olga_Dphi_Projections");
    f2->cd("Olga_Dphi_Projections");
    
    for (int i = 0; i < rawsignal.size(); i++){
        
        if(i!=0)cout << "Writting 1D Deta projections: " << (float)100*i/rawsignal.size() << " %" << endl;

        //Delta Phi projections

        //near side
        TH1 *rawsignal_Deta_hist_nearside_olga= ProjectRegionDeltaEta(rawsignal[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        rawsignal_Deta_hist_nearside_olga->Write(Form("signal_raw_nearside_olga_Deta_c%i",i));
        TH1 *mixing_Deta_hist_nearside_olga= ProjectRegionDeltaEta(mixing[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        mixing_Deta_hist_nearside_olga->Write(Form("mixing_nearside_olga_Deta_c%i",i));
        TH1 *normsignal_Deta_hist_nearside_olga= ProjectRegionDeltaEta(normsignal[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        normsignal_Deta_hist_nearside_olga->Write(Form("normsignal_raw_nearside_olga_Deta_c%i",i));
        TH1 *normmixing_Deta_hist_nearside_olga= ProjectRegionDeltaEta(normmixing[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        normmixing_Deta_hist_nearside_olga->Write(Form("normmixing_nearside_olga_Deta_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Deta_hist_nearside_olga= ProjectRegionDeltaEta(rawsignal_after_mixing_correction[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        rawsignal_after_mixing_correction_Deta_hist_nearside_olga->Write(Form("signal_signalovermixing_nearside_olga_Deta_c%i",i));
        TH1 *signal_seagull_corrected_correction_Deta_hist_nearside_olga= ProjectRegionDeltaEta(signal_after_seagull_correction[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        signal_seagull_corrected_correction_Deta_hist_nearside_olga->Write(Form("signal_seagull_corrected_nearside_olga_Deta_c%i",i));

        TH1 *raw_seagull_Deta_hist_nearside_olga = ProjectRegionDeltaEta(raw_seagull[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        raw_seagull_Deta_hist_nearside_olga->Write(Form("raw_seagull_nearside_olga_Deta_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Deta_hist_nearside_olga= ProjectRegionDeltaEta(signal_after_sideband_bkg_correction[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        signal_after_sideband_bkg_correction_Deta_hist_nearside_olga->Write(Form("signal_sideband_bkg_corrected_nearside_olga_Deta_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Deta_hist_nearside_olga= ProjectRegionDeltaEta(signal_after_seagull_and_sideband_bkg_correction[i],olga_nearside.first,olga_nearside.second,"",dobinnorm); 
        signal_after_seagull_and_sideband_bkg_correction_Deta_hist_nearside_olga->Write(Form("signal_fully_corrected_nearside_olga_Deta_c%i",i));

        //away side       
        TH1 *rawsignal_Deta_hist_awayside_olga= ProjectRegionDeltaEta(rawsignal[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        rawsignal_Deta_hist_awayside_olga->Write(Form("rawsignal_awayside_olga_Deta_c%i",i));
        TH1 *mixing_Deta_hist_awayside_olga= ProjectRegionDeltaEta(mixing[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        mixing_Deta_hist_awayside_olga->Write(Form("mixing_awayside_olga_Deta_c%i",i));
        TH1 *normsignal_Deta_hist_awayside_olga= ProjectRegionDeltaEta(normsignal[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        normsignal_Deta_hist_awayside_olga->Write(Form("normrawsignal_awayside_olga_Deta_c%i",i));
        TH1 *normmixing_Deta_hist_awayside_olga= ProjectRegionDeltaEta(normmixing[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        normmixing_Deta_hist_awayside_olga->Write(Form("normmixing_awayside_olga_Deta_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Deta_hist_awayside_olga= ProjectRegionDeltaEta(rawsignal_after_mixing_correction[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        rawsignal_after_mixing_correction_Deta_hist_awayside_olga->Write(Form("rawsignal_over_mixing_awayside_olga_Deta_c%i",i));
        TH1 *signal_seagull_corrected_correction_Deta_hist_awayside_olga= ProjectRegionDeltaEta(signal_after_seagull_correction[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        signal_seagull_corrected_correction_Deta_hist_awayside_olga->Write(Form("signal_seagull_corrected_awayside_olga_Deta_c%i",i));

        TH1 *raw_seagull_Deta_hist_awayside_olga = ProjectRegionDeltaEta(raw_seagull[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        raw_seagull_Deta_hist_awayside_olga->Write(Form("raw_seagull_awayside_olga_Deta_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Deta_hist_awayside_olga= ProjectRegionDeltaEta(signal_after_sideband_bkg_correction[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        signal_after_sideband_bkg_correction_Deta_hist_awayside_olga->Write(Form("signal_bkg_corrected_awayside_olga_Deta_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Deta_hist_awayside_olga= ProjectRegionDeltaEta(signal_after_seagull_and_sideband_bkg_correction[i],olga_awayside.first,olga_awayside.second,"",dobinnorm); 
        signal_after_seagull_and_sideband_bkg_correction_Deta_hist_awayside_olga->Write(Form("signal_seagull_and_bkg_corrected_awayside_olga_Deta_c%i",i));

        //between peaks   
        TH1 *rawsignal_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(rawsignal[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        rawsignal_Deta_hist_betweenpeaks_olga->Write(Form("rawsignal_betweenpeaks_olga_Deta_c%i",i));
        TH1 *mixing_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(mixing[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        mixing_Deta_hist_betweenpeaks_olga->Write(Form("mixing_betweenpeaks_olga_Deta_c%i",i));
        TH1 *normsignal_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(normsignal[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        normsignal_Deta_hist_betweenpeaks_olga->Write(Form("normrawsignal_betweenpeaks_olga_Deta_c%i",i));
        TH1 *normmixing_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(normmixing[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        normmixing_Deta_hist_betweenpeaks_olga->Write(Form("normmixing_betweenpeaks_olga_Deta_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(rawsignal_after_mixing_correction[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        rawsignal_after_mixing_correction_Deta_hist_betweenpeaks_olga->Write(Form("rawsignal_over_mixing_betweenpeaks_olga_Deta_c%i",i));
        TH1 *signal_seagull_corrected_correction_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(signal_after_seagull_correction[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        signal_seagull_corrected_correction_Deta_hist_betweenpeaks_olga->Write(Form("signal_seagull_corrected_betweenpeaks_olga_Deta_c%i",i));

        TH1 *raw_seagull_Deta_hist_betweenpeaks_olga = ProjectRegionDeltaEta(raw_seagull[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        raw_seagull_Deta_hist_betweenpeaks_olga->Write(Form("raw_seagull_betweenpeaks_olga_Deta_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(signal_after_sideband_bkg_correction[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        signal_after_sideband_bkg_correction_Deta_hist_betweenpeaks_olga->Write(Form("signal_bkg_corrected_betweenpeaks_olga_Deta_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Deta_hist_betweenpeaks_olga= ProjectRegionDeltaEta(signal_after_seagull_and_sideband_bkg_correction[i],olga_betweenpeaks.first,olga_betweenpeaks.second,"",dobinnorm); 
        signal_after_seagull_and_sideband_bkg_correction_Deta_hist_betweenpeaks_olga->Write(Form("signal_seagull_and_bkg_corrected_betweenpeaks_olga_Deta_c%i",i));
 
    }
    
    cout << "=========================================================" << endl;
    f2->mkdir("Jussi_Dphi_Projections");
    f2->cd("Jussi_Dphi_Projections");
    

    for (int i = 0; i < rawsignal.size(); i++){

        //near side
        TH1 *rawsignal_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(rawsignal[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        rawsignal_Deta_hist_nearside_jussi->Write(Form("rawsignal_nearside_jussi_Deta_c%i",i));
        TH1 *mixing_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(mixing[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        mixing_Deta_hist_nearside_jussi->Write(Form("mixing_nearside_jussi_Deta_c%i",i));
        TH1 *normsignal_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(normsignal[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        normsignal_Deta_hist_nearside_jussi->Write(Form("normrawsignal_nearside_jussi_Deta_c%i",i));
        TH1 *normmixing_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(normmixing[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        normmixing_Deta_hist_nearside_jussi->Write(Form("normmixing_nearside_jussi_Deta_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(rawsignal_after_mixing_correction[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        rawsignal_after_mixing_correction_Deta_hist_nearside_jussi->Write(Form("rawsignal_over_mixing_nearside_jussi_Deta_c%i",i));
        TH1 *signal_seagull_corrected_correction_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(signal_after_seagull_correction[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        signal_seagull_corrected_correction_Deta_hist_nearside_jussi->Write(Form("signal_seagull_corrected_nearside_jussi_Deta_c%i",i));

        TH1 *raw_seagull_Deta_hist_nearside_jussi = ProjectRegionDeltaEta(raw_seagull[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        raw_seagull_Deta_hist_nearside_jussi->Write(Form("raw_seagull_nearside_jussi_Deta_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(signal_after_sideband_bkg_correction[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        signal_after_sideband_bkg_correction_Deta_hist_nearside_jussi->Write(Form("signal_bkg_corrected_nearside_jussi_Deta_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Deta_hist_nearside_jussi= ProjectRegionDeltaEta(signal_after_seagull_and_sideband_bkg_correction[i],jussi_nearside.first,jussi_nearside.second,"",dobinnorm); 
        signal_after_seagull_and_sideband_bkg_correction_Deta_hist_nearside_jussi->Write(Form("signal_seagull_and_bkg_corrected_nearside_jussi_Deta_c%i",i));

        //away side       
        TH1 *rawsignal_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(rawsignal[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        rawsignal_Deta_hist_awayside_jussi->Write(Form("rawsignal_awayside_jussi_Deta_c%i",i));
        TH1 *mixing_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(mixing[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        mixing_Deta_hist_awayside_jussi->Write(Form("mixing_awayside_jussi_Deta_c%i",i));
        TH1 *normsignal_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(normsignal[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        normsignal_Deta_hist_awayside_jussi->Write(Form("normrawsignal_awayside_jussi_Deta_c%i",i));
        TH1 *normmixing_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(normmixing[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        normmixing_Deta_hist_awayside_jussi->Write(Form("normmixing_awayside_jussi_Deta_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(rawsignal_after_mixing_correction[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        rawsignal_after_mixing_correction_Deta_hist_awayside_jussi->Write(Form("rawsignal_over_mixing_awayside_jussi_Deta_c%i",i));
        TH1 *signal_seagull_corrected_correction_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(signal_after_seagull_correction[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        signal_seagull_corrected_correction_Deta_hist_awayside_jussi->Write(Form("signal_seagull_corrected_awayside_jussi_Deta_c%i",i));

        TH1 *raw_seagull_Deta_hist_awayside_jussi = ProjectRegionDeltaEta(raw_seagull[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        raw_seagull_Deta_hist_awayside_jussi->Write(Form("raw_seagull_awayside_jussi_Deta_c%i",i));

        TH1 *signal_after_sideband_bkg_correction_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(signal_after_sideband_bkg_correction[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        signal_after_sideband_bkg_correction_Deta_hist_awayside_jussi->Write(Form("signal_bkg_corrected_awayside_jussi_Deta_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Deta_hist_awayside_jussi= ProjectRegionDeltaEta(signal_after_seagull_and_sideband_bkg_correction[i],jussi_awayside.first,jussi_awayside.second,"",dobinnorm); 
        signal_after_seagull_and_sideband_bkg_correction_Deta_hist_awayside_jussi->Write(Form("signal_seagull_and_bkg_corrected_awayside_jussi_Deta_c%i",i));

        //between peaks   
        TH1 *rawsignal_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(rawsignal[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        rawsignal_Deta_hist_betweenpeaks_jussi->Write(Form("rawsignal_betweenpeaks_jussi_Deta_c%i",i));
        TH1 *mixing_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(mixing[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        mixing_Deta_hist_betweenpeaks_jussi->Write(Form("mixing_betweenpeaks_jussi_Deta_c%i",i));
        TH1 *normsignal_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(normsignal[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        normsignal_Deta_hist_betweenpeaks_jussi->Write(Form("normrawsignal_betweenpeaks_jussi_Deta_c%i",i));
        TH1 *normmixing_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(normmixing[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        normmixing_Deta_hist_betweenpeaks_jussi->Write(Form("normmixing_betweenpeaks_jussi_Deta_c%i",i));
        TH1 *rawsignal_after_mixing_correction_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(rawsignal_after_mixing_correction[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        rawsignal_after_mixing_correction_Deta_hist_betweenpeaks_jussi->Write(Form("rawsignal_over_mixing_betweenpeaks_jussi_Deta_c%i",i));
        TH1 *signal_seagull_corrected_correction_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(signal_after_seagull_correction[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        signal_seagull_corrected_correction_Deta_hist_betweenpeaks_jussi->Write(Form("signal_seagull_corrected_betweenpeaks_jussi_Deta_c%i",i));

        TH1 *raw_seagull_Deta_hist_betweenpeaks_jussi = ProjectRegionDeltaEta(raw_seagull[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        raw_seagull_Deta_hist_betweenpeaks_jussi->Write(Form("raw_seagull_betweenpeaks_jussi_Deta_c%i",i));
        
        TH1 *signal_after_sideband_bkg_correction_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(signal_after_sideband_bkg_correction[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        signal_after_sideband_bkg_correction_Deta_hist_betweenpeaks_jussi->Write(Form("signal_bkg_corrected_betweenpeaks_jussi_Deta_c%i",i));
        TH1 *signal_after_seagull_and_sideband_bkg_correction_Deta_hist_betweenpeaks_jussi= ProjectRegionDeltaEta(signal_after_seagull_and_sideband_bkg_correction[i],jussi_betweenpeaks.first,jussi_betweenpeaks.second,"",dobinnorm); 
        signal_after_seagull_and_sideband_bkg_correction_Deta_hist_betweenpeaks_jussi->Write(Form("signal_seagull_and_bkg_corrected_betweenpeaks_jussi_Deta_c%i",i));
 
    }

    
    if (!proc_sig_over_mix.empty()) {   // method B components (reco V0): signal/mixing before the sideband subtraction, sideband terms
        f2->mkdir("process_2D");
        f2->cd("process_2D");
        for (size_t i = 0; i < proc_sig_over_mix.size(); i++){
            proc_sig_over_mix[i]->Write(Form("sig_over_mix_c%zu",i));
            proc_norm_bkg[i]->Write(Form("normbkg_c%zu",i));
            proc_norm_mix_bkg[i]->Write(Form("normmixingbkg_c%zu",i));
            proc_bkg_over_mix[i]->Write(Form("bkg_over_mix_c%zu",i));
        }
    }
    f2->mkdir("jetshape");
    f2->cd("jetshape");
    for (int i = 0; i < rawsignal.size(); i++){
        jet_shape_histo_raw[i]->Write(Form("JetShape_raw_c%i",i));
        jet_shape_histo_signal_over_mixing[i]->Write(Form("JetShape_signalovermixing_c%i",i));
        jet_shape_histo_seagull_corrected[i]->Write(Form("JetShape_seagull_corrected_c%i",i));
        jet_shape_histo_sideband_bkg_corrected[i]->Write(Form("JetShape_sideband_corrected_c%i",i));
        jet_shape_histo_fully_corrected[i]->Write(Form("JetShape_fully_corrected_c%i",i));
    }

    f2->mkdir("Integral_vs_pt");
    f2->cd("Integral_vs_pt");
    Integral_vs_pt_raw->Write("Integral_vs_pt_raw");
    Integral_vs_pt_signal_over_mixing->Write("Integral_vs_pt_signalovermixing");
    Integral_vs_pt_seagull_corrected->Write("Integral_vs_pt_seagull_corrected");
    Integral_vs_pt_sideband_bkg_corrected->Write("Integral_vs_pt_sideband_corrected");
    Integral_vs_pt_fully_corrected->Write("Integral_vs_pt_fully_corrected");
    Integral_vs_pt_bin_widths_scale->Write();


    f2->mkdir("Integral_vs_pt_innercone");
    f2->cd("Integral_vs_pt_innercone");
    Integral_vs_pt_innercone_raw->Write("Integral_vs_pt_innercone_raw");
    Integral_vs_pt_innercone_signal_over_mixing->Write("Integral_vs_pt_innercone_signalovermixing");
    Integral_vs_pt_innercone_seagull_corrected->Write("Integral_vs_pt_innercone_seagull_corrected");
    Integral_vs_pt_innercone_sideband_bkg_corrected->Write("Integral_vs_pt_innercone_sideband_corrected");
    Integral_vs_pt_innercone_fully_corrected->Write("Integral_vs_pt_innercone_fully_corrected");

    f2->Close();
    
    cout << "---------------------------------------- DONE ----------------------------------------" << endl;

}