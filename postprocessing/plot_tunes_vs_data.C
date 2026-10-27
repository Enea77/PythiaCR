// plot_tunes_vs_data.C -- in-cone (Lambda+anti-Lambda)/2K0S and Delta r profiles: MB data vs PYTHIA tunes.
// Inputs are jetshape.cpp outputs (copied unchanged from v0-jet-analysis/postprocessing), same definitions as
// plot_incone.C: stage "sideband_corrected" (Delta eta UE subtracted, no seagull; nominal), Delta r < 0.4
// (Integral_vs_pt_innercone), ratio = 0.5 Lambda / K0S, 3-10 GeV = sum of the per-pT-bin yields (SumBins),
// Delta r = sum of the jet shapes of all V0 pT bins (1.5-10 GeV). Data: statistical errors only.
//
//   root -l -b -q 'plot_tunes_vs_data.C("<data Lam js>,<data K0s js>","<label>:<Lam js>:<K0s js>;<label>:...","plots/tunes_vs_data")'
// (see run_tunes_vs_data.sh)
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

namespace ptd {
const TString kStage = "sideband_corrected";
const int kCol[8] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet - 3, kGray + 2};
const int kMark[8] = {24, 25, 26, 32, 27, 28, 30, 46};

// as plot_incone.C (GetYield, SumBins, GetDR)
TH1D* GetYield(const TString& file, const TString& name) {
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { std::cout << "missing " << file << std::endl; return nullptr; }
    TH1* h = (TH1*)f->Get("Integral_vs_pt_innercone/Integral_vs_pt_innercone_" + kStage);
    TH1D* o = h ? (TH1D*)h->Clone(name) : nullptr;
    if (o) o->SetDirectory(0); else std::cout << "missing yield in " << file << std::endl;
    f->Close();
    return o;
}
TH1D* SumBins(TH1D* h, double lo, double hi, const TString& name) {
    if (!h) return nullptr;
    TH1D* o = new TH1D(name, "", 1, lo, hi); o->SetDirectory(0);
    double v = 0, e2 = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        if (h->GetBinCenter(b) < lo || h->GetBinCenter(b) > hi) continue;
        const double f = h->GetBinWidth(b) / h->GetBinWidth(1);
        v += h->GetBinContent(b) * f; e2 += std::pow(h->GetBinError(b) * f, 2);
    }
    o->SetBinContent(1, v); o->SetBinError(1, std::sqrt(e2));
    return o;
}
TH1D* GetDR(const TString& file, const TString& name) {
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) return nullptr;
    TH1D* s = nullptr;
    for (int i = 0;; ++i) {
        TH1* h = (TH1*)f->Get(Form("jetshape/JetShape_%s_c%d", kStage.Data(), i));
        if (!h) break;
        if (!s) { s = (TH1D*)h->Clone(name); s->SetDirectory(0); } else s->Add(h);
    }
    f->Close();
    return s;
}
TH1D* Ratio(TH1D* l, TH1D* k, const TString& name) {
    if (!l || !k) return nullptr;
    TH1D* r = (TH1D*)l->Clone(name); r->SetDirectory(0); r->Divide(k); r->Scale(0.5);
    return r;
}

struct Sample { TString label; TH1D *r = nullptr, *rI = nullptr, *dr = nullptr, *lam = nullptr, *k0s = nullptr; };

Sample Load(const TString& label, const TString& fl, const TString& fk, int i) {
    Sample s; s.label = label;
    s.lam = GetYield(fl, Form("lam_%d", i)); s.k0s = GetYield(fk, Form("k0s_%d", i));
    s.r = Ratio(s.lam, s.k0s, Form("r_%d", i));
    s.rI = Ratio(SumBins(s.lam, 3, 10, Form("lI_%d", i)), SumBins(s.k0s, 3, 10, Form("kI_%d", i)), Form("rI_%d", i));
    s.dr = Ratio(GetDR(fl, Form("dl_%d", i)), GetDR(fk, Form("dk_%d", i)), Form("dr_%d", i));
    return s;
}

// top: data points + MC; bottom: MC / data. Right-hand panel (optional): the 3-10 GeV bin.
void Draw(const Sample& data, const std::vector<Sample>& mc, TH1D* Sample::*what, TH1D* Sample::*whatI,
          const TString& xtitle, const TString& ytitle, const TString& out, double ymin, double ymax) {
    TCanvas c("c", "", whatI ? 1000 : 750, 750);
    const double xs = whatI ? 0.78 : 1.0;
    TPad top("top", "", 0, 0.32, xs, 1), bot("bot", "", 0, 0, xs, 0.32);
    top.SetBottomMargin(0.02); top.SetLeftMargin(0.14); bot.SetTopMargin(0.02); bot.SetBottomMargin(0.3); bot.SetLeftMargin(0.14);
    top.Draw(); bot.Draw();
    TH1D* d = data.*what;
    top.cd();
    TH1D* fr = (TH1D*)d->Clone("frame_" + out); fr->Reset(); fr->SetMinimum(ymin); fr->SetMaximum(ymin + (ymax - ymin) / 0.62);
    fr->GetYaxis()->SetTitle(ytitle); fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.2);
    fr->Draw("axis");
    TLegend leg(0.16, 0.74, 0.92, 0.93); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.03); leg.SetNColumns(2);
    for (size_t i = 0; i < mc.size(); ++i) {
        TH1D* h = mc[i].*what; if (!h) continue;
        h->SetLineColor(kCol[i % 8]); h->SetMarkerColor(kCol[i % 8]); h->SetMarkerStyle(kMark[i % 8]); h->SetLineWidth(2);
        h->Draw("e1 x0 same");
        leg.AddEntry(h, mc[i].label, "lp");
    }
    d->SetMarkerStyle(20); d->SetMarkerSize(1.3); d->SetLineColor(kBlack); d->SetMarkerColor(kBlack); d->SetLineWidth(2);
    d->Draw("e1 x0 same");
    leg.AddEntry(d, data.label, "lp");
    leg.Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
    tx.DrawLatex(0.17, 0.69, "pPb 8.16 TeV / PYTHIA 8 pp");
    tx.DrawLatex(0.17, 0.645, "jets p_{T} > 60 GeV, |#eta| < 2");
    tx.DrawLatex(0.17, 0.60, "UE subtracted");
    bot.cd();
    // ratio range from the lowest / highest MC/data point (including its error), 15% padding
    std::vector<TH1D*> rs;
    double lo = 1e30, hi = -1e30;
    for (size_t i = 0; i < mc.size(); ++i) {
        TH1D* h = mc[i].*what; if (!h) continue;
        TH1D* r = (TH1D*)h->Clone(Form("%s_over_data_%zu", out.Data(), i)); r->Divide(d);
        for (int b = 1; b <= r->GetNbinsX(); ++b) {
            if (r->GetBinContent(b) <= 0 || d->GetBinContent(b) <= 0) continue;
            lo = std::min(lo, r->GetBinContent(b) - r->GetBinError(b)); hi = std::max(hi, r->GetBinContent(b) + r->GetBinError(b));
        }
        rs.push_back(r);
    }
    if (lo > hi) { lo = 0.9; hi = 1.1; }
    lo = std::min(lo, 1.0); hi = std::max(hi, 1.0);
    const double pad = 0.15 * (hi - lo);
    TH1D* bf = (TH1D*)d->Clone("bframe_" + out); bf->Reset(); bf->SetMinimum(lo - pad); bf->SetMaximum(hi + pad);
    bf->GetXaxis()->SetTitle(xtitle); bf->GetYaxis()->SetTitle("MC / data"); bf->GetYaxis()->SetNdivisions(505);
    for (auto ax : {bf->GetXaxis(), bf->GetYaxis()}) { ax->SetLabelSize(0.1); ax->SetTitleSize(0.11); }
    bf->GetYaxis()->SetTitleOffset(0.55); bf->GetXaxis()->SetTitleOffset(1.1);
    bf->Draw("axis");
    TLine one(bf->GetXaxis()->GetXmin(), 1, bf->GetXaxis()->GetXmax(), 1); one.SetLineStyle(2); one.Draw();
    for (auto r : rs) r->Draw("e1 x0 same");
    if (whatI) {   // right-hand panel: 3-10 GeV, one column per sample
        c.cd();
        TPad rp("rp", "", xs, 0.32, 1, 1); rp.SetLeftMargin(0.02); rp.SetRightMargin(0.25); rp.SetBottomMargin(0.02); rp.Draw(); rp.cd();
        const int n = mc.size() + 1;
        TH1D* f2 = new TH1D("f2_" + out, "", n, 0, n); f2->SetMinimum(ymin); f2->SetMaximum(ymin + (ymax - ymin) / 0.62);   // same range as the main panel
        f2->GetXaxis()->SetLabelSize(0); f2->GetYaxis()->SetLabelSize(0); f2->Draw("axis");
        auto pt = [&](TH1D* h, double x, int col, int mk) {
            TGraphErrors* g = new TGraphErrors(1); g->SetPoint(0, x, h->GetBinContent(1)); g->SetPointError(0, 0, h->GetBinError(1));
            g->SetMarkerStyle(mk); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(1.3); g->SetLineWidth(2); g->Draw("p same");
        };
        pt(data.*whatI, 0.5, kBlack, 20);
        for (size_t i = 0; i < mc.size(); ++i) if (mc[i].*whatI) pt(mc[i].*whatI, i + 1.5, kCol[i % 8], kMark[i % 8]);
        TLatex t2; t2.SetNDC(); t2.SetTextSize(0.09); t2.DrawLatex(0.08, 0.92, "3-10 GeV");
        c.cd();
        c.SaveAs(out + ".png"); c.SaveAs(out + ".pdf");
        return;
    }
    c.SaveAs(out + ".png"); c.SaveAs(out + ".pdf");
}
} // namespace ptd

void plot_tunes_vs_data(TString dataFiles, TString mcList, TString outDir = "plots/tunes_vs_data") {
    using namespace ptd;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    TH1::AddDirectory(kFALSE);
    gSystem->mkdir(outDir, true);
    TObjArray* df = dataFiles.Tokenize(",");
    Sample data = Load("Data MB (stat. only)", ((TObjString*)df->At(0))->GetString(), ((TObjString*)df->At(1))->GetString(), 0);
    std::vector<Sample> mc;
    TObjArray* ml = mcList.Tokenize(";");
    for (int i = 0; i < ml->GetEntries(); ++i) {
        TObjArray* p = ((TObjString*)ml->At(i))->GetString().Tokenize(":");
        mc.push_back(Load(((TObjString*)p->At(0))->GetString(), ((TObjString*)p->At(1))->GetString(), ((TObjString*)p->At(2))->GetString(), i + 1));
    }
    if (!data.r) { std::cout << "no data ratio" << std::endl; return; }

    Draw(data, mc, &Sample::r, &Sample::rI, "p_{T}^{V0} (GeV)", "in-cone (#Lambda+#bar{#Lambda})/2K^{0}_{S}, #Deltar < 0.4",
         outDir + "/InCone_Ratio", 0.0, 0.42);
    Draw(data, mc, &Sample::dr, nullptr, "#Deltar", "(#Lambda+#bar{#Lambda})/2K^{0}_{S}, 1.5 < p_{T}^{V0} < 10 GeV",
         outDir + "/DeltaR_Ratio", 0.0, 0.75);
    Draw(data, mc, &Sample::lam, nullptr, "p_{T}^{V0} (GeV)", "in-cone #Lambda+#bar{#Lambda} per jet (jetshape.cpp units)",
         outDir + "/InCone_Yield_Lam", 0.0, 1.6 * data.lam->GetMaximum());
    Draw(data, mc, &Sample::k0s, nullptr, "p_{T}^{V0} (GeV)", "in-cone K^{0}_{S} per jet (jetshape.cpp units)",
         outDir + "/InCone_Yield_K0s", 0.0, 1.6 * data.k0s->GetMaximum());

    // table
    std::ofstream txt((outDir + "/values.txt").Data());
    auto row = [&](const Sample& s) {
        TString l = Form("%-34s 3-10 GeV: %.4f +- %.4f |", s.label.Data(), s.rI->GetBinContent(1), s.rI->GetBinError(1));
        for (int b = 1; b <= s.r->GetNbinsX(); ++b) l += Form(" %.3f", s.r->GetBinContent(b));
        if (data.rI && &s != &data) l += Form(" | 3-10 MC/data %.3f", s.rI->GetBinContent(1) / data.rI->GetBinContent(1));
        std::cout << l << std::endl; txt << l << "\n";
    };
    TString hdr = Form("%-34s %-28s| per pT bin:", "in-cone (L+AL)/2K0S, dr < 0.4", "");
    for (int b = 1; b <= data.r->GetNbinsX(); ++b) hdr += Form(" %g-%g", data.r->GetBinLowEdge(b), data.r->GetXaxis()->GetBinUpEdge(b));
    std::cout << hdr << std::endl; txt << hdr << "\n";
    row(data);
    for (const auto& s : mc) if (s.rI) row(s);
    txt.close();
}
