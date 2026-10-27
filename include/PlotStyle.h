// PlotStyle.h
// Common plotting style for the analysis (same look as histo_plots / histo_plots_with_ratios in
// JetAnalyzer/make_plots.h): DSL colours, open thick markers 71/72/75..., ticks on all four sides,
// centred axis titles, top and ratio pads of equal height sharing the x axis (no gap).
#ifndef PLOT_STYLE_H
#define PLOT_STYLE_H

#include <vector>
#include <string>
#include "TColor.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"

namespace PS {

inline int Color(int i) {   // same order as DSL_colors in make_plots.h: 0 black, 4 dark red, 5 dark green, 6 dark blue, 8 orange ...
    static const int rgb[][3] = {{0,0,0},{230,0,0},{0,190,0},{0,200,255},{150,0,0},{0,100,0},{0,0,200},{230,180,30},{250,100,0},{80,80,85},{85,40,80},{250,20,250}};
    const int n = sizeof(rgb) / sizeof(rgb[0]);
    return TColor::GetColor(rgb[i % n][0], rgb[i % n][1], rgb[i % n][2]);
}
// series i: 0 black 71, 1 dark red 72, 2 dark green 75, 3 dark blue 76, 4 orange 77 ...
inline int SeriesColor(int i)  { return i == 0 ? Color(0) : Color(i + 3); }
inline int SeriesMarker(int i) { return i == 0 ? 71 : (i == 1 ? 72 : 73 + i); }

inline void StyleSeries(TH1* h, int i, double msize = 1.1) {
    h->SetMarkerStyle(SeriesMarker(i)); h->SetMarkerSize(msize);
    h->SetMarkerColor(SeriesColor(i));  h->SetLineColor(SeriesColor(i));
    h->SetLineWidth(2); h->SetLineStyle(1);
}

inline void SetGlobal() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
}

// Splits the current pad into a top pad and a ratio pad of equal height with a shared x axis.
struct RatioPads { TPad* top; TPad* bot; };
inline RatioPads MakeRatioPads(const char* tag, double split = 0.5, double left = 0.2, double right = 0.03) {
    TVirtualPad* mother = gPad;
    TPad* t = new TPad(Form("top_%s", tag), "", 0.0, split, 1.0, 1.0);
    TPad* b = new TPad(Form("bot_%s", tag), "", 0.0, 0.0, 1.0, split);
    for (TPad* p : {t, b}) { p->SetTickx(1); p->SetTicky(1); p->SetLeftMargin(left); p->SetRightMargin(right); p->SetFillStyle(4000); }
    t->SetTopMargin(0.1); t->SetBottomMargin(0.0);
    b->SetTopMargin(0.0); b->SetBottomMargin(0.3);
    mother->cd(); t->Draw(); mother->cd(); b->Draw();
    return {t, b};
}

// axis styling for a top pad histogram (text sizes are fractions of the pad height)
inline void StyleTopAxes(TH1* h, const char* ytitle, double ts = 0.085) {
    h->GetYaxis()->SetTitle(ytitle); h->GetYaxis()->CenterTitle(1);
    h->GetYaxis()->SetTitleSize(ts * 1.1); h->GetYaxis()->SetLabelSize(ts * 0.85);
    h->GetYaxis()->SetTitleOffset(1.0);    h->GetYaxis()->SetLabelOffset(0.006);
    h->GetYaxis()->SetMaxDigits(3);        h->GetYaxis()->SetTickLength(0.025);
    h->GetXaxis()->SetLabelSize(0);        h->GetXaxis()->SetTitleSize(0);
    h->GetXaxis()->SetTickLength(0.035);
}
inline void StyleRatioAxes(TH1* h, const char* xtitle, const char* ytitle, double lo, double hi, double ts = 0.085) {
    h->GetYaxis()->SetRangeUser(lo, hi);
    h->GetXaxis()->SetTitle(xtitle); h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->CenterTitle(1);   h->GetYaxis()->CenterTitle(1);
    h->GetXaxis()->SetTitleSize(ts * 1.2); h->GetYaxis()->SetTitleSize(ts * 1.1);
    h->GetXaxis()->SetLabelSize(ts * 0.85); h->GetYaxis()->SetLabelSize(ts * 0.85);
    h->GetXaxis()->SetTitleOffset(1.05);    h->GetYaxis()->SetTitleOffset(1.0);
    h->GetXaxis()->SetLabelOffset(0.008);   h->GetYaxis()->SetLabelOffset(0.006);
    h->GetYaxis()->SetNdivisions(505);      h->GetXaxis()->SetTickLength(0.035);
    h->GetYaxis()->SetTickLength(0.025);    h->GetYaxis()->SetNoExponent(true);
}

inline void DrawLine(double x0, double x1, double y, int style = 2) {
    TLine* l = new TLine(x0, y, x1, y); l->SetLineStyle(style); l->SetLineColor(kGray + 2); l->Draw();
}

inline void Text(double x, double y, const char* s, double size = 0.08, int font = 42, int align = 12) {
    TLatex* t = new TLatex(); t->SetNDC(); t->SetTextFont(font); t->SetTextSize(size); t->SetTextAlign(align);
    t->DrawLatex(x, y, s);
}
inline void CMSLabel(double x, double y, const char* extra = "Simulation", double size = 0.08) {
    Text(x, y, Form("#font[61]{CMS} #font[52]{%s}", extra), size, 42, 12);
}

// creates a directory and its parents, one level at a time (gSystem->mkdir(path, kTRUE) fails for relative
// paths in ROOT on Windows)
inline void MkDirs(const TString& path) {
    TString p = path; p.ReplaceAll("\\", "/");
    const bool abs = p.BeginsWith("/");
    TObjArray* parts = p.Tokenize("/");
    TString cur;
    for (int i = 0; i < parts->GetEntries(); ++i) {
        const TString s = ((TObjString*)parts->At(i))->GetString();
        cur = (i == 0 && !abs) ? s : cur + "/" + s;
        if (gSystem->AccessPathName(cur)) gSystem->MakeDirectory(cur);
    }
    delete parts;
}

inline void Save(TCanvas* c, const TString& base) {   // .pdf and .png
    c->SaveAs(base + ".pdf");
    c->SaveAs(base + ".png");
}

} // namespace PS

#endif
