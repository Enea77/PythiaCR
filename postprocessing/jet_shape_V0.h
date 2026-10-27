#include "call_libraries.h"
#define PI ( 4.*atan(1.) )
#define ME_peak 0.119

//define projections
TH1 *ProjectRegionDeltaEta(const TH2 *deltaPhiDeltaEtaHistogram, const float minDeltaPhi, const float maxDeltaPhi, const char *newName, const bool donormbybin){

    int lowDeltaPhiBin = deltaPhiDeltaEtaHistogram->GetXaxis()->FindBin(minDeltaPhi + 0.0001);
    int highDeltaPhiBin = deltaPhiDeltaEtaHistogram->GetXaxis()->FindBin(maxDeltaPhi - 0.0001);

    int fnBinsProjectedOver = highDeltaPhiBin - lowDeltaPhiBin + 1;

    char histogramName[200];
    sprintf(histogramName, "%s%s", deltaPhiDeltaEtaHistogram->GetName(), newName);
    TH1 *projectedDeltaEta = deltaPhiDeltaEtaHistogram->ProjectionY(histogramName, lowDeltaPhiBin, highDeltaPhiBin);

    if(donormbybin) projectedDeltaEta->Scale(1.0 / fnBinsProjectedOver);
    projectedDeltaEta->Scale(deltaPhiDeltaEtaHistogram->GetXaxis()->GetBinWidth(1));

    return projectedDeltaEta;
}

TH1 *ProjectRegionDeltaPhi(const TH2 *deltaPhiDeltaEtaHistogram, const float minDeltaEta, const float maxDeltaEta, const char *newName, const bool oneSide, const bool donormbybin){

    bool oneRegion = (TMath::Abs(minDeltaEta) < 0.01 || TMath::Abs(minDeltaEta) > 5.);

    int lowNegativeDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-maxDeltaEta + 0.0001);
    int highNegativeDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-minDeltaEta - 0.0001);
    int lowPositiveDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(minDeltaEta + 0.0001);
    int highPositiveDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(maxDeltaEta - 0.0001);

    int fnBinsProjectedOver = highNegativeDeltaEtaBin - lowNegativeDeltaEtaBin + highPositiveDeltaEtaBin - lowPositiveDeltaEtaBin + 2;
    if (oneRegion)
        fnBinsProjectedOver = fnBinsProjectedOver - highNegativeDeltaEtaBin + lowPositiveDeltaEtaBin - 1;
    if (oneSide)
        fnBinsProjectedOver = highPositiveDeltaEtaBin - lowPositiveDeltaEtaBin - 1;

    char histogramName[200];
    sprintf(histogramName, "%s%s", deltaPhiDeltaEtaHistogram->GetName(), newName);
    TH1 *projectedDeltaPhi;
    if (oneRegion)
    {
        projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName, lowNegativeDeltaEtaBin, highPositiveDeltaEtaBin);
    }
    else if (oneSide)
    {
        projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName, lowPositiveDeltaEtaBin, highPositiveDeltaEtaBin);
    }
    else
    {
        projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName, lowNegativeDeltaEtaBin, highNegativeDeltaEtaBin);
        projectedDeltaPhi->Add(deltaPhiDeltaEtaHistogram->ProjectionX("dummyName", lowPositiveDeltaEtaBin, highPositiveDeltaEtaBin));
    }

    if(donormbybin) projectedDeltaPhi->Scale(1.0 / fnBinsProjectedOver);
    projectedDeltaPhi->Scale(deltaPhiDeltaEtaHistogram->GetYaxis()->GetBinWidth(1));

    return projectedDeltaPhi;
}

// IMPROVEMENT: Added default bounds parameters to synchronize with background_subtraction
TH1 *ProjectBackgroundDeltaPhi(TH2 *deltaPhiDeltaEtaHistogram, float bkgEtaMin = 1.5, float bkgEtaMax = 2.5){
    bool fOneBackgroundRegion = false;
    bool donormbybin = true;
    return ProjectRegionDeltaPhi(deltaPhiDeltaEtaHistogram, bkgEtaMin, bkgEtaMax, "BackgroundDeltaPhi", fOneBackgroundRegion, donormbybin);
}

//raw signal normalization
TH2 *raw_signal_normalization(TH2 *signal, TH1* ptbin, double fMaximumDeltaEta){
  float njets = ptbin->Integral("width");
  cout<<"Number of jets for normalization: "<<njets<<endl;

  char newName[100];
  sprintf(newName,"%sCorrected",signal->GetName());
  TH2D* correctedHistogram = (TH2D*) signal->Clone(newName);

  int minimumDeltaEtaBin = correctedHistogram->GetYaxis()->FindBin(-fMaximumDeltaEta-0.001);
  int maximumDeltaEtaBin = correctedHistogram->GetYaxis()->FindBin(fMaximumDeltaEta+0.001);
  for(int iDeltaPhi = 1; iDeltaPhi <= correctedHistogram->GetNbinsX(); iDeltaPhi++){
    for(int iDeltaEta = 1; iDeltaEta <= minimumDeltaEtaBin; iDeltaEta++){
      correctedHistogram->SetBinContent(iDeltaPhi,iDeltaEta,0);
      correctedHistogram->SetBinError(iDeltaPhi,iDeltaEta,0);
    }
    for(int iDeltaEta = maximumDeltaEtaBin; iDeltaEta <= correctedHistogram->GetNbinsY(); iDeltaEta++){
      correctedHistogram->SetBinContent(iDeltaPhi,iDeltaEta,0);
      correctedHistogram->SetBinError(iDeltaPhi,iDeltaEta,0);
    }
  }
  
    correctedHistogram->Scale(1./njets); 
    return correctedHistogram;
} 

//mixing scaling for histogram
double mixing_scale(TH2 *mixing_histogram, const bool find_peak){

    TH1 *hDeltaEtaMixed;
    int nBinsProjectedOver;
    double binSum;
    double maxSum = 0;
    int peakFinderSideWidth = 16; 
    float fMixedEventFitRegionLow;
    float fMixedEventFitRegionHigh;
    
    if (find_peak)
    {
        for (int iDeltaPhi = 1; iDeltaPhi < mixing_histogram->GetNbinsX() + 2 - peakFinderSideWidth; iDeltaPhi++)
        {
            for (int iDeltaEta = 1; iDeltaEta < mixing_histogram->GetNbinsY(); iDeltaEta++)
            {
                binSum = 0;
                for (int iPeakFinderPhi = 0; iPeakFinderPhi < peakFinderSideWidth; iPeakFinderPhi++)
                {
                    for (int iPeakFinderEta = 0; iPeakFinderEta < peakFinderSideWidth; iPeakFinderEta++)
                    {
                        binSum += mixing_histogram->GetBinContent(iDeltaPhi + iPeakFinderPhi, iDeltaEta + iPeakFinderEta);
                    }
                }
                if (binSum > maxSum) maxSum = binSum;
            }
        }
        binSum = maxSum;
        nBinsProjectedOver = peakFinderSideWidth * peakFinderSideWidth;
    }
    else
    {
        hDeltaEtaMixed = mixing_histogram->ProjectionY("MixedDeltaEtaProjection", 1, mixing_histogram->GetNbinsX());
        nBinsProjectedOver = mixing_histogram->GetNbinsX();
        hDeltaEtaMixed->Fit("pol0", "Q0", "", fMixedEventFitRegionLow, fMixedEventFitRegionHigh);
        binSum = hDeltaEtaMixed->GetFunction("pol0")->GetParameter(0);
    }

  return (binSum / nBinsProjectedOver);
}

//return mixing normalized histogram
//return mixing normalized histogram
TH2 *mixing_normalization(TH2 *leadingMixedEventHistogram, const bool find_peak, const bool smoothandavoidpeaks, double fMixedEventFitRegion){
  char newName[100];
  double normalizationScale = mixing_scale(leadingMixedEventHistogram,find_peak);
  
  // Normalize the mixed event histogram
  sprintf(newName,"%sNormalized",leadingMixedEventHistogram->GetName());
  TH2D* fNormalizedMixedEventHistogram = (TH2D*)leadingMixedEventHistogram->Clone(newName);
  fNormalizedMixedEventHistogram->Scale(1.0/normalizationScale);
  
  if(smoothandavoidpeaks){
    TH1D *etaProjection = leadingMixedEventHistogram->ProjectionY("_eta");
    int maxBinEta = etaProjection->GetMaximumBin();
    float maxPointEta = etaProjection->GetBinCenter(maxBinEta);
    int lowFitBin = etaProjection->FindBin(maxPointEta-fMixedEventFitRegion);
    int highFitBin = etaProjection->FindBin(maxPointEta+fMixedEventFitRegion);
    double mean = 0;
    for(int i = lowFitBin; i <= highFitBin; i++) mean += etaProjection->GetBinContent(i);
    mean /= ((highFitBin-lowFitBin+1)*leadingMixedEventHistogram->GetNbinsX());
    etaProjection->Scale(1.0/mean);
    
    // Fill the smoothed projection, but protect the peak
    for(int iEta = 1; iEta <= leadingMixedEventHistogram->GetNbinsY(); iEta++){
      double binEta = leadingMixedEventHistogram->GetYaxis()->GetBinCenter(iEta);
      
      for(int iPhi = 1; iPhi <= leadingMixedEventHistogram->GetNbinsX(); iPhi++){
        if(TMath::Abs(binEta) < 0.2){ 
          // Protect the central peak region to perfectly trust the raw yield
          fNormalizedMixedEventHistogram->SetBinContent(iPhi, iEta, 1.0);
          fNormalizedMixedEventHistogram->SetBinError(iPhi, iEta, 0.0);
        } else {
          // Use the smoothed eta projection for the rest of the acceptance
          fNormalizedMixedEventHistogram->SetBinContent(iPhi, iEta, etaProjection->GetBinContent(iEta)/leadingMixedEventHistogram->GetNbinsX());
          fNormalizedMixedEventHistogram->SetBinError(iPhi, iEta, etaProjection->GetBinError(iEta)/sqrt(leadingMixedEventHistogram->GetNbinsX()));
        }
      }
    }
  } else {
    // Even if we are not smoothing, we still protect the peak from statistical fluctuations
    for(int iEta = 1; iEta <= leadingMixedEventHistogram->GetNbinsY(); iEta++){
      double binEta = leadingMixedEventHistogram->GetYaxis()->GetBinCenter(iEta);
      
      if(TMath::Abs(binEta) < 0.2){ 
        for(int iPhi = 1; iPhi <= leadingMixedEventHistogram->GetNbinsX(); iPhi++){
          fNormalizedMixedEventHistogram->SetBinContent(iPhi, iEta, 1.0);
          fNormalizedMixedEventHistogram->SetBinError(iPhi, iEta, 0.0);
        }
      }
    }
  }
  
  TH2D *noErrorMixedEventHistogram = (TH2D*) fNormalizedMixedEventHistogram->Clone("mixingButNoErrors");
  
  return noErrorMixedEventHistogram;
}

//take the ratio raw signal over mixing for first correction
TH2 *corrected_via_mixing_events(TH2 *signal, TH2 *mixing){
    signal->Divide(signal, mixing);
    return signal;
}

//=================== Seagull correction =================== 
double seagullPoly1(double *x, double *par){
    if (x[0] < par[3] && x[0] > (-1 * par[3])) return par[0];
    if (x[0] <= (-1 * par[3])) return par[0] + par[1] * (x[0] + par[3]);
    return par[0] + par[2] * (x[0] - par[3]);
}

double seagullPoly2(double *x, double *par){
    if (x[0] < par[3] && x[0] > (-1 * par[3])) return par[0];
    return par[0] + par[1] * x[0] * x[0] + par[2] * x[0];
}

double seagullExp(double *x, double *par){
    if (x[0] >= 0) return par[0] + par[1] * TMath::Exp(x[0] * par[2]);
    return par[0] + par[3] * TMath::Exp(x[0] * par[4]);
}

double seagullPolyExp(double *x, double *par){
    if (x[0] >= 0) return par[0] + par[1] * TMath::Exp(x[0] * par[2]) + par[3] * x[0] * x[0];
    return par[0] + par[4] * TMath::Exp(x[0] * par[5]) + par[6] * x[0] * x[0];
}

TH2 *do_seagull_correction(const TH2 *mixedEventCorrectedHistogram, const int normalizationMethod, const int vetoFlag, const double middleArea, int plotnumber){

    float fMinBackgroundDeltaPhi = TMath::Pi()/2-0.1;
    float fMaxBackgroundDeltaPhi = TMath::Pi()/2+0.1;
    int minDeltaPhiBin = mixedEventCorrectedHistogram->GetXaxis()->FindBin(fMinBackgroundDeltaPhi + 0.0001);
    int maxDeltaPhiBin = mixedEventCorrectedHistogram->GetXaxis()->FindBin(fMaxBackgroundDeltaPhi + 0.0001);
    TH1 *fBackgroundEtaProjection = mixedEventCorrectedHistogram->ProjectionY("SeagullEtaProjection", minDeltaPhiBin, maxDeltaPhiBin);

    float fSeagullRebin = 4.0;
    fBackgroundEtaProjection->Scale(1.0 / (maxDeltaPhiBin - minDeltaPhiBin + 1));
    fBackgroundEtaProjection->Scale(mixedEventCorrectedHistogram->GetXaxis()->GetBinWidth(1));
    fBackgroundEtaProjection->Rebin(fSeagullRebin);
    fBackgroundEtaProjection->Scale(1.0 / fSeagullRebin);

    TH2 *seagullCorrectedHistogram = (TH2 *)mixedEventCorrectedHistogram->Clone();

    float fSeagullChi2Limit = 1.15;
    if (vetoFlag != 1)
    {
        // Reverted to strict bounds to avoid empty acceptance regions
        fBackgroundEtaProjection->Fit("pol0", "Q0", "", -2.5, 2.5);
        TF1 *fSeagullFit = fBackgroundEtaProjection->GetFunction("pol0");
        double chi2PerNdf = fSeagullFit->GetChisquare() / fSeagullFit->GetNDF();
        if (chi2PerNdf < fSeagullChi2Limit || vetoFlag == 2)
            return seagullCorrectedHistogram;

        fBackgroundEtaProjection->RecursiveRemove(fSeagullFit);
    }

    if (normalizationMethod == 1 || normalizationMethod == 3 || normalizationMethod == 4)
    {
        double negativeContent, negativeError;
        double positiveContent, positiveError;
        double averageContent, averageError;
        int nDeltaEtaBins = fBackgroundEtaProjection->GetNbinsX();
        for (int iBin = 1; iBin <= nDeltaEtaBins / 2; iBin++)
        {
            negativeContent = fBackgroundEtaProjection->GetBinContent(iBin);
            negativeError = fBackgroundEtaProjection->GetBinError(iBin);
            positiveContent = fBackgroundEtaProjection->GetBinContent(nDeltaEtaBins - iBin + 1);
            positiveError = fBackgroundEtaProjection->GetBinError(nDeltaEtaBins - iBin + 1);
            averageContent = (positiveContent + negativeContent) / 2.0;
            averageError = (positiveError + negativeError) / 2.0;

            fBackgroundEtaProjection->SetBinContent(iBin, averageContent);
            fBackgroundEtaProjection->SetBinError(iBin, averageError);
            fBackgroundEtaProjection->SetBinContent(nDeltaEtaBins - iBin + 1, averageContent);
            fBackgroundEtaProjection->SetBinError(nDeltaEtaBins - iBin + 1, averageError);
        }
    }

    int nFitParameters[] = {4, 3, 4, 3, 4, 5, 7};

    // Reverted bounds for polynomial fits
    TF1 *fSeagullFit = new TF1("seagullFit", seagullPoly2, -2.5, 2.5, 4);
    double initialLevel = fBackgroundEtaProjection->GetBinContent(fBackgroundEtaProjection->FindBin(0));
    fSeagullFit->SetParameters(initialLevel, 0, 0);
    fSeagullFit->FixParameter(3, middleArea);
    fBackgroundEtaProjection->Fit(fSeagullFit, "Q0", "", -2.5, 2.5);

    double backgroundLevel = fSeagullFit->GetParameter(0);
    double maxFitRange = 2.5;
    double minFitRange = -2.5;
    if (normalizationMethod == 3) maxFitRange = 2.0; 
    if (normalizationMethod == 1 || normalizationMethod == 3 || normalizationMethod == 4) minFitRange = 0.0;

    if (normalizationMethod == 1 || normalizationMethod == 3 || normalizationMethod == 5)
    {
        fBackgroundEtaProjection->RecursiveRemove(fSeagullFit);
        fSeagullFit = new TF1("seagullFitExp", seagullExp, -2.5, 2.5, nFitParameters[normalizationMethod]);
        initialLevel = fBackgroundEtaProjection->GetBinContent(fBackgroundEtaProjection->FindBin(0));
        fSeagullFit->SetParameters(initialLevel, -1, -1, -1, 1); 
        fBackgroundEtaProjection->Fit(fSeagullFit, "Q0", "", minFitRange, maxFitRange);
        backgroundLevel = fSeagullFit->GetParameter(0);
    }

    if (normalizationMethod == 4 || normalizationMethod == 6)
    {
        fBackgroundEtaProjection->RecursiveRemove(fSeagullFit);
        fSeagullFit = new TF1("seagullFitPolyExp", seagullPolyExp, -2.5, 2.5, nFitParameters[normalizationMethod]);
        initialLevel = fBackgroundEtaProjection->GetBinContent(fBackgroundEtaProjection->FindBin(0));
        fSeagullFit->SetParameters(initialLevel, -1, -1, 0, -1, 1, 0); 
        fBackgroundEtaProjection->Fit(fSeagullFit, "Q0", "", minFitRange, maxFitRange);
        backgroundLevel = fSeagullFit->GetParameter(0);
    }

    if (normalizationMethod == 2)
    {
        fBackgroundEtaProjection->RecursiveRemove(fSeagullFit);
        fSeagullFit = new TF1("seagullFitPoly1", seagullPoly1, -2.5, 2.5, 4);
        fSeagullFit->SetParameters(initialLevel, 0, 0);
        fSeagullFit->FixParameter(3, middleArea);
        fBackgroundEtaProjection->Fit(fSeagullFit, "Q0", "", -2.5, 2.5);
    }

    double binEta;
    double seagullCorrection;
    double binContent;
    double binError;

    for (int iEta = 1; iEta <= seagullCorrectedHistogram->GetNbinsY(); iEta++)
    {
        binEta = seagullCorrectedHistogram->GetYaxis()->GetBinCenter(iEta);
        if (normalizationMethod == 1 || normalizationMethod == 3 || normalizationMethod == 4)
            binEta = TMath::Abs(binEta);
        seagullCorrection = backgroundLevel / fSeagullFit->Eval(binEta);

        for (int iPhi = 1; iPhi <= seagullCorrectedHistogram->GetNbinsX(); iPhi++)
        {
            binContent = seagullCorrectedHistogram->GetBinContent(iPhi, iEta);
            binError = seagullCorrectedHistogram->GetBinError(iPhi, iEta);
            seagullCorrectedHistogram->SetBinContent(iPhi, iEta, binContent * seagullCorrection);
            seagullCorrectedHistogram->SetBinError(iPhi, iEta, binError * seagullCorrection);
        }
    }

    // (the fit plots Seagull_<i>.png/.pdf written here were removed: overwritten by every call, see plot_process.C Seagull_Check)

    return seagullCorrectedHistogram;
}

//Apply sideband background subtraction
// IMPROVEMENT: Added default parameters to link with ProjectBackgroundDeltaPhi 
TH2 *background_subtraction(TH2 *leadingHistogramWithBackground, TH2 *subleadingHistogramWithBackground, double maxDeltaEta, bool isInclusive, float bkgEtaMin = 1.5, float bkgEtaMax = 2.5){

    TH1 *backgroundDeltaPhiLeading = ProjectBackgroundDeltaPhi(leadingHistogramWithBackground, bkgEtaMin, bkgEtaMax);
    TH1 *backgroundDeltaPhiSubleading = ProjectBackgroundDeltaPhi(subleadingHistogramWithBackground, bkgEtaMin, bkgEtaMax);

    char histogramName[200];  
    double deltaPhiValueNear; 
    double deltaPhiErrorNear; 
    double deltaPhiValueAway; 
    double deltaPhiErrorAway; 

    double leadingOverlap = 0;
    double leadingOverlapError = 0;
    double subleadingOverlap = 0;
    double subleadingOverlapError = 0;

    sprintf(histogramName, "%sBackground", leadingHistogramWithBackground->GetName());
    TH2 *fBackgroundDistribution = (TH2 *)leadingHistogramWithBackground->Clone(histogramName);
    int nDeltaPhiBins = fBackgroundDistribution->GetNbinsX(); 

    sprintf(histogramName, "%sBackgroundOverlap", leadingHistogramWithBackground->GetName());
    TH2 *fBackgroundOverlap = (TH2 *)leadingHistogramWithBackground->Clone(histogramName);

    double binWidthDeltaEta = leadingHistogramWithBackground->GetYaxis()->GetBinWidth(1);
    int offset = isInclusive ? nDeltaPhiBins / 2 : 0; 

    int minFilledDeltaEtaBin = leadingHistogramWithBackground->GetYaxis()->FindBin(-maxDeltaEta + 0.001);
    int maxFilledDeltaEtaBin = leadingHistogramWithBackground->GetYaxis()->FindBin(maxDeltaEta - 0.001);

    int fnOverlapBins = 3;
    int fnBinsProjectedOver = 0;
    for (int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins / 2; iDeltaPhi++){
        deltaPhiValueNear = backgroundDeltaPhiLeading->GetBinContent(iDeltaPhi);
        deltaPhiErrorNear = backgroundDeltaPhiLeading->GetBinError(iDeltaPhi);
        deltaPhiValueAway = backgroundDeltaPhiSubleading->GetBinContent(iDeltaPhi + offset);
        deltaPhiErrorAway = backgroundDeltaPhiSubleading->GetBinError(iDeltaPhi + offset);

        if (iDeltaPhi <= fnOverlapBins){
            subleadingOverlap += deltaPhiValueAway;
            subleadingOverlapError += deltaPhiErrorAway;
        }
        if (iDeltaPhi > (nDeltaPhiBins / 2 - fnOverlapBins)){
            leadingOverlap += deltaPhiValueNear;
            leadingOverlapError += deltaPhiErrorNear;
        }

        for (int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
            if (iDeltaEta >= minFilledDeltaEtaBin && iDeltaEta <= maxFilledDeltaEtaBin){
                fBackgroundDistribution->SetBinContent(iDeltaPhi, iDeltaEta, deltaPhiValueNear / binWidthDeltaEta);
                fBackgroundDistribution->SetBinError(iDeltaPhi, iDeltaEta, deltaPhiErrorNear * TMath::Sqrt(fnBinsProjectedOver) / binWidthDeltaEta);
                fBackgroundDistribution->SetBinContent(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, deltaPhiValueAway / binWidthDeltaEta);
                fBackgroundDistribution->SetBinError(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, deltaPhiErrorAway * TMath::Sqrt(fnBinsProjectedOver) / binWidthDeltaEta);
            }

            fBackgroundOverlap->SetBinContent(iDeltaPhi, iDeltaEta, 0);
            fBackgroundOverlap->SetBinError(iDeltaPhi, iDeltaEta, 0);
            fBackgroundOverlap->SetBinContent(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, 0);
            fBackgroundOverlap->SetBinError(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, 0);
        }
    }

    for (int iDeltaPhi = nDeltaPhiBins / 2 + 1; iDeltaPhi <= nDeltaPhiBins / 2 + fnOverlapBins; iDeltaPhi++){
        deltaPhiValueNear = backgroundDeltaPhiLeading->GetBinContent(iDeltaPhi);
        deltaPhiErrorNear = backgroundDeltaPhiLeading->GetBinError(iDeltaPhi);

        deltaPhiValueAway = backgroundDeltaPhiSubleading->GetBinContent(iDeltaPhi + nDeltaPhiBins / 2 - fnOverlapBins);
        deltaPhiErrorAway = backgroundDeltaPhiSubleading->GetBinError(iDeltaPhi + nDeltaPhiBins / 2 - fnOverlapBins);

        leadingOverlap += deltaPhiValueNear;
        leadingOverlapError += deltaPhiErrorNear;
        subleadingOverlap += deltaPhiValueAway;
        subleadingOverlapError += deltaPhiErrorAway;

        for (int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
            if (iDeltaEta >= minFilledDeltaEtaBin && iDeltaEta <= maxFilledDeltaEtaBin){
                fBackgroundOverlap->SetBinContent(iDeltaPhi, iDeltaEta, deltaPhiValueNear / binWidthDeltaEta);
                fBackgroundOverlap->SetBinError(iDeltaPhi, iDeltaEta, deltaPhiErrorNear * TMath::Sqrt(fnBinsProjectedOver) / binWidthDeltaEta);
                fBackgroundOverlap->SetBinContent(iDeltaPhi - fnOverlapBins, iDeltaEta, deltaPhiValueAway / binWidthDeltaEta);
                fBackgroundOverlap->SetBinError(iDeltaPhi - fnOverlapBins, iDeltaEta, deltaPhiErrorAway * TMath::Sqrt(fnBinsProjectedOver) / binWidthDeltaEta);
            }
        } 
    } 

    double leadingScalingFactor = 1;
    double subleadingScalingFactor = 1;
    bool fAdjustBackground = false;
    if (fAdjustBackground && !isInclusive){
        double leadingToSubleadingRatio = leadingOverlap / subleadingOverlap;
        double subleadingToLeadingRatio = subleadingOverlap / leadingOverlap;
        double leadingToSubleadingRatioError = TMath::Sqrt(TMath::Power(leadingOverlapError / subleadingOverlap, 2) + TMath::Power(leadingOverlap * subleadingOverlapError / TMath::Power(subleadingOverlap, 2), 2));
        double scaledLeadingContent, scaledLeadingError;
        double scaledSubleadingContent, scaledSubleadingError;
        if (leadingToSubleadingRatioError < 1){
            if (leadingToSubleadingRatio > 1){
                leadingScalingFactor = 1 + ((leadingToSubleadingRatio - 1) / 2);
                subleadingScalingFactor = 1 - ((1 - subleadingToLeadingRatio) / 2);
            }else{
                leadingScalingFactor = 1 - ((1 - leadingToSubleadingRatio) / 2);
                subleadingScalingFactor = 1 + ((subleadingToLeadingRatio - 1) / 2);
            }

            for (int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins / 2; iDeltaPhi++){
                for (int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
                    if (iDeltaEta >= minFilledDeltaEtaBin && iDeltaEta <= maxFilledDeltaEtaBin){
                        scaledLeadingContent = fBackgroundDistribution->GetBinContent(iDeltaPhi, iDeltaEta) / leadingScalingFactor;
                        scaledLeadingError = fBackgroundDistribution->GetBinError(iDeltaPhi, iDeltaEta) / leadingScalingFactor;
                        scaledSubleadingContent = fBackgroundDistribution->GetBinContent(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta) / subleadingScalingFactor;
                        scaledSubleadingError = fBackgroundDistribution->GetBinError(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta) / subleadingScalingFactor;
                        fBackgroundDistribution->SetBinContent(iDeltaPhi, iDeltaEta, scaledLeadingContent);
                        fBackgroundDistribution->SetBinError(iDeltaPhi, iDeltaEta, scaledLeadingError);
                        fBackgroundDistribution->SetBinContent(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, scaledSubleadingContent);
                        fBackgroundDistribution->SetBinError(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, scaledSubleadingError);

                        scaledLeadingContent = fBackgroundOverlap->GetBinContent(iDeltaPhi, iDeltaEta) / subleadingScalingFactor;
                        scaledLeadingError = fBackgroundOverlap->GetBinError(iDeltaPhi, iDeltaEta) / subleadingScalingFactor;
                        scaledSubleadingContent = fBackgroundOverlap->GetBinContent(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta) / leadingScalingFactor;
                        scaledSubleadingError = fBackgroundOverlap->GetBinError(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta) / leadingScalingFactor;
                        fBackgroundOverlap->SetBinContent(iDeltaPhi, iDeltaEta, scaledLeadingContent);
                        fBackgroundOverlap->SetBinError(iDeltaPhi, iDeltaEta, scaledLeadingError);
                        fBackgroundOverlap->SetBinContent(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, scaledSubleadingContent);
                        fBackgroundOverlap->SetBinError(iDeltaPhi + nDeltaPhiBins / 2, iDeltaEta, scaledSubleadingError);
                    } 
                } 
            } 
        } 
    } 

    sprintf(histogramName, "%sBackgroundSubtracted", leadingHistogramWithBackground->GetName());
    TH2 *backgroundSubtractedHistogram = (TH2 *)leadingHistogramWithBackground->Clone(histogramName);
    backgroundSubtractedHistogram->Add(fBackgroundDistribution, -1);

    // IMPROVEMENT: Automatically utilize passed synchronized parameters
    double fBackgroundErrorScalingFactor = TMath::Sqrt(maxDeltaEta / (bkgEtaMax - bkgEtaMin));

    return backgroundSubtractedHistogram;
}

//Get jet shapes
TH1 *get_jet_shape(TH2 *backgroundSubtractedHistogram, TString extraword){
    
    //const int nRBins = 16;                                                                                                
    //double RBins[nRBins + 1] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5}; 

    //const int nRBins = 11;    double RBins[nRBins + 1] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.25, 1.5}; 

    const int nRBins = 6;    double RBins[nRBins + 1] = {0.0, 0.1, 0.2, 0.4, 0.6, 1.0, 1.5}; 

    char histogramName[200];
    TH1 *jetShapeHistogram = new TH1D(Form("%s_JetShape_%s",backgroundSubtractedHistogram->GetName(),extraword.Data()), Form("%s_JetShape_%s",backgroundSubtractedHistogram->GetName(),extraword.Data()), nRBins, RBins);

    sprintf(histogramName, "%sJetShapeCounts", backgroundSubtractedHistogram->GetName());
    TH1 *fhJetShapeCounts = (TH1 *)jetShapeHistogram->Clone(histogramName);

    sprintf(histogramName, "%sRbinMap", backgroundSubtractedHistogram->GetName());
    TH2 *fhJetShapeBinMap = (TH2 *)backgroundSubtractedHistogram->Clone(histogramName);

    for (int iPhiBin = 0; iPhiBin < fhJetShapeBinMap->GetNbinsX(); iPhiBin++)
    {
        for (int iEtaBin = 0; iEtaBin < fhJetShapeBinMap->GetNbinsY(); iEtaBin++)
        {
            fhJetShapeBinMap->SetBinContent(iPhiBin, iEtaBin, 0);
            fhJetShapeBinMap->SetBinError(iPhiBin, iEtaBin, 0);
        }
    }

    double binWidthDeltaPhi = backgroundSubtractedHistogram->GetXaxis()->GetBinWidth(1);
    double binWidthDeltaEta = backgroundSubtractedHistogram->GetYaxis()->GetBinWidth(1);
    double oneBinArea = binWidthDeltaPhi * binWidthDeltaEta;

    int minBinPhi = backgroundSubtractedHistogram->GetXaxis()->FindBin((-TMath::Pi() / 2.0) + 0.0001);
    int maxBinPhi = backgroundSubtractedHistogram->GetXaxis()->FindBin((TMath::Pi() / 2.0) - 0.0001);
    int minBinEta = backgroundSubtractedHistogram->GetYaxis()->FindBin((-TMath::Pi() / 2.0) + 0.0001);
    int maxBinEta = backgroundSubtractedHistogram->GetYaxis()->FindBin((TMath::Pi() / 2.0) - 0.0001);

    double deltaEta;   
    double deltaPhi;   
    double R;          
    int Rbin;          
    double newContent; 
    double newError;   
    double oldContent; 
    double oldError;   
    double oldCounts;  

    for (int iPhiBin = minBinPhi; iPhiBin <= maxBinPhi; iPhiBin++)
    {
        deltaPhi = backgroundSubtractedHistogram->GetXaxis()->GetBinCenter(iPhiBin); 

        for (int iEtaBin = minBinEta; iEtaBin <= maxBinEta; iEtaBin++)
        {
            deltaEta = backgroundSubtractedHistogram->GetYaxis()->GetBinCenter(iEtaBin); 

            R = TMath::Sqrt(TMath::Power(deltaPhi, 2) + TMath::Power(deltaEta, 2));
            Rbin = jetShapeHistogram->FindBin(R);

            newContent = backgroundSubtractedHistogram->GetBinContent(iPhiBin, iEtaBin) * oneBinArea;
            newError = backgroundSubtractedHistogram->GetBinError(iPhiBin, iEtaBin) * oneBinArea;

            oldContent = jetShapeHistogram->GetBinContent(Rbin);
            oldError = jetShapeHistogram->GetBinError(Rbin);
            oldCounts = fhJetShapeCounts->GetBinContent(Rbin);

            jetShapeHistogram->SetBinContent(Rbin, newContent + oldContent);
            jetShapeHistogram->SetBinError(Rbin, TMath::Sqrt(TMath::Power(newError, 2) + TMath::Power(oldError, 2)));

            fhJetShapeCounts->SetBinContent(Rbin, oldCounts + 1);
            fhJetShapeCounts->SetBinError(Rbin, 0);

            fhJetShapeBinMap->SetBinContent(iPhiBin, iEtaBin, Rbin);
        } 
    }    

    bool fJetShapeNormalizationMethod = kTRUE;
    if (fJetShapeNormalizationMethod == kTRUE)
    {
        jetShapeHistogram->Scale(1.0, "width");

        double totalBinArea;
        double ringArea;

        for (int iRBin = 1; iRBin <= jetShapeHistogram->GetNbinsX(); iRBin++)
        {
            totalBinArea = fhJetShapeCounts->GetBinContent(iRBin) * oneBinArea;
            ringArea = TMath::Pi() * RBins[iRBin] * RBins[iRBin] - TMath::Pi() * RBins[iRBin - 1] * RBins[iRBin - 1];
            jetShapeHistogram->SetBinContent(iRBin, jetShapeHistogram->GetBinContent(iRBin) * (ringArea / totalBinArea));
            jetShapeHistogram->SetBinError(iRBin, jetShapeHistogram->GetBinError(iRBin) * (ringArea / totalBinArea));
        }
    }

    return jetShapeHistogram;
}

TH1* Get_yield_vs_pt(vector<TH1*> jet_shapes_bypt, float dr_end = 0.399, float dr_start = 0.01)
{
    float ptbins_edges[5] = {1,2,3,4,10};
    Double_t Integral_error = 0;
    TH1D* Yield_vs_pt = new TH1D("Yield_vs_pt","Yield_vs_pt",4,ptbins_edges); Yield_vs_pt->Sumw2();
    TH1D* Yield_vs_pt_norm = new TH1D("Yield_vs_pt_norm","Yield_vs_pt_norm",4,ptbins_edges); Yield_vs_pt_norm->Sumw2();
    if (jet_shapes_bypt.size() < 4) {
        printf("Error: only %.1llu jet shapes were included \n",jet_shapes_bypt.size());
        return(Yield_vs_pt);
    }

    int first_bin = jet_shapes_bypt[0]->FindBin(dr_start);
    int last_bin =  jet_shapes_bypt[0]->FindBin(dr_end);

    for (int i=0; i<4; i++){
        Yield_vs_pt->SetBinContent(i+1,jet_shapes_bypt[i]->IntegralAndError(first_bin,last_bin, Integral_error,"width"));
        Yield_vs_pt->SetBinError(i+1,Integral_error);
        Yield_vs_pt_norm->SetBinContent(i+1, Yield_vs_pt_norm->GetBinWidth(i+1)/Yield_vs_pt_norm->GetBinWidth(1));
    }

    return(Yield_vs_pt);
}