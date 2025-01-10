import re
import ROOT
import matplotlib.pyplot as plt
import numpy as np
import argparse

def plotHistogram(binCenters, normedSAFYields, normedSAFUncs, plotType, region, rebin_factor=1):
    if rebin_factor > 1:
        def rebin_data(data, factor):
            return np.mean(data.reshape(-1, factor), axis=1)

        def rebin_uncertainties(data, factor):
            return np.sqrt(np.sum((data.reshape(-1, factor)) ** 2, axis=1)) / factor

        binCenters = rebin_data(binCenters, rebin_factor)
        normedSAFYields = rebin_data(normedSAFYields, rebin_factor)
        normedSAFUncs = rebin_uncertainties(normedSAFUncs, rebin_factor)

    plt.figure(figsize=(10, 7))
    plt.errorbar(binCenters, normedSAFYields, yerr=normedSAFUncs, fmt='o', color='b', label=plotType, capsize=5, capthick=2)
    plt.xlabel('Bin Center')
    plt.ylabel('Normalized Yield')
    plt.title(plotType)
    plt.ylim(0, 1.2 * max(normedSAFYields))

    plt.savefig('Hadronic_'+args.arg1 + '/' + plotType + '_' + region + '.pdf')
    plt.savefig('Hadronic_'+args.arg1 + '/' + plotType + '_' + region + '.png', dpi=300)
    plt.close()

def plotRatio(inclusiveRootBins, inclusiveRootYields, inclusiveRootUncs, normedSAFYields, normedSAFUncs, plotType, region):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})

    ax1.errorbar(inclusiveRootBins, inclusiveRootYields, yerr=inclusiveRootUncs, fmt='o', color='b', label='Coffea (Normalized to RunII Lumi and Xsec)', capsize=5, capthick=2)
    ax1.errorbar(inclusiveRootBins, normedSAFYields, yerr=normedSAFUncs, fmt='o', color='r', label='MadAnalysis (Normalized to Coffea)', capsize=5, capthick=2)

    ax1.set_ylim(0, 1.2 * np.max(inclusiveRootYields))
    ax1.set_xlabel('Bin Center')
    ax1.set_ylabel('Yield')
    ax1.set_title(plotType)
    ax1.legend()

    ratio = np.array(normedSAFYields) / np.array(inclusiveRootYields)
    ratioUncertainty = ratio * np.sqrt((np.array(normedSAFUncs) / np.array(normedSAFYields))**2 + (np.array(normedSAFUncs) / np.array(normedSAFYields))**2)
    ax2.errorbar(inclusiveRootBins, ratio, yerr=ratioUncertainty, fmt='o', color='k', capsize=5, capthick=2)
    ax2.axhline(1, color='gray', linestyle='--')

    ax2.set_ylim(0.1, 2.5)
    ax2.set_xlabel('Bin Center')
    ax2.set_ylabel('MadAnalysis / Coffea')

    xLimits = [np.min(inclusiveRootBins), np.max(inclusiveRootBins)]
    ax1.set_xlim(xLimits)
    ax2.set_xlim(xLimits)

    plt.tight_layout()
    plt.savefig('Hadronic_'+args.arg1 + '/' + plotType + '_' + region + '.pdf')
    plt.savefig('Hadronic_'+args.arg1 + '/' + plotType + '_' + region + '.png', dpi=300)

def plotLeptonRatio(binCenters, muonYields, muonUncertainties, elecYields, elecUncertainties, kinematic, region, rebin_factor=1):
    if rebin_factor > 1:
        # Rebin the data by grouping according to the rebin_factor
        def rebin_data(data, factor):
            return np.mean(data.reshape(-1, factor), axis=1)

        def rebin_uncertainties(data, factor):
            return np.sqrt(np.sum((data.reshape(-1, factor)) ** 2, axis=1)) / factor

        binCenters = rebin_data(binCenters, rebin_factor)
        muonYields = rebin_data(muonYields, rebin_factor)
        muonUncertainties = rebin_uncertainties(muonUncertainties, rebin_factor)
        elecYields = rebin_data(elecYields, rebin_factor)
        elecUncertainties = rebin_uncertainties(elecUncertainties, rebin_factor)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})

    ax1.errorbar(binCenters, muonYields, yerr=muonUncertainties, fmt='o', color='b', label='Muon', capsize=5)
    ax1.errorbar(binCenters, elecYields, yerr=elecUncertainties, fmt='o', color='r', label='Electron', capsize=5)
    ax1.set_ylabel('Normalized Yield')
    ax1.set_title('Comparing Lepton ' + kinematic)
    ax1.legend()
    ax1.set_ylim(0, 1.2 * max(muonYields))

    ratio = np.divide(muonYields, elecYields, out=np.zeros_like(muonYields), where=elecYields != 0)
    ratio_uncertainty = ratio * np.sqrt((muonUncertainties / muonYields) ** 2 + (elecUncertainties / elecYields) ** 2)

    ax2.errorbar(binCenters, ratio, yerr=ratio_uncertainty, fmt='o', color='k', capsize=5)
    ax2.axhline(1, color='gray', linestyle='--')
    ax2.set_xlabel('Bin Center')
    ax2.set_ylabel('Muon / Electron Ratio')
    ax2.set_ylim(0.2, 4.0)

    xLimits = [np.min(binCenters), np.max(binCenters)]
    ax1.set_xlim(xLimits)
    ax2.set_xlim(xLimits)

    # Save the plot
    plt.tight_layout()
    plt.savefig('Hadronic_'+args.arg1 + '/' + 'LeptonComparison'+kinematic + '_' + region + '.pdf')
    plt.savefig('Hadronic_'+args.arg1 + '/' + 'LeptonComparison'+kinematic + '_' + region + '.png', dpi=300)
    plt.close()

def readRootFile(rootPath, histogramName):
    rootFile = ROOT.TFile(rootPath, "READ")
    histogram = rootFile.Get(histogramName)

    if not histogram:
        print(f"Histogram '{histogramName}' not found.")
        return None, None, None

    bin_centers = []
    yields = []
    uncertainties = []

    for i in range(1, histogram.GetNbinsX() + 1):
        bin_center = histogram.GetBinCenter(i)
        yield_value = histogram.GetBinContent(i)
        uncertainty_value = histogram.GetBinError(i)

        bin_centers.append(bin_center)
        yields.append(yield_value)
        uncertainties.append(uncertainty_value)

    return bin_centers, yields, uncertainties

def readSAFFile(safContent, subsectionName):
    pattern = rf'<Histo>[\s\S]*?<Description>\s*"{subsectionName}"[\s\S]*?<Data>([\s\S]*?)</Data>'
    match = re.search(pattern, safContent)

    if not match:
        print(f"Histogram '{subsectionName}' not found in SAF file.")
        return []

    dataLines = match.group(1).strip().splitlines()

    # Exclude the first line (underflow) and last line (overflow)
    if len(dataLines) > 2:
        dataLines = dataLines[1:-1]
    
    data = [float(line.split()[0]) for line in dataLines]

    return data

def readFiles(rootPath, histogramName, safContent, subsectionName):
    rootBinCenters, rootYields, rootUncs = readRootFile(rootPath, histogramName)
    safYields = readSAFFile(safContent, subsectionName)
    return rootBinCenters, rootYields, rootUncs, safYields

def rebinData(data, binsToMerge):
    rebinnedData = [sum(data[i:i+binsToMerge]) for i in range(0, len(data), binsToMerge)]
    binCenters = [i * binsToMerge + binsToMerge / 2 for i in range(len(rebinnedData))]
    return binCenters, rebinnedData

def main(arg1, arg2):
    inputFileSets = {
        "SR": {
            "leadclustertracks":        ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks.root",    f"leadclustertracks_SUEP_hadronic_mS125_{arg2}",             "LeadAk15NTRACKS",    4),
            "leadjet_pt":               ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadjet_pt.root",           f"leadjet_pt_SUEP_hadronic_mS125_{arg2}",                    "LeadAk4PT",          10),
            "leadclustertracks_A":      ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_SR.root", f"leadclustertracks_onecluster_SUEP_hadronic_mS125_{arg2}",  "ABCD_A",             5),
            "leadclustertracks_B1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_B1.root", f"leadclustertracks_B1_SUEP_hadronic_mS125_{arg2}",          "ABCD_B1",            5),
            "leadclustertracks_B2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_B2.root", f"leadclustertracks_B2_SUEP_hadronic_mS125_{arg2}",          "ABCD_B2",            5),
            "leadclustertracks_C1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_C1.root", f"leadclustertracks_C1_SUEP_hadronic_mS125_{arg2}",          "ABCD_C1",            1),
            "leadclustertracks_C2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_C2.root", f"leadclustertracks_C2_SUEP_hadronic_mS125_{arg2}",          "ABCD_C2",            1),
            "leadclustertracks_D1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_D1.root", f"leadclustertracks_D1_SUEP_hadronic_mS125_{arg2}",          "ABCD_D1",            1),
            "leadclustertracks_D2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_D2.root", f"leadclustertracks_D2_SUEP_hadronic_mS125_{arg2}",          "ABCD_D2",            1),
            "leadclustertracks_E1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_E1.root", f"leadclustertracks_E1_SUEP_hadronic_mS125_{arg2}",          "ABCD_E1",            1),
            "leadclustertracks_E2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/SR_RunII_CommonBounds_Fixed/leadclustertracks_E2.root", f"leadclustertracks_E2_SUEP_hadronic_mS125_{arg2}",          "ABCD_E2",            1),
            "ZMass":                    ("None", "None", "ZMass",           (61, 60.0, 120.0)),
            "ZpT":                      ("None", "None", "ZpT",             (501, 0.0, 500.0)),
            "ZETA":                     ("None", "None", "ZETA",            (61, -3.0, 3.0)),
            "ZPHI":                     ("None", "None", "ZPHI",            (61, -3.14, 3.14)),
            "LeadLepPT":                ("None", "None", "LeadLepPT",       (301, 0.0, 300.0)),
            "SubleadLepPT":             ("None", "None", "SubleadLepPT",    (301, 0.0, 300.0)),
            "LeadLepETA":               ("None", "None", "LeadLepETA",      (61, -3.0, 3.0)),
            "SubleadLepETA":            ("None", "None", "SubleadLepETA",   (61, -3.0, 3.0)),
            "LeadLepPHI":               ("None", "None", "LeadLepPHI",      (61, -3.14, 3.14)),
            "SubleadLepPHI":            ("None", "None", "SubleadLepPHI",   (61, -3.14, 3.14)),
            "LeadAk4ETA":               ("None", "None", "LeadAk4ETA",      (61, -3.0, 3.0)),
            "LeadAk4PHI":               ("None", "None", "LeadAk4PHI",      (61, -3.14, 3.14)),
            "LeadAk15PT":               ("None", "None", "LeadAk15PT",      (501, 0.0, 500.0)),
            "LeadAk15ETA":              ("None", "None", "LeadAk15ETA",     (61, -3.0, 3.0)),
            "LeadAk15PHI":              ("None", "None", "LeadAk15PHI",     (61, -3.14, 3.14)),
            "LeadAk15Mass":             ("None", "None", "LeadAk15Mass",    (401, 0, 400)),
            "LeadMuonPT":               ("None", "None", "LeadMuonPT",      (301, 0.0, 300.0)),
            "SubleadMuonPT":            ("None", "None", "SubleadMuonPT",   (301, 0.0, 300.0)),
            "LeadMuonETA":              ("None", "None", "LeadMuonETA",     (61, -3.0, 3.0)),
            "SubleadMuonETA":           ("None", "None", "SubleadMuonETA",  (61, -3.0, 3.0)),
            "LeadMuonPHI":              ("None", "None", "LeadMuonPHI",     (61, -3.14, 3.14)),
            "SubleadMuonPHI":           ("None", "None", "SubleadMuonPHI",  (61, -3.14, 3.14)),
            "LeadElecPT":               ("None", "None", "LeadElecPT",      (301, 0.0, 300.0)),
            "SubleadElecPT":            ("None", "None", "SubleadElecPT",   (301, 0.0, 300.0)),
            "LeadElecETA":              ("None", "None", "LeadElecETA",     (61, -3.0, 3.0)),
            "SubleadElecETA":           ("None", "None", "SubleadElecETA",  (61, -3.0, 3.0)),
            "LeadElecPHI":              ("None", "None", "LeadElecPHI",     (61, -3.14, 3.14)),
            "SubleadElecPHI":           ("None", "None", "SubleadElecPHI",  (61, -3.14, 3.14)),
            "safFile":                 (f"histos_{arg1}_Ak15Mad.saf")
        },

        "CRDY": {
            "leadclustertracks":        ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks.root",    f"leadclustertracks_SUEP_hadronic_mS125_{arg2}",             "LeadAk15NTRACKS",    4),
            "leadjet_pt":               ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadjet_pt.root",           f"leadjet_pt_SUEP_hadronic_mS125_{arg2}",                    "LeadAk4PT",          10),
            "leadclustertracks_A":      ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_SR.root", f"leadclustertracks_onecluster_SUEP_hadronic_mS125_{arg2}",  "ABCD_A",             5),
            "leadclustertracks_B1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_B1.root", f"leadclustertracks_B1_SUEP_hadronic_mS125_{arg2}",          "ABCD_B1",            5),
            "leadclustertracks_B2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_B2.root", f"leadclustertracks_B2_SUEP_hadronic_mS125_{arg2}",          "ABCD_B2",            5),
            "leadclustertracks_C1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_C1.root", f"leadclustertracks_C1_SUEP_hadronic_mS125_{arg2}",          "ABCD_C1",            1),
            "leadclustertracks_C2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_C2.root", f"leadclustertracks_C2_SUEP_hadronic_mS125_{arg2}",          "ABCD_C2",            1),
            "leadclustertracks_D1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_D1.root", f"leadclustertracks_D1_SUEP_hadronic_mS125_{arg2}",          "ABCD_D1",            1),
            "leadclustertracks_D2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_D2.root", f"leadclustertracks_D2_SUEP_hadronic_mS125_{arg2}",          "ABCD_D2",            1),
            "leadclustertracks_E1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_E1.root", f"leadclustertracks_E1_SUEP_hadronic_mS125_{arg2}",          "ABCD_E1",            1),
            "leadclustertracks_E2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRDY_RunII_CommonBounds_Fixed/leadclustertracks_E2.root", f"leadclustertracks_E2_SUEP_hadronic_mS125_{arg2}",          "ABCD_E2",            1),
            "ZMass":                    ("None", "None", "ZMass",           (61, 60.0, 120.0)),
            "ZpT":                      ("None", "None", "ZpT",             (501, 0.0, 500.0)),
            "ZETA":                     ("None", "None", "ZETA",            (61, -3.0, 3.0)),
            "ZPHI":                     ("None", "None", "ZPHI",            (61, -3.14, 3.14)),
            "LeadLepPT":                ("None", "None", "LeadLepPT",       (301, 0.0, 300.0)),
            "SubleadLepPT":             ("None", "None", "SubleadLepPT",    (301, 0.0, 300.0)),
            "LeadLepETA":               ("None", "None", "LeadLepETA",      (61, -3.0, 3.0)),
            "SubleadLepETA":            ("None", "None", "SubleadLepETA",   (61, -3.0, 3.0)),
            "LeadLepPHI":               ("None", "None", "LeadLepPHI",      (61, -3.14, 3.14)),
            "SubleadLepPHI":            ("None", "None", "SubleadLepPHI",   (61, -3.14, 3.14)),
            "LeadAk4ETA":               ("None", "None", "LeadAk4ETA",      (61, -3.0, 3.0)),
            "LeadAk4PHI":               ("None", "None", "LeadAk4PHI",      (61, -3.14, 3.14)),
            "LeadAk15PT":               ("None", "None", "LeadAk15PT",      (501, 0.0, 500.0)),
            "LeadAk15ETA":              ("None", "None", "LeadAk15ETA",     (61, -3.0, 3.0)),
            "LeadAk15PHI":              ("None", "None", "LeadAk15PHI",     (61, -3.14, 3.14)),
            "LeadAk15Mass":             ("None", "None", "LeadAk15Mass",    (401, 0, 400)),
            "LeadMuonPT":               ("None", "None", "LeadMuonPT",      (301, 0.0, 300.0)),
            "SubleadMuonPT":            ("None", "None", "SubleadMuonPT",   (301, 0.0, 300.0)),
            "LeadMuonETA":              ("None", "None", "LeadMuonETA",     (61, -3.0, 3.0)),
            "SubleadMuonETA":           ("None", "None", "SubleadMuonETA",  (61, -3.0, 3.0)),
            "LeadMuonPHI":              ("None", "None", "LeadMuonPHI",     (61, -3.14, 3.14)),
            "SubleadMuonPHI":           ("None", "None", "SubleadMuonPHI",  (61, -3.14, 3.14)),
            "LeadElecPT":               ("None", "None", "LeadElecPT",      (301, 0.0, 300.0)),
            "SubleadElecPT":            ("None", "None", "SubleadElecPT",   (301, 0.0, 300.0)),
            "LeadElecETA":              ("None", "None", "LeadElecETA",     (61, -3.0, 3.0)),
            "SubleadElecETA":           ("None", "None", "SubleadElecETA",  (61, -3.0, 3.0)),
            "LeadElecPHI":              ("None", "None", "LeadElecPHI",     (61, -3.14, 3.14)),
            "SubleadElecPHI":           ("None", "None", "SubleadElecPHI",  (61, -3.14, 3.14)),
            "safFile":                 (f"histos_{arg1}_Ak15Mad_CRDY.saf")
        },

        "CRTT": {
            "leadclustertracks":        ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks.root",    f"leadclustertracks_SUEP_hadronic_mS125_{arg2}",             "LeadAk15NTRACKS",    4),
            "leadjet_pt":               ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadjet_pt.root",           f"leadjet_pt_SUEP_hadronic_mS125_{arg2}",                    "LeadAk4PT",          10),
            "leadclustertracks_A":      ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_SR.root", f"leadclustertracks_onecluster_SUEP_hadronic_mS125_{arg2}",  "ABCD_A",             5),
            "leadclustertracks_B1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_B1.root", f"leadclustertracks_B1_SUEP_hadronic_mS125_{arg2}",          "ABCD_B1",            5),
            "leadclustertracks_B2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_B2.root", f"leadclustertracks_B2_SUEP_hadronic_mS125_{arg2}",          "ABCD_B2",            5),
            "leadclustertracks_C1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_C1.root", f"leadclustertracks_C1_SUEP_hadronic_mS125_{arg2}",          "ABCD_C1",            1),
            "leadclustertracks_C2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_C2.root", f"leadclustertracks_C2_SUEP_hadronic_mS125_{arg2}",          "ABCD_C2",            1),
            "leadclustertracks_D1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_D1.root", f"leadclustertracks_D1_SUEP_hadronic_mS125_{arg2}",          "ABCD_D1",            1),
            "leadclustertracks_D2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_D2.root", f"leadclustertracks_D2_SUEP_hadronic_mS125_{arg2}",          "ABCD_D2",            1),
            "leadclustertracks_E1":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_E1.root", f"leadclustertracks_E1_SUEP_hadronic_mS125_{arg2}",          "ABCD_E1",            1),
            "leadclustertracks_E2":     ("/eos/user/g/gdecastr/SUEPCoffea_dask/Plots/CRTT_RunII_CommonBounds_Fixed/leadclustertracks_E2.root", f"leadclustertracks_E2_SUEP_hadronic_mS125_{arg2}",          "ABCD_E2",            1),
            "ZMass":                    ("None", "None", "ZMass",           (61, 60.0, 120.0)),
            "ZpT":                      ("None", "None", "ZpT",             (501, 0.0, 500.0)),
            "ZETA":                     ("None", "None", "ZETA",            (61, -3.0, 3.0)),
            "ZPHI":                     ("None", "None", "ZPHI",            (61, -3.14, 3.14)),
            "LeadLepPT":                ("None", "None", "LeadLepPT",       (301, 0.0, 300.0)),
            "SubleadLepPT":             ("None", "None", "SubleadLepPT",    (301, 0.0, 300.0)),
            "LeadLepETA":               ("None", "None", "LeadLepETA",      (61, -3.0, 3.0)),
            "SubleadLepETA":            ("None", "None", "SubleadLepETA",   (61, -3.0, 3.0)),
            "LeadLepPHI":               ("None", "None", "LeadLepPHI",      (61, -3.14, 3.14)),
            "SubleadLepPHI":            ("None", "None", "SubleadLepPHI",   (61, -3.14, 3.14)),
            "LeadAk4ETA":               ("None", "None", "LeadAk4ETA",      (61, -3.0, 3.0)),
            "LeadAk4PHI":               ("None", "None", "LeadAk4PHI",      (61, -3.14, 3.14)),
            "LeadAk15PT":               ("None", "None", "LeadAk15PT",      (501, 0.0, 500.0)),
            "LeadAk15ETA":              ("None", "None", "LeadAk15ETA",     (61, -3.0, 3.0)),
            "LeadAk15PHI":              ("None", "None", "LeadAk15PHI",     (61, -3.14, 3.14)),
            "LeadAk15Mass":             ("None", "None", "LeadAk15Mass",    (401, 0, 400)),
            "LeadMuonPT":               ("None", "None", "LeadMuonPT",      (301, 0.0, 300.0)),
            "SubleadMuonPT":            ("None", "None", "SubleadMuonPT",   (301, 0.0, 300.0)),
            "LeadMuonETA":              ("None", "None", "LeadMuonETA",     (61, -3.0, 3.0)),
            "SubleadMuonETA":           ("None", "None", "SubleadMuonETA",  (61, -3.0, 3.0)),
            "LeadMuonPHI":              ("None", "None", "LeadMuonPHI",     (61, -3.14, 3.14)),
            "SubleadMuonPHI":           ("None", "None", "SubleadMuonPHI",  (61, -3.14, 3.14)),
            "LeadElecPT":               ("None", "None", "LeadElecPT",      (301, 0.0, 300.0)),
            "SubleadElecPT":            ("None", "None", "SubleadElecPT",   (301, 0.0, 300.0)),
            "LeadElecETA":              ("None", "None", "LeadElecETA",     (61, -3.0, 3.0)),
            "SubleadElecETA":           ("None", "None", "SubleadElecETA",  (61, -3.0, 3.0)),
            "LeadElecPHI":              ("None", "None", "LeadElecPHI",     (61, -3.14, 3.14)),
            "SubleadElecPHI":           ("None", "None", "SubleadElecPHI",  (61, -3.14, 3.14)),
            "safFile":                 (f"histos_{arg1}_Ak15Mad_CRTT.saf")
        },
    }

    for region, plottingInfo in inputFileSets.items():
        with open(plottingInfo["safFile"], 'r') as safFile:
            safContent = safFile.read()

        ABCD_histograms = {}
        Comparison_histograms = {}
        Solo_histograms = {}

        inclusiveTracks = plottingInfo["leadclustertracks"]
        inclusiveRootBins, inclusiveRootYields, inclusiveRootUncs, SAFYields = readFiles(inclusiveTracks[0], inclusiveTracks[1], safContent, inclusiveTracks[2])

        normFactor = np.sum(inclusiveRootYields) / np.sum(SAFYields)
        rebinnedSAFYields = rebinData(SAFYields, inclusiveTracks[3])

        normedSAFYields = np.array(rebinnedSAFYields) * normFactor
        normedSAFUncs = np.sqrt(rebinnedSAFYields) * normFactor

        Comparison_histograms["leadclustertracks"] = [inclusiveRootBins, inclusiveRootYields, inclusiveRootUncs, normedSAFYields, normedSAFUncs]

        for plotType, plotDetails in plottingInfo.items():
            if plotType == "safFile":
                continue
            filePath, treeName, safSubsection, rebin = plotDetails
            if filePath.endswith(".root"):
                rootBins, rootYields, rootUncs, SAFYields = readFiles(filePath, treeName, safContent, safSubsection)
                rebinnedSAFYields = rebinData(SAFYields, rebin)
                normedSAFYields = np.array(rebinnedSAFYields) * normFactor
                normedSAFUncs = np.sqrt(rebinnedSAFYields) * normFactor
                if 'ABCD' in safSubsection:
                    ABCD_histograms[plotType] = [rootBins, rootYields, rootUncs, normedSAFYields, normedSAFUncs]
                else:
                    Comparison_histograms[plotType] = [rootBins, rootYields, rootUncs, normedSAFYields, normedSAFUncs]
                continue
            else:
                SAFYields = readSAFFile(safContent, safSubsection)
                normedSAFYields = np.array(SAFYields) * normFactor
                normedSAFUncs = np.sqrt(SAFYields) * normFactor
                Solo_histograms[plotType] = [normedSAFYields, normedSAFUncs, rebin]

        for plotType, dataToPlot in Solo_histograms.items():
            normedSAFYields, normedSAFUncs, rebin = dataToPlot
            binEdges = np.linspace(rebin[1], rebin[2], rebin[0])
            binCenters = (binEdges[:-1] + binEdges[1:]) / 2
            plotHistogram(binCenters, normedSAFYields, normedSAFUncs, plotType, region, 2)

        for plotType, dataToPlot in Comparison_histograms.items():
            inclusiveRootBins, inclusiveRootYields, inclusiveRootUncs, normedSAFYields, normedSAFUncs = dataToPlot
            plotRatio(inclusiveRootBins, inclusiveRootYields, inclusiveRootUncs, normedSAFYields[1], normedSAFUncs[1], plotType, region)

        for kinematic in ["PT", "ETA", "PHI"]:
            elecSAFYields, elecSAFUncs, rebin = Solo_histograms['LeadElec'+kinematic]
            muonSAFYields, muonSAFUncs, rebin = Solo_histograms['LeadMuon'+kinematic]
            binEdges = np.linspace(rebin[1], rebin[2], rebin[0])
            binCenters = (binEdges[:-1] + binEdges[1:]) / 2       
            plotLeptonRatio(binCenters, muonSAFYields, muonSAFUncs, elecSAFYields, elecSAFUncs, kinematic, region, 5)
        
        manualOrder = ["leadclustertracks_E2", "leadclustertracks_E1", "leadclustertracks_B2", "leadclustertracks_D2", "leadclustertracks_D1", "leadclustertracks_B1", "leadclustertracks_C2", "leadclustertracks_C1", "leadclustertracks_A"]
        coolRootYields, coolRootUncs, coolSAFYields, coolSAFUncs = [], [], [], []
        for plotType in manualOrder:
            rootBins, rootYields, rootUncs, SAFYields, SAFUncs = ABCD_histograms[plotType]
            SAFYields = SAFYields[1]
            SAFUncs = SAFUncs[1]
            if ("B2" not in plotType) and ("B1" not in plotType) and ("A" not in plotType):
                coolRootYields.append(rootYields[0])
                coolRootUncs.append(rootUncs[0])
                coolSAFYields.append(SAFYields[0])
                coolSAFUncs.append(SAFUncs[0])
            else:
                for i in range(7):
                    coolSAFYields.append(SAFYields[i])
                    coolSAFUncs.append(SAFUncs[i])
                coolSAFYields.append(sum(SAFYields[7:]))
                coolSAFUncs.append(np.sqrt(sum([unc**2 for unc in SAFUncs[7:]])))
                for i in range(len(rootYields)):
                    coolRootYields.append(rootYields[i])
                    coolRootUncs.append(rootUncs[i])

        plotRatio(np.linspace(0,30,30), coolRootYields, coolRootUncs, coolSAFYields, coolSAFUncs, "COOL", region)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for processing ROOT and SAF files.")
    parser.add_argument("arg1", type=str, help="Tag for SAF files and output folder naming")
    parser.add_argument("arg2", type=str, help="Tag for histogram name in ROOT files")

    args = parser.parse_args()

    main(args.arg1, args.arg2)