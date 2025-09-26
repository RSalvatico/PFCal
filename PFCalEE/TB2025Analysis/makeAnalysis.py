import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob

ROOT.gROOT.ProcessLineSync(".L ../userlib/src/HGCSSSamplingSection.cc+")


def plotAbsorberE(input_files, output_dir, conversion_factor, energy, setup):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if setup == '140':
        custom_length = 25
    elif setup == '141':
        custom_length = 22
    else:
        custom_length = 27

    n_events = 0
    super_total_absorberE = []
    super_total_measuredE = []
    total_absorberE = [0] * custom_length
    total_measuredE = [0] * custom_length
    layers_absorberE = {x : [] for x in range(1, custom_length + 1)}
    layers_measuredE = {x : [] for x in range(1, custom_length + 1)}

    for file in input_files:
        print(f"Processing file {file}")
        fIn = ROOT.TFile.Open(file)
        tIn = fIn.Get("HGCSSTree")
        n_events += tIn.GetEntries()
        for entry in range(tIn.GetEntries()):
            tIn.GetEntry(entry)
            samplingSectionVec = tIn.HGCSSSamplingSectionVec
            n_layers = len(samplingSectionVec)
            absorberE = []
            measuredE = []
            for i in range(0, n_layers):
                absorberE.append(samplingSectionVec[i].absorberE())
                layers_absorberE[i+1].append(samplingSectionVec[i].absorberE())

                if i == 0:
                    layers_absorberE[i+1].append(samplingSectionVec[i].absorberE()+samplingSectionVec[i+1].absorberE()+samplingSectionVec[i+2].absorberE())
                    measuredE.append((samplingSectionVec[i].measuredE()+samplingSectionVec[i+1].measuredE()+samplingSectionVec[i+2].measuredE()) * conversion_factor / 2.0 * (samplingSectionVec[i].voldEdx()+samplingSectionVec[i+1].voldEdx()+samplingSectionVec[i+2].voldEdx()))
                    layers_measuredE[i+1].append( (samplingSectionVec[i].measuredE()+samplingSectionVec[i+1].measuredE()) * conversion_factor / 2.0 * samplingSectionVec[i].voldEdx() )
                    continue
                #elif i%2 == 1:
                #    measuredE.append((samplingSectionVec[i].measuredE()+samplingSectionVec[i+1].measuredE()) * conversion_factor / 2.0 * (samplingSectionVec[i].voldEdx()+samplingSectionVec[i+1].voldEdx()))
                #    layers_measuredE[i+1].append( (samplingSectionVec[i].measuredE()+samplingSectionVec[i+1].measuredE()) * conversion_factor / 2.0 * (samplingSectionVec[i].voldEdx()+samplingSectionVec[i+1].voldEdx()) )
                #elif i%2 == 0:
                #    measuredE.append(0)
                #    layers_measuredE[i+1].append(0)

            super_total_absorberE.append(sum(absorberE))
            super_total_measuredE.append(sum(measuredE))
            for i in range(0, n_layers-1):
                total_absorberE[i] += absorberE[i]
                total_measuredE[i] += measuredE[i]


    for layer in layers_absorberE:
        plt.figure()
        hist_range = (0, 70000) if layer==1 else (0, 15000)
        plt.hist(layers_absorberE[layer], bins=80, range=hist_range, histtype='step', 
             color='red', label='Absorber Energy', linewidth=2, alpha=0.8)
        plt.hist(layers_measuredE[layer], bins=80, range=hist_range, histtype='step', 
             color='blue', label='Calculated Energy', linewidth=2, alpha=0.8)
        plt.xlabel('Energy [MeV]')
        plt.ylabel('Entries')
    
        stats_text = (f"Absorber: μ={np.mean(layers_absorberE[layer]):.1f}, "
                      f"σ={np.std(layers_absorberE[layer]):.1f}\n"
                      f"Calculated: μ={np.mean(layers_measuredE[layer]):.1f}, "
                      f"σ={np.std(layers_measuredE[layer]):.1f}")

        plt.text(0.05, 0.95, stats_text,
                 transform=plt.gca().transAxes, verticalalignment='top',
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))
        plt.text(0.05, 0.75, f"Layer {layer}",
                 transform=plt.gca().transAxes, verticalalignment='top',
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))

        plt.legend(loc='upper right')
        plt.grid(True, alpha=0.3)
        print(f"Saving absorberE_layer{layer}_{setup}_{energy}GeV.png/pdf in {output_dir}")
        plt.savefig(os.path.join(output_dir, f'absorberE_layer{layer}_{setup}_{energy}GeV.png'))
        plt.savefig(os.path.join(output_dir, f'absorberE_layer{layer}_{setup}_{energy}GeV.pdf'))
        plt.close()

    plt.figure()
    plt.hist(super_total_absorberE, bins=80, range=(0, 95000), histtype='step', color='red', label='Absorber Energy', linewidth=2, alpha=0.8)
    plt.hist(super_total_measuredE, bins=80, range=(0, 95000), histtype='step', color='blue', label='Calculated Energy', linewidth=2, alpha=0.8)
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Entries')

    stats_text = (f"Absorber: μ={np.mean(super_total_absorberE):.1f}, "
                  f"σ={np.std(super_total_absorberE):.1f}\n"
                  f"Calculated: μ={np.mean(super_total_measuredE):.1f}, "
                  f"σ={np.std(super_total_measuredE):.1f}")

    plt.text(0.05, 0.95, stats_text,
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))
    plt.legend(loc='upper right')
    plt.grid(True, alpha=0.3)
    print(f"Saving super_total_absorberE_{setup}_{energy}GeV.png/pdf in {output_dir}")
    plt.savefig(os.path.join(output_dir, f'super_total_absorberE_{setup}_{energy}GeV.png'))
    plt.savefig(os.path.join(output_dir, f'super_total_absorberE_{setup}_{energy}GeV.pdf'))
    plt.close()

    # Third tipe of plot
    total_absorberE = [x / n_events for x in total_absorberE]
    total_measuredE = [x / n_events for x in total_measuredE]
    ratio = [ (total_absorberE[i]/total_measuredE[i]) if total_measuredE[i] != 0 else 0 for i in range(len(total_measuredE))]
    #print(total_absorberE)
    #print(total_measuredE)
    #print(ratio)
    #print(super_total_measuredE)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, 
                                   gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.1})
    
    ax1.plot(range(1, n_layers), total_absorberE, 'ro', label='Absorber Energy', markersize=6)
    ax1.plot(range(1, n_layers), total_measuredE, 'b^', label='Calculated Energy', markersize=6)
    ax1.set_ylabel('Energy [MeV]')
    ax1.set_ylim(top=50000)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(range(1, n_layers), ratio, 'gs', markersize=6)
    ax2.set_xlabel('Layer Number')
    ax2.set_ylabel('Meas/Calc')
    ax2.set_ylim(0.8, 2.0)  
    ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5) 
    ax2.grid(True, alpha=0.3)

    print(f"Saving absorberE_{setup}_{energy}GeV.png/pdf in {output_dir}")
    plt.savefig(os.path.join(output_dir, f'absorberE_{setup}_{energy}GeV.png'))
    plt.savefig(os.path.join(output_dir, f'absorberE_{setup}_{energy}GeV.pdf'))

    print("All done!")


def find_files(input_dir, pattern):

    root_files = glob.glob(os.path.join(input_dir, f'*{pattern}*.root'))
    return root_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare energies released in the absorber via calculation and pure G4 simulation')
    parser.add_argument('-i', '--input-dir', default='', required=True, help='List of input ROOT files.')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory for results.')
    parser.add_argument('-p', '--pattern', default='', help='Pattern to match input ROOT files.')
    parser.add_argument('-e', '--energy', default='', help='Beam energy in GeV.')
    parser.add_argument('-s', '--setup', default='', help='TB setup number.')
    args = parser.parse_args()

    conversion_Mev_fC = 1000000./3.6*0.00016
    conversion_fC_MIP = 1./5.15
    conversion_factor = conversion_Mev_fC * conversion_fC_MIP
    plotAbsorberE(find_files(args.input_dir, args.pattern), args.output_dir, conversion_factor, args.energy, args.setup)