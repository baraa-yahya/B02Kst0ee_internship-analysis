import sys
import os
import pandas as pd
import numpy as np
import uproot as ur
import logging
import argparse
import awkward as ak

sys.path.append('/sps/lhcb/yahya/setup/ewp-bd2ksteeangular-legacy/samples_and_selection/kpipi_reweighting')
import utils_wrkpipi

parser = argparse.ArgumentParser(description='Configuration of the parameters')
parser.add_argument("-f", "--MC_file_path", dest="MC_file_path", required=True, help="Path to the directory containing the MC root files")    
parser.add_argument("-m", "--mags", dest="mags", nargs='+', required=True, help="Choose magnet polarities", choices=['MagDown', "MagUp"])
parser.add_argument("-y", "--years", dest="years", nargs='+', required=True, help="Choose years", choices=['2011', "2012", '2015', '2016', '2017', '2018'])

args = parser.parse_args()
MC_file_path = args.MC_file_path
mags = args.mags
years = args.years

mc_branches = ['runNumber', 'eventNumber', 'nCandidate', 'Bp_TRUEETA', 'Bp_TRUEPHI', 'Bp_TRUEP_E', 
               'K_TRUEP_E', 'K_TRUEP_X', 'K_TRUEP_Y', 'K_TRUEP_Z', 'K_TRUEETA', 'K_TRUEPHI',
               'Pi1_TRUEP_E', 'Pi1_TRUEP_X', 'Pi1_TRUEP_Y', 'Pi1_TRUEP_Z', 'Pi1_TRUEETA', 'Pi1_TRUEPHI',
               'Pi2_TRUEP_E', 'Pi2_TRUEP_X', 'Pi2_TRUEP_Y', 'Pi2_TRUEP_Z', 'Pi2_TRUEETA', 'Pi2_TRUEPHI',
               'E1_TRUEP_E', 'E1_TRUEP_X', 'E1_TRUEP_Y', 'E1_TRUEP_Z', 'E1_TRUEETA', 'E1_TRUEPHI',
               'E2_TRUEP_E', 'E2_TRUEP_X', 'E2_TRUEP_Y', 'E2_TRUEP_Z', 'E2_TRUEETA', 'E2_TRUEPHI','kpipi_weight']
DT_branches = [
    'runNumber','eventNumber',
    'B_eta', 'B_phi', 'B_TRUEP_E', 'ctl_lhcb', 'ctk_lhcb', 'phi_lhcb', 'B_DTF_PVandB_B_M_0', 'B_DTF_PVandB_Kstar_M_0'
]

for year in years:
    for mag in mags:
        ifile_MCDT = f"/sps/lhcb/yahya/Internship/test_Bu/Bu2Kpipiee-MCDT-{year}-{mag}-RKpipiWeights-with_angles.root"
        MC_file_name = os.path.join(MC_file_path, f"Bu2Kpipiee-MC-{year}-{mag}-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput.root")
        ifnout = MC_file_name.replace(".root", "-test2.root")
        
        dfMCDT = ur.open(ifile_MCDT)["MCDecayTree"].arrays(expressions=mc_branches, library='pd')
        logging.warning("In case of multiple candidates in MCDecayTree, select the one that matches the event in the DecayTree")

        mult_cand_events = dfMCDT.query('nCandidate != 0')
        dfDT = ur.open(MC_file_name)["DecayTree"].arrays(expressions=DT_branches, library='pd')

        for i in range(len(mult_cand_events)):
            matched_indices = np.where((dfDT["eventNumber"]==mult_cand_events["eventNumber"].iloc[i]) & (dfDT["runNumber"]==mult_cand_events["runNumber"].iloc[i]))[0]

            if len(matched_indices) == 0:
                dfMCDT = dfMCDT.query(f"~((eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]}))", engine='python')
            elif len(matched_indices) == 1:
                matched_MCDT = dfMCDT.query(f"(eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]})", engine='python')
                matched_DT = dfDT.query(f"(eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]})", engine='python')

                true_eta_0 = matched_MCDT["Bp_TRUEETA"].iloc[0]
                true_eta_1 = matched_MCDT["Bp_TRUEETA"].iloc[1]
                true_phi_0 = matched_MCDT["Bp_TRUEPHI"].iloc[0]
                true_phi_1 = matched_MCDT["Bp_TRUEPHI"].iloc[1]
                eta = matched_DT['B_eta'].iloc[0]
                phi = matched_DT['B_phi'].iloc[0]

                delta_0 = abs(true_eta_0-eta) + abs(true_phi_0-phi)
                delta_1 = abs(true_eta_1-eta) + abs(true_phi_1-phi)

                if delta_0 < delta_1:
                    dfMCDT = dfMCDT.query(f"~((eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]}) & (nCandidate==1))", engine='python')
                elif delta_1 < delta_0:
                    dfMCDT = dfMCDT.query(f"~((eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]}) & (nCandidate==0))", engine='python')
                else:
                    print("Case not implemented, improve code")
                    sys.exit(1)
            elif len(matched_indices) > 1:
                matched_MCDT = dfMCDT.query(f"(eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]})", engine='python')
                matched_DT = dfDT.query(f"(eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]})", engine='python')

                true_eta_0 = matched_MCDT["Bp_TRUEETA"].iloc[0]
                true_eta_1 = matched_MCDT["Bp_TRUEETA"].iloc[1]
                true_phi_0 = matched_MCDT["Bp_TRUEPHI"].iloc[0]
                true_phi_1 = matched_MCDT["Bp_TRUEPHI"].iloc[1]

                choices = []
                for j in range(len(matched_indices)):
                    eta = matched_DT['B_eta'].iloc[j]
                    phi = matched_DT['B_phi'].iloc[j]

                    delta_0 = abs(true_eta_0-eta) + abs(true_phi_0-phi)
                    delta_1 = abs(true_eta_1-eta) + abs(true_phi_1-phi)

                    if delta_0 < delta_1:
                        choices.append(0)
                    elif delta_1 < delta_0:
                        choices.append(1)
                    else:
                        print("Case not implemented, improve code")
                        sys.exit(1)

                if (len(set(choices)) == 1):
                    if 0 in set(choices):
                        dfMCDT = dfMCDT.query(f"~((eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]}) & (nCandidate==1))", engine='python')
                    elif 1 in set(choices):
                        dfMCDT = dfMCDT.query(f"~((eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]}) & (nCandidate==0))", engine='python')
                    else:
                        print("Problem")
                        sys.exit(1)
                else:
                    dfMCDT = dfMCDT.query(f"~((eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]}) & (nCandidate==1))", engine='python')
                    remaining_index = dfMCDT.query(f"((eventNumber=={mult_cand_events['eventNumber'].iloc[i]}) & (runNumber=={mult_cand_events['runNumber'].iloc[i]}))", engine='python').index[0]
                    dfMCDT.at[remaining_index, "kpipi_weight"] = 0

        dfMCDT = dfMCDT[~dfMCDT.duplicated(subset=None, keep='first')]
        dfMCDT = dfMCDT.drop(columns=["nCandidate", 'Bp_TRUEETA', 'Bp_TRUEPHI', 'Bp_TRUEP_E'])

        logging.warning("Matching weights of MCDT to DT")

        # Iterate over the DataFrame chunks
        for df_i_chunk in ur.iterate(MC_file_name+":DecayTree", step_size=50000, library='pd'):
            matched_df_i = df_i_chunk.join(dfMCDT.set_index(["runNumber", "eventNumber"]), on=["runNumber", "eventNumber"], how='left', lsuffix='_left', rsuffix='_NEW')
            assert matched_df_i.shape[0] == df_i_chunk.shape[0]

            # Convert DataFrame columns to lists for uproot
            arrays = {name: matched_df_i[name].values for name in matched_df_i.columns}

            # Create a new ROOT file and add the DecayTree branch to it
            with ur.recreate(ifnout) as f:
                f["DecayTree"] = arrays
