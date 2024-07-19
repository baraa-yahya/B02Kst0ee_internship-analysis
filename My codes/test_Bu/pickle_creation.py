import os
import pandas as pd
import uproot
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

env = "/sps/lhcb/yahya/Internship/test_Bu/new_root_files"

def read_root_files_to_df(file_paths, branches):
    dfs = []
    for file_path in file_paths:
        with uproot.open(file_path) as file:
            tree = file["DecayTree"]
            arrays = tree.arrays(branches, library="pd")
            # Convert arrays to a DataFrame
            df = pd.DataFrame(arrays)
            dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def produce_pkl(decay, mode, years):
    print("Decay:               {}".format(decay))
    print("Mode:                {}".format(mode))
    print("Years:               {}".format(years))

    selection = "(GenericPresel==1) & (GenericPresel_Additional==1)& (MeerkatPresel_Tight==1) & (PIDPresel==1) & (TighterKst0Presel==1) & (TriggerPresel==1) & (VetoesPresel==1) & (VetoesPresel_Additional==1) & (CloneVeto==1) & (q2_NB>14.0)"
    branches = [
        "GenericPresel",
        "GenericPresel_Additional",
        "MeerkatPresel_Tight",
        "PIDPresel",
        "TighterKst0Presel",
        "TriggerPresel",
        "VetoesPresel",
        "VetoesPresel_Additional",
        "VetoesPresel_Bplus",
        "CloneVeto",
        "XGBOutput",
        "L1_L0ElectronDecision_TOS",
        "L2_L0ElectronDecision_TOS",
        "L1_L0MuonDecision_TOS",
        "L2_L0MuonDecision_TOS",
        "L1_L0Calo_ECAL_realET",
        "L2_L0Calo_ECAL_realET",
        "q2_DTF_PV",
        "B_DTF_PV_B_M_0",
        "B_DTF_PVandJpsi_B_M_0",
        "B_DTF_PVandPsi2S_B_M_0",
        "B_DTF_PV_Kstar_M_0",
        "B_L0Global_TIS",
        "L0En",
        "ctl_lhcb",
        "ctk_lhcb",
        "phi_lhcb",
        "XGBOutput_q2BDT",
        "XGBOutput_q2BDT_jpsi",
        "Year",
        "q2_NB"
    ]

    if mode == "MC":
        branches.extend(["w_Track_Trig_Kin", "TruthMatch", "B0_BKGCAT", "L0I"])
        selection += " & (TruthMatch==1)"
        if decay == "Bu2Kpipiee":
            logging.info("Extending selection with additional_selection_bu2kpipiee")
            branches.extend(["kpipi_weight_left", "kpipi_weight_right"])

    magnets = ["MagUp", "MagDown"]

    dfs = []
    for year in years:
        l_data = []
        if mode == "Data":
            data_path = os.path.join(env, "Data", decay, year, "Flagged")
        elif mode == "MC":
            data_path = os.path.join(env, "MC", decay, year, "Flagged")

        for mag in magnets:
            if mode == "Data":
                file_name = f"{decay}-Data-{year}-{mag}-Translated-Flagged-XGBOutput-q2BDTOutput.root"
            elif mode == "MC":
                file_name = f"{decay}-MC-{year}-{mag}-RXv9-ForPIDCalib-TrackEffCorr-nTracksCorr-PIDMeerkat-HLTaligned-Flagged-L0TrigEffCorr-HLTTrigEffCorr-KinRew2D-MCweights-XGBOutput-q2BDTOutput-test.root"

            l_data.append(os.path.join(data_path, file_name))

        filtered_df = read_root_files_to_df(l_data, branches)

        # Apply selection criteria
        filtered_df = filtered_df.query(selection)
        
        dfs.append(filtered_df)

    combined_df = pd.concat(dfs, ignore_index=True)

    # Ensure the pickle directory exists
    os.makedirs("pkl", exist_ok=True)

    # Save combined DataFrame to a pickle file
    pkl_name = f"{decay}_all_years_{mode}_new_root_files_no_q2_cut"
    pkl_file = os.path.join("pkl", f"{pkl_name}.pkl")
    combined_df.to_pickle(pkl_file)
    print("Producing combined pkl file {}".format(pkl_file))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input configuration for producing combined pickle files')
    parser.add_argument('-d', '--decay', dest='decay', required=True,
                        help='Select the decay, e.g. B02Kst0Jpsi2ee, B02Kst0Jpsi2eeSS, B02Kst0MmEp, B02Kst0MpEm, etc.')
    parser.add_argument('-y', '--years', dest='years', nargs='+', required=True,
                        help='Select the years, e.g. 2011 2012 2015 2016 2017 2018')
    parser.add_argument('-m', '--mode', dest='mode', required=True, choices=['Data', 'MC'],
                        help='Select the mode: Data or MC')

    args = parser.parse_args()
    decay = args.decay
    mode = args.mode
    years = args.years

    produce_pkl(decay, mode, years)
