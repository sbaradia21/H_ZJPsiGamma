#!/bin/sh  
python runAnalysis.py --yaml analysis_2016PreSigMC_config.yml --outdirtag ZSignal_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
python runAnalysis.py --yaml analysis_2016PostSigMC_config.yml --outdirtag ZSignal_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
python runAnalysis.py --yaml analysis_2017SigMC_config.yml --outdirtag ZSignal_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
python runAnalysis.py --yaml analysis_2018SigMC_config.yml --outdirtag ZSignal_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
