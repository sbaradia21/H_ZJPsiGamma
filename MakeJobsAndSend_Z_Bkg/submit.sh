#!/bin/sh  
python runAnalysis.py --yaml analysis_2016PreBkgMC_config.yml --outdirtag ZBkg_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
python runAnalysis.py --yaml analysis_2016PostBkgMC_config.yml --outdirtag ZBkg_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
python runAnalysis.py --yaml analysis_2017BkgMC_config.yml --outdirtag ZBkg_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
python runAnalysis.py --yaml analysis_2018BkgMC_config.yml --outdirtag ZBkg_Nominal_CutFlowCorrected_ElectronVeto_SFCorr --send
