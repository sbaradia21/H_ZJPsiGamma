#!/bin/sh  
python runAnalysis.py --yaml analysis_2016PreData_config.yml --outdirtag ZData_Nominal_PxSeed --send
python runAnalysis.py --yaml analysis_2016PostData_config.yml --outdirtag ZData_Nominal_PxSeed --send
python runAnalysis.py --yaml analysis_2017Data_config.yml --outdirtag ZData_Nominal_PxSeed --send
python runAnalysis.py --yaml analysis_2018Data_config.yml --outdirtag ZData_Nominal_PxSeed --send
