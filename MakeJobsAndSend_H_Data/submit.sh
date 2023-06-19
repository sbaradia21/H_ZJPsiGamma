#!/bin/sh  
python runAnalysis.py --yaml analysis_2016Pre_config.yml --outdirtag HData_Nominal --send
python runAnalysis.py --yaml analysis_2016Post_config.yml --outdirtag HData_Nominal --send
python runAnalysis.py --yaml analysis_2017_config.yml --outdirtag HData_Nominal --send
python runAnalysis.py --yaml analysis_2018_config.yml --outdirtag HData_Nominal --send
