1_Planck_Pre_Launch.camb:
	-@bash $(SCRIPT_DIR)/run_params.sh $(PARAMETER_DIR)/1_Planck_Pre_Launch.camb.ini $(PARAMETER_ANALYSIS_DIR)/1_Planck_Pre_Launch.camb.ini
2_DESGC_Planck.camb:
	-@bash $(SCRIPT_DIR)/run_params.sh $(PARAMETER_DIR)/2_DESGC_Planck.camb.ini $(PARAMETER_ANALYSIS_DIR)/2_DESGC_Planck.camb.ini
all_targets: 1_Planck_Pre_Launch.camb 2_DESGC_Planck.camb 