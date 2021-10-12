.ONESHELL:
	conda activate r_env

main:
	r < 015_data_prep.R --no-save
	r < 019_predict_phenological_early.R --no-save
	r < 02_predict_phenological.R --no-save
predict1:
	r < 015_data_prep.R --no-save
	r < 021_produce_action_units.R --no-save
	r < 03_helpers.R --no-save
	r < 03_predict_data_driven.R --no-save
	r < 032_predict_data_driven_specgd.R --no-save
predict2:
	r < 035_predict_data_driven_hourly.R --no-save
earlytest:
	r < 01_hustle.R --no-save
