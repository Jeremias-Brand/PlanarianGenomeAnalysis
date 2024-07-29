# Conservation scoring

Details on how to prepare the liftover are in `conserved_elements`.
Once the liftover is done use `conserved_elements/granges_liftover_full_peak_size.R` to score conservation of each *Schmidtea mediterranea* ATAC-seq peak.
Summarize the results using `conserved_elements/summarize_liftover_full.R`.

Finally, plot the results and perform significance testing using `scripts/plot_and_test_conservation.R`.

