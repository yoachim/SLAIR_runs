# SLAIR_runs
LSST observing cadence experiments


| Directory name | notes |
|------ | ------ |
| baseline | A basline-ish run. The 4 anounced DDFs, regions selected in large blocks, pairs in g,r,i. |
| baseline_mix | Similar to baseline, pairs in g+r, r+i, i+g. |
| rolling | Rolling cadence splitting WFD region into North and South equal areas, Pairs in g,r,i.|
| roll_mix | Rolling cadence on WFD area, north south split 75%/25%. Pairs in different filters. |
| cadence_mix | Add a basis function that drives 3-day cadence in g,r,i in WFD region. Paris in different filters. |
| cadence_roll_75_mix | use the cadence basis function and rolling cadence together |
| roll_mix_100 | Use cadence basis function and rolling cadence, with the roll being a full on/off |
| cadence_mix_wfd | same as cadence_mix, but with only the WFD regions |
| tight_mask | use a restrictive alt-az mask to force merdian scanning, try to force cadence (didn't work well) |
| tight_mask_simple | Use a tight alt-az mask, only y in twilight. No 5-sigma depth used for filter selection |
| tms_drive | tight mask, only y in twilight, no 5-sigma depth used for filter, added basis function to reward 2.1-5 day cadence in g,r,i |
| tms_roll | Like tight_mask_simple, adding rolling cadence | 
| year_1 | Work on a survey that does a good job in year 1 closing sky and gathering templates |

Results at:  https://lsst-web.ncsa.illinois.edu/sim-data/beta_slair_surveys/