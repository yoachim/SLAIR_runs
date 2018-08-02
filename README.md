# SLAIR_runs
LSST observing cadence experiments


| Directory name | notes |
|------ | ------ |
| baseline | A basline-ish run. The 4 anounced DDFs, regions selected in large blocks, pairs in g,r,i. |
| baseline_mix | Similar to baseline, pairs in g+r, r+i, i+g. |
| rolling | Rolling cadence splitting WFD region into North and South equal areas, Pairs in g,r,i.|
| roll_mix | Rolling cadence on WFD area, north south split 75%/25%. Pairs in different filters. |
| cadence_mix | Add a basis fucntion that drives 3-day cadence in g,r,i in WFD region. Paris in different filters. |
| cadence_roll_75_mix | use the cadence basis function and rolling cadence together |
