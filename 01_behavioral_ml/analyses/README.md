# analyses

Jupyter notebooks for behavioral analysis of MoSeq motion segmentation data in CS and WMS stroke mouse models.

## Notebooks

| Notebook | Description |
|----------|-------------|
| `Classifier.ipynb` | Logistic regression classifiers on MoSeq syllable usage to distinguish treatment groups. This contains these analyses: syllable cluster classification, fingerprint (syllable + scalar) classification, and scalars vs. MoSeq comparison. |
| `Correlation.ipynb` | Spearman correlations between syllable usage and six behavioral assay scores (rotarod, wire hang, gridwalk, Y-maze, NORT, open field) with FDR correction. Outputs a syllable × metric heatmap and cluster-level correlation matrix. |
| `Generate-Syllable-Heatmap-CS-WMS.ipynb` | Heatmap of 7DPI/BSL syllable usage ratios in young males, comparing CS and WMS groups across all 38 syllables. |
| `Scalar-Plots.ipynb` | Bar plots of mean velocity and total distance traveled per animal at BSL, 7DPI, and 30DPI, stratified by treatment and sex. |
| `Stroke_Size_Quantification_Plot.ipynb` | Bar plot of histological stroke injury volume (mm^3) across CS/WMS × age × sex groups. |
| `Traditional-Behavior-Tests.ipynb` | Plots and Mann-Whitney comparisons for six standard assays: rotarod, wire hanging, gridwalk, Y-maze, NORT, and open field. |
| `Transition_Probabilities.ipynb` | Per-animal 38×38 syllable transition matrices. Compares probability of transitioning into motor-deficit (LVHD) syllables in CS mice and cognitive-deficit (LVLD) syllables in WMS mice, split by age and sex. |