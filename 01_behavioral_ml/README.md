# 01_behavioral_ml

Behavioral analysis of MoSeq and traditional behavioral assays in CS and WMS stroke mouse models.

## Structure

### `analyses/`
Core analysis notebooks for the main CS/WMS cohort. See [`analyses/README.md`](analyses/README.md) for details.

| Notebook | Description |
|----------|-------------|
| `Classifier.ipynb` | Logistic regression classifiers on MoSeq syllable usage to distinguish treatment groups. This contains these analyses: syllable cluster classification, fingerprint (syllable + scalar) classification, and scalars vs. MoSeq comparison. |
| `Correlate-MoSeq-And-Traditional-Behavior.ipynb` | Spearman correlations between syllable usage and six behavioral assay scores (rotarod, wire hang, gridwalk, Y-maze, NORT, open field) with FDR correction. Outputs a syllable × metric heatmap and cluster-level correlation matrix. |
| `Generate-Syllable-Heatmap-CS-WMS.ipynb` | Heatmap of 7DPI/BSL syllable usage ratios in young males, comparing CS and WMS groups across all 38 syllables. |
| `Scalar-Comparison.ipynb` | Bar plots of mean velocity and total distance traveled per animal at BSL, 7DPI, and 30DPI, stratified by treatment and sex. |
| `Stroke_Size_Quantification_Plot.ipynb` | Bar plot of histological stroke injury volume (mm^3) across CS/WMS × age × sex groups. |
| `Standard-Behavioral-Assays.ipynb` | Plots and Mann-Whitney comparisons for six standard assays: rotarod, wire hanging, gridwalk, Y-maze, NORT, and open field. |
| `Transition_Probabilities.ipynb` | Per-animal 38×38 syllable transition matrices. Compares probability of transitioning into motor-deficit (LVHD) syllables in CS mice and cognitive-deficit (LVLD) syllables in WMS mice, split by age and sex. |

---

### `nort/`
SVM-based classifier for automated NORT (Novel Object Recognition Test) scoring from DeepLabCut pose data.

| File | Description |
|------|-------------|
| `Novel-Object-Recognition-Test.ipynb` | Training pipeline — extracts geometric features from DLC pose estimates and trains an RBF-SVM to detect object exploration bouts |
| `Apply-NORT-Classifier.ipynb` | Inference pipeline — applies the trained SVM to novel videos and computes discrimination index per animal |
| `GUI.py` | Helper GUI for manual annotation review |

For convenience purposes, these functions are bundled as a Python package and are located [here](https://github.com/llorente-lab/nort).

---

### `nsg_vs_bl6/`
Companion analyses comparing NSG (immunodeficient) vs. C57BL/6 mice at baseline to characterize strain-level behavioral differences.

| Notebook | Description |
|----------|-------------|
| `Clean-Dataframe-And-Plot-Cognitive-Motor-Clusters.ipynb` | Syllable clustering and behavioral group comparisons for NSG vs. C57BL/6 |
| `Transition_probabilities_NSG_bsl_vs_C57_bsl_CS7dpi.ipynb` | Syllable transition probability comparisons across strains and stroke conditions |