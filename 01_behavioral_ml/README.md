# 01_behavioral_ml

Behavioral analysis of MoSeq motion segmentation data and traditional behavioral assays in CS and WMS stroke mouse models.

## Structure

### `analyses/`
Core analysis notebooks for the main CS/WMS cohort. See [`analyses/README.md`](analyses/README.md) for details.

| Notebook | Description |
|----------|-------------|
| `Classifier.ipynb` | Logistic regression classifiers on syllable usage, fingerprint features, and scalar summaries |
| `Correlation.ipynb` | Spearman correlations between syllable usage and behavioral assay scores |
| `Generate-Syllable-Heatmap-CS-WMS.ipynb` | Heatmap of 7DPI/BSL syllable usage ratios (CS vs. WMS, young males) |
| `Scalar-Plots.ipynb` | Velocity and distance traveled plots across timepoints |
| `Stroke_Size_Quantification_Plot.ipynb` | Stroke injury volume bar plot by group |
| `Traditional-Behavior-Tests.ipynb` | Rotarod, wire hanging, gridwalk, Y-maze, NORT, open field |
| `Transition_Probabilities.ipynb` | Syllable transition matrices; motor and cognitive deficit transition probabilities |

---

### `nort/`
SVM-based classifier for automated NORT (Novel Object Recognition Test) scoring from DeepLabCut pose data.

| File | Description |
|------|-------------|
| `nort.ipynb` | Training pipeline — extracts geometric features from DLC pose estimates and trains an RBF-SVM to detect object exploration bouts |
| `Run_Results.ipynb` | Inference pipeline — applies the trained SVM to novel videos and computes discrimination index per animal |
| `gui.py` | Helper GUI for manual annotation review |

---

### `nsg_vs_bl6/`
Companion analyses comparing NSG (immunodeficient) vs. C57BL/6 mice at baseline to characterize strain-level behavioral differences.

| Notebook | Description |
|----------|-------------|
| `Clean-Dataframe-And-Plot-Cognitive-Motor-Clusters-2.ipynb` | Syllable clustering and behavioral group comparisons for NSG vs. C57BL/6 |
| `Transition_probabilities_NSG_bsl_vs_C57_bsl_CS7dpi-2.ipynb` | Syllable transition probability comparisons across strains and stroke conditions |
