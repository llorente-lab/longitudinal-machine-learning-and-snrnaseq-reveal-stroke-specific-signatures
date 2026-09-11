# Welcome to the official repository for the manuscript "Longitudinal Machine Learning and snRNA-Seq Reveal Stroke-Specific Behavioral and Transcriptional Signatures in a Mouse Model of Stroke".

## 🔬 Overview & Scientific Rationale

Ischemic stroke is a leading cause of death and long-term disability worldwide, accounting for over 690,000 new cases each year. The clinical presentation and underlying pathology diverge sharply between stroke subtypes:
*   **Cortical Stroke (CS):** Patients typically present with prominent acute sensorimotor impairments, such as hemiparesis and speech deficits.
*   **White Matter Stroke (WMS):** Represents a highly prevalent subtype that often manifests as subtle, progressively worsening cognitive deficits rather than overt motor paralysis.

### The Preclinical Bottleneck
Traditional rodent behavioral assessments are constrained, low-dimensional, and rely on subjective scoring metrics. Consequently, they lack the comprehensive, unbiased resolution required to capture these nuanced, naturalistic motor and cognitive functional deficits over time.

### Our Approach
To overcome these limitations, this study integrates longitudinal machine learning, single-nucleus transcriptomics, and high-resolution microglial morphometry to bridge the gap between behavioral phenotypes, immune-mediated molecular programs, and localized structural remodeling in the brain:

1. **Unbiased Behavioral Tracking (MoSeq):** We implement **Motion Sequencing (MoSeq)**, an unsupervised machine-learning pipeline analyzing 3D depth recordings of freely moving mice to decode sub-second behavioral motifs ("syllables") longitudinally.
2. **Stroke Subtype-Specific Trajectories:** 
   * **CS mice** exhibit motor impairments that partially resolve during the chronic phase.
   * **WMS mice** develop progressively worsening cognitive-like deficits. 
   * These behavioral deficits are consistent across age and sex, but are **most pronounced in aged, female mice**.
3. **Parallel Molecular Signatures (snRNA-seq):** Single-nucleus transcriptomic profiling of ischemic tissue reveals parallel cellular cascades:
   * **CS tissue** undergoes a dynamic temporal shift, transitioning from a highly active inflammatory state to a homeostatic gene expression profile.
   * **WMS tissue** exhibits a chronic and sustained downregulation of essential myelination and synaptic organization pathways.
4. **Immunological & Morphological Causal Mapping:** We leverage comparative studies between immunocompetent (C57BL/6J) and immunodeficient (NSG) mice to establish a causal link between neuroinflammation, microglial activation, and behavioral dysfunction:
   * **Microglial Morphometry & Clustering** processes over 50 shape, skeleton, and graph-theoretical features to distinguish cellular activation states. PCA reveals that immunocompetent C57 microglia maintain a highly ramified surveillance state, whereas immunodeficient NSG microglia display a distinct hyper-ramified morphology—an intermediate state of activation marked by smaller volumes, lower betweenness centrality, and shorter branch skeletons.
   * **Behavioral Subspace Coupling** demonstrates that distinct immunological baseline states occupy unique PC behavioral subspaces. MoSeq successfully captures these downstream behavioral changes (increased motor-related and decreased cognitive-related syllable usage in NSG mice), which are further exacerbated by stroke injury.

By mapping these unique, subtype-specific behavioral trajectories directly to underlying cellular morphometrics and transcriptional programs, this study provides a crucial scientific foundation for the development of subtype-tailored stroke therapies.

## 📂 Repository Structure
To facilitate reproducibility, this repository is organized into dedicated, self-contained directories for each distinct phase of the analysis:

* `01_behavioral_ml/`: Contains the Motion Sequencing (MoSeq) pipeline and machine learning code to process raw 3D depth video data, train the unsupervised autoregressive hidden Markov model (AR-HMM), extract sub-second motor and cognitive behavioral syllables, and run supervised classifiers to distinguish stroke subtypes (CS vs. WMS) and recovery stages across age and sex.
* `02_snrna/`: Contains the single-nucleus RNA-sequencing (snRNA-seq) analysis scripts. It handles quality control, integration across cohorts, cell type identification (excitatory/inhibitory neurons, astrocytes, oligodendrocytes, OPCs, microglia/macrophages, and endothelial cells), differential expression testing between acute (7 DPI) and chronic (30 DPI) phases, pathway enrichment, and CellChat-based ligand-receptor interactome modeling.
* `03_microglia_morphology`: Houses the morphological and skeleton-based analysis scripts for microglia. It processes cell body and arbor metrics (e.g., cell volume, node count, betweenness centrality, skeleton length), conducts PCA to isolate key structural components, and computes comparative statistics between immunocompetent (C57BL/6J) and immunodeficient (NSG) mice. 
