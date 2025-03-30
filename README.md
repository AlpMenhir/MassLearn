# MassLearn

**MassLearn** is a lightweight tool for untargeted metabolomics data processing, designed to streamline feature extraction, noise removal, and statistical analysis.

---

## 🚀 Overview

MassLearn is built to work with raw `.raw` files from **Waters Xevo G2-XS qToF** instruments. It supports preprocessing, feature grouping, and statistical analysis, making it easier to explore and visualize complex metabolomics datasets.

---

## 🔧 Pipeline Summary

1. **Input Files:**
   - `.raw` files from Waters instruments
   - Conversion to open format using **MSConvert**
   - Feature detection using **MZmine**

2. **Noise Filtering:**
   - Removes noise traces before MZmine processing

3. **Feature Grouping:**
   - Groups features across samples using MS/MS data and retention time
   - Features from the same precursor ion are grouped together

4. **Project File Generation:**
   - Internally creates a project structure for organized data handling

5. **Statistical Analysis:**
   - Automated PCA and volcano plot generation
   - Parametric and non-parametric tests on feature group abundance

---

## 📊 Output

- Cleaned feature lists
- Grouped feature tables
- PCA plots
- Volcano plots
- Statistical test results

---

## 📎 Requirements

- MZmine
- MSConvert (from ProteoWizard)
- Python 3.x with required packages (see `requirements.txt`)

---

## 📁 Folder Structure (Example)

