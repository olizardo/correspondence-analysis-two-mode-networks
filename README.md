# The Correspondence Analysis of Two-Mode Networks Revisited

This repository contains the complete replication materials—data, code, manuscripts, and presentation slides—for the paper **"The Correspondence Analysis of Two-Mode Networks Revisited"** by Omar Lizardo.

The study re-examines the mathematical and conceptual linkages between several dual/two-mode (bipartite) network analysis methods, demonstrating that they are mathematically isomorphic or closely aligned when mapping dual relationships:
1. **Correspondence Analysis (CA)** (specifically under dual projection and Faust-style direct scaling)
2. **Bonacichs simultaneous group and individual scaling**
3. **The Method of Reflections** (developed by Hidalgo & Hausmann in the economic complexity literature)
4. **SimRank** (the iterative, graph-theoretic similarity algorithm)

The classic **Southern Women bipartite network** (18 women and 14 social events) is used as the empirical case study to illustrate these linkages.

---

## 📂 Repository Structure

The project directory is structured as follows:

```text
├── analysis.qmd                    # Main reproducible Quarto notebook (primary entry point)
├── ca.bib                          # Bibliography file for the LaTeX manuscript
├── correspondence-analysis-two-mode-networks.Rproj  # RStudio/Positron project file
├── main-final.tex                  # Complete, self-contained LaTeX manuscript
├── main-R1.tex                     # Revision 1 LaTeX manuscript
├── Functions/                      # Custom mathematical models and helper scripts
│   ├── reflections.R               # Iterative Method of Reflections (Hidalgo & Hausmann)
│   ├── SimRank.R                   # Iterative bipartite SimRank algorithm (Jeh & Widom)
│   └── ref.long.dat.R              # Reshapes reflections data to long format for ggbump
├── Plots/                          # Directory where all manuscript figures are saved
├── lit/                            # Reference literature and key foundation papers (PDFs)
├── Reviews and Response/           # Peer-review letters and journal response documents
├── slides/                         # Quarto (.qmd) and HTML presentation slide decks
└── old/                            # Archive of historical drafts, LaTeX fragments, and old scripts
```

---

## 🛠️ Prerequisites & Installation

To run the reproducibility workflow, you will need **R** and the **Quarto** CLI (pre-installed in Positron and RStudio).

### 1. Required R Packages
Ensure you have the required R packages installed. You can install them by running the following command in your R console:

```R
install.packages(c(
  "conflicted", "cowplot", "ggcorrplot", "ggbump", "ggplot2", 
  "ggpubr", "ggrepel", "here", "factoextra", "igraph", 
  "kableExtra", "networkdata", "pals", "patchwork", "stringr"
))
```



---

## 🚀 How to Reproduce the Findings

1. **Clone the Repository**: Clone this repository to your local machine using git or download it as a ZIP file.
2. **Open the Project**: Open the `.Rproj` or `.R` file in your editor (e.g., RStudio or Positron). This ensures paths are resolved correctly relative to the project root.
3. **Install Dependencies**: Ensure the packages listed above are installed.
4. **Run the Computational Pipeline**:
   * **Using the Command Line (Quarto CLI)**:
     ```bash
     quarto render analysis.qmd
     ```
   * **Using R**:
     ```R
     quarto::quarto_render("analysis.qmd")
     ```
   * Or run interactively in your IDE. This step ensures that all tables and figures are updated directly from the code, guaranteeing that the numbers in the paper are exactly what the code computes.

---

## 📝 License & Citation

If you use the materials or code in this repository, please cite the paper:

```bibtex
@article{lizardo2025correspondence,
  title={The Correspondence Analysis of two-mode networks revisited},
  author={Lizardo, Omar},
  journal={Social Networks},
  volume={83},
  pages={134--151},
  year={2025},
  publisher={Elsevier}
}
```
