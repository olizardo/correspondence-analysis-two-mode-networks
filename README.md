# The Correspondence Analysis of Two-Mode Networks Revisited

This repository contains the complete replication materials—data, code, manuscripts, and presentation slides—for the paper **"The Correspondence Analysis of Two-Mode Networks Revisited"** by Omar Lizardo. 

The study re-examines the mathematical and conceptual linkages between several dual/two-mode (bipartite) network analysis methods, demonstrating that they are mathematically isomorphic or closely aligned when mapping dual relationships:
1. **Correspondence Analysis (CA)** (specifically under dual projection and Faust-style direct scaling)
2. **Bonacich's simultaneous group and individual scaling**
3. **The Method of Reflections** (developed by Hidalgo & Hausmann in the economic complexity literature)
4. **SimRank** (the iterative, graph-theoretic similarity algorithm)

The classic **Southern Women bipartite network** (18 women and 14 social events) is used as the empirical case study to illustrate these linkages.

---

## Repository Structure

The project directory is structured as follows:

```text
├── analysis.qmd                    # Main reproducible Quarto notebook (primary entry point)
├── ca.bib                          # Bibliography file for the LaTeX manuscript
├── correspondence-analysis-two-mode-networks.Rproj  # RStudio/Positron project file
├── main-final.tex                  # Complete, self-contained LaTeX manuscript
├── main-R1.tex                     # Revision 1 LaTeX manuscript
│
├── Functions/                      # Custom mathematical models and helper scripts
│   ├── reflections.R               # Iterative Method of Reflections (Hidalgo & Hausmann)
│   ├── SimRank.R                   # Iterative bipartite SimRank algorithm (Jeh & Widom)
│   └── ref.long.dat.R              # Reshapes reflections data to long format for ggbump
│
├── Plots/                          # Directory where all manuscript figures are saved
│
├── lit/                            # Reference literature and key foundation papers (PDFs)
│
├── Reviews and Response/           # Peer-review letters and journal response documents
│
├── slides/                         # Quarto (.qmd) and HTML presentation slide decks
│
└── old/                            # Archive of historical drafts, LaTeX fragments, and old scripts
```

---

## Getting Started & Prerequisites

To run the reproducibility workflow, you will need **R** and the **Quarto** CLI (pre-installed in Positron and RStudio).

### Required R Packages

Ensure you have the required R packages installed. You can install them by running the following command in your R console:

```r
install.packages(c(
  "conflicted", "cowplot", "ggcorrplot", "ggbump", "ggplot2", 
  "ggpubr", "ggrepel", "here", "factoextra", "igraph", 
  "kableExtra", "networkdata", "pals", "patchwork", "stringr"
))
```

---

## How to Reproduce the Findings

1. Clone or download this repository.
2. Open the R project file `correspondence-analysis-two-mode-networks.Rproj` in Positron or RStudio. This automatically sets the working directory to the project root.
3. Open `analysis.qmd` in the editor.
4. Render the notebook:
   - In **Positron**: Click the **Render** button in the upper-right of the editor, or use the Command Palette (`Ctrl+Shift+P` / `Cmd+Shift+P` and search `Quarto: Render`).
   - In **RStudio**: Click the **Render** button at the top of the file editor.
   - Or run in the terminal:
     ```bash
     quarto render analysis.qmd
     ```
5. Rendering compiles the document and regenerates all figures, saving them directly into the `Plots/` directory.

---

## License & Citation

If you use the materials or code in this repository, please cite the paper:

```bibtex
@article{lizardo2026correspondence,
  title={The Correspondence Analysis of Two-Mode Networks Revisited},
  author={Lizardo, Omar},
  year={2026},
  journal={Social Networks}
}
```
