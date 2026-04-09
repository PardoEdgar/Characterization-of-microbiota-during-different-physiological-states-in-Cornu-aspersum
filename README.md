# Characterization of microbiota during different physiological states in *Cornu aspersum* (Coursework Project for undergraduate studies, 2024)

This repository contains the processed data and code associated with the research on gut microbiota characterization in the land snail *Cornu aspersum* across physiological states (active, aestivation, post-aestivation). In this brief study, we performed quality control with NanoStat for sequences obtained using Nanopore, taxonomic assignment with Kraken2, and statistical analyses in R (diversity, PCA, co-occurrence networks).

## Methods
- **Quality control**: NanoStat
- **Taxonomic classification**: Kraken2
- **Statistical analysis**: R (tidyverse, vegan, phyloseq, igraph) including alpha/beta diversity, PCA, and co-occurrence networks
- **Data visualization**: R (ggplot2)

## Physiological states
- Active
- Aestivation
- Post-aestivation

## Repository contents

| File | Description |
|------|-------------|
| `research_document.pdf` | Full research document |
| `Data_Bacteria_Reads.xlsx.csv` | Dataset used for the analyses (CSV format) |
| `microbiota_analysis.R` | R script containing the code for data processing, statistical analysis, and visualization |

## Requirements

- **Operating system**: Linux/Ubuntu
- **Sequence data**: FastQ files, Biom file, tsv file and csv table file.
- **Tools**: NanoStat, Kraken2
- **R version**: 4.4.2
- **R packages**: tidyverse, vegan, phyloseq, igraph

## Acknowledgments

I am grateful to Dr. Luz Helena Patiño for her invaluable support, lectures, and guidance throughout this research as part of my coursework. I also thank Universidad del Rosario for providing the resources, facilities, and equipment that made this work possible during my funded undergraduate studies. Special thanks go to Centro de Investigaciones Microbiológicas y Biotecnológicas (CIMBIUR) at Universidad del Rosario for their technical support.

## License

- **Code:** MIT License  
- **Research document, figures, and data:** Creative Commons Attribution 4.0 International (CC BY 4.0)

You are free to use, modify, and distribute the code under the MIT License.  
The research content (text, figures, and datasets) may be shared and adapted for any purpose, provided appropriate credit is given.

For more details:
- MIT License: https://opensource.org/licenses/MIT  
- CC BY 4.0: https://creativecommons.org/licenses/by/4.0/
