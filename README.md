# CellCraft Plugins

This repository contains plugin configurations and information for [CellCraft](https://github.com/cxinsys/cellcraft), a web-based visual programming application for gene regulatory network (GRN) inference.

## Overview

CellCraft integrates multiple GRN reconstruction tools through its modular plugin system. This repository serves as a central hub for managing and documenting all available plugins that can be used within the CellCraft environment.

## Available Plugins

| Plugin | Description | Paper | GitHub |
|--------|-------------|-------|--------|
| **TENET** | Transfer Entropy-based Network Reconstruction - Employs transfer entropy from information theory to infer causal relationships in gene regulatory networks from single-cell RNA-seq data | [Kim et al., NAR 2021](https://academic.oup.com/nar/article/49/1/e1/5973444?login=false) | [neocaleb/TENET](https://github.com/neocaleb/TENET) |
| **SCODE** | Efficient Regulatory Network Inference During Differentiation - Infers regulatory networks from single-cell RNA-seq data during cell differentiation using ordinary differential equations | [Matsumoto et al., Bioinformatics 2017](https://doi.org/10.1093/bioinformatics/btx194) | [hmatsu1226/SCODE](https://github.com/hmatsu1226/SCODE) |
| **SCRIBE** | Causal Network Inference Using Single-Cell Expression Dynamics - Uses restricted directed information to infer causal regulatory networks from coupled single-cell expression dynamics | [Qiu et al., Cell Systems 2020](https://www.sciencedirect.com/science/article/pii/S2405471220300363?via%3Dihub) | [cole-trapnell-lab/Scribe](https://github.com/cole-trapnell-lab/Scribe) |
| **LEAP** | Pseudotime-based Co-expression Network Construction - Constructs gene co-expression networks for single-cell RNA-seq data using pseudotime ordering | [Specht & Li, Bioinformatics 2016](https://doi.org/10.1093/bioinformatics/btw729) | [cran/LEAP](https://github.com/cran/LEAP) |
| **GRNBoost2** | Scalable Gene Regulatory Network Inference - Provides efficient and scalable inference of gene regulatory networks using gradient boosting | [Moerman et al., Bioinformatics 2018](https://doi.org/10.1093/bioinformatics/bty916) | [aertslab/GRNBoost](https://github.com/aertslab/GRNBoost) |
| **GENIE3** | Tree-based Network Inference - Infers regulatory networks from expression data using tree-based methods (Random Forests or Extra-Trees) | [Huynh-Thu et al., PLOS ONE 2010](https://doi.org/10.1371/journal.pone.0012776) | [vahuynh/GENIE3](https://github.com/vahuynh/GENIE3) |

## Contributing

If you would like to contribute a new plugin or improve existing ones, please:

1. Fork this repository
2. Create a feature branch
3. Submit a pull request with a clear description of your changes
