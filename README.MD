# mGeneChart 📊🧫

## Overview 

This pipeline interprets results from ABRicate and mobileelementfinder CSV and .tab/.tabular outputs to analyze the colocation of antimicrobial resistance (AMR) genes and mobile genetic elements (MGEs). It is tailored for these two output formats but can be adapted for others.

## Purpose

The script identifies proximal instances between AMR genes and mobile elements, creating a numerical representation of the risk of gene mobility. While plasmid contigs are assumed mobile, the pipeline focuses on predicting the mobility of AMR genes within chromosomal structures by observing their proximity to MGEs.

![Results of E.coli](Ecoli_Results_Sample.png)


## Features

- **Colocation Analysis:** Parses and cross-checks AMR gene indices with IS segments to assess mobility risk.
- **Visualization:** Returns results using MatPlotLib, with Biopython formatting for prokaryotic genes.
- **Transposase Identification:** Maps transposase genes using TE identifiers, leveraging existing labels from IS and MGE Finder.
- **Extensible:** Designed to incorporate additional factors such as TEs, ICEs, and integrons into mobility prediction.

## Resources & Imports

- **ABRicate** (CARD database) for AMR gene identification
- **mobileelementfinder** for locating insertion sequences (IS)
- **Integron Finder** (Galaxy local)
- **MatPlotLib** and **Biopython** for visualization

## Glossary

- **IS:** Insertion Sequence
- **cn:** Composite Transposon
- **Tn:** (Unit) Transposon
- **MITE:** Miniature Inverted Repeat
- **ICE:** Integrative Conjugative Element
- **CIME:** Cis-mobilizable Element
- **IME:** Integrative Mobilizable Element

## Future Directions

- Integration of additional TE screeners
- Enhanced mapping of transposase proximity
- Support for more databases and output formats

---
