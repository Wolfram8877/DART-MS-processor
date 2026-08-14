# DART-MS Processor

An automated data processing pipeline for Direct Analysis in Real Time Mass Spectrometry (DART-MS) datasets. This tool streamlines the extraction, filtering, and visualization of raw mass spectrometry data, bridging the gap between high-throughput data acquisition and actionable analytical results.

https://github.com/user-attachments/assets/39a2c30a-32ca-4138-a4c5-1916de1447f8

## The Problem

While DART-MS allows for rapid, ambient ionization and high-throughput data acquisition—often analyzing samples in mere seconds—the subsequent data processing remains a severe bottleneck in analytical chemistry workflows. 

This is especially true in polymer additive analysis, where spectra are highly congested and complex:

<img width="645" height="316" alt="image" src="https://github.com/user-attachments/assets/35467add-29c8-4c8d-9787-1f59aeba81d4" />

Researchers relying on manual processing face several critical challenges:
* **The Throughput Mismatch:** The fundamental advantage of DART-MS (rapid analysis) is entirely negated by the hours required to manually process the resulting exports.
* **Data Complexity:** Raw mass spectrometry outputs require rigorous background subtraction, peak alignment, and thresholding. Performing this manually across dozens or hundreds of samples is tedious and subjective.
* **The Spreadsheet Bottleneck:** Relying on manual data extraction, repetitive copy-pasting, and formatting across multiple spreadsheets introduces a high risk of human error, compromises data integrity, and ultimately destroys analytical reproducibility.
* **Database Limitations for Polymers:** Unlike pesticide analysis, there is a severe lack of dedicated spectral databases for polymer additives. Manually querying general repositories (like PubChem) by chemical formula yields overwhelming and unreliable results due to structural isomers and the absence of polymer-specific context.

## The Solution

The DART-MS Processor replaces the manual data treatment bottleneck with a robust, programmatic workflow. By ingesting instrument exports, this pipeline automates the most time-consuming steps of mass spectrometry data analysis. It ensures that peak extraction and dataset alignment are handled consistently, guaranteeing strict analytical reproducibility.

To overcome the lack of specialized databases, the program automatically contextualizes findings by providing critical information for all candidate compounds: common polymer associations, potential contaminants, adduct formations, and industrial uses. It synthesizes this data to calculate a **confidence score**, reliably proposing the most probable candidate for each peak.

Finally, the processor features a multi-sample screening tool, allowing researchers to instantly visualize and map the distribution of specific compounds across large datasets.

## Core Capabilities

* **Batch Automation:** Processes large-scale sets of DART-MS data files simultaneously, completely eliminating file-by-file manual entry.
* **Data Parsing & Peak Filtering:** Cleans raw output, automatically handling noise reduction and peak extraction based on user-defined parameters.
* **Isotopic Analysis & Intelligent Scoring:** Evaluates mass accuracy (ppm error), isotopic profiles, adduct formations, and polymer associations to assign a robust confidence score to candidate molecules, removing the guesswork from identification.
* **Automated Structural Visualization:** The tool automatically fetches and displays the chemical structure, CAS registry number, and known industrial uses for each candidate, providing immediate visual and chemical context.
* **Cross-Sample Screening:** Aggregates results from multiple analyses to easily track and map the presence of specific compounds across various samples.
* **Excel-Ready Export:** Outputs clean, formatted, and analysis-ready data directly into Excel spreadsheets, eliminating tedious manual copy-pasting and facilitating immediate review.
