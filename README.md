# DART-MS Processor

An automated data processing pipeline for Direct Analysis in Real Time Mass Spectrometry (DART-MS) datasets. This tool streamlines the extraction, filtering, and visualization of raw mass spectrometry data, bridging the gap between high-throughput data acquisition and actionable analytical results.

https://github.com/user-attachments/assets/39a2c30a-32ca-4138-a4c5-1916de1447f8

## The Problem

While DART-MS allows for rapid, ambient ionization and high-throughput data acquisition—often analyzing samples in mere seconds—the subsequent data processing remains a severe bottleneck in analytical chemistry workflows. 

Especially with polymer additive analysis where the spectrums are fat : 

<img width="1290" height="633" alt="image" src="https://github.com/user-attachments/assets/35467add-29c8-4c8d-9787-1f59aeba81d4" />


The manual processing face several critical challenges:
* **The Throughput Mismatch:** The fundamental advantage of DART-MS (rapid analysis) is entirely negated by the hours required to manually process the resulting exports.
* **Data Complexity:** Raw mass spectrometry outputs require rigorous background subtraction, peak alignment, and thresholding. Performing this manually across dozens or hundreds of samples is tedious and subjective.
* **The Spreadsheet Bottleneck:** Relying on manual data extraction, repetitive copy-pasting, and formatting across multiple spreadsheets introduces a high risk of human error, compromises data integrity, and ultimately destroys analytical reproducibility.
* **The Lack of plastic database for DART-MS:** Unlike pesticides plastic really lack of specific database for polymer additives and searching manually on general database like Pubchem with the right chemical formula is unreliable cause to the many options and the lack of information concerning each compounds related to polymer.


## The Solution

The DART-MS Processor replaces the manual data treatment bottleneck with a robust, programmatic workflow. By ingesting instrument exports, this pipeline automates the most time-consuming steps of mass spectrometry data analysis. It ensures that peak extraction and dataset alignment are handled consistently every single time. This guarantees strict analytical reproducibility and frees up researchers to focus on data interpretation rather than data formatting. The program gives all the useful information to conclude the analysis : common polymer presence, possible contaminations, adducts information and industrial uses of all candidates. To sum all the given informations, the program gives a score to propose the most probable candidate to the peak.

The 2nd use of the processor is to make a screening with the results of many samples to see easily which compounds are present in which samples.

## Core Capabilities

* **Batch Automation:** Processes large-scale sets of DART-MS data files simultaneously, completely eliminating file-by-file manual entry.
* **Data Parsing & Normalization:** Cleans and structures raw output into standardized formats, automatically handling noise reduction and peak filtering based on user-defined parameters.
* **Automated Visualization:** Generates instant graphical representations of mass spectra to facilitate rapid quality control and dataset comparison.
* **Standardized Export:** Outputs clean, analysis-ready data into universally accepted formats (CSV, Excel) for seamless integration into downstream statistical software or laboratory information management systems (LIMS).
