# Files

Folder structure
- `data/`: raw data (measured values)
- `data_aux/`: parameters, constants, configuration files
- `data_proc/`: intermediate files generated during data processing
- `data_ref/`: reference data
- `results/`: analysis results

# Experiments

Folder names correspond to date of data generation (i.e., sequencing).

<!--- Markdown table formatting notes: left align text; use underscores to indicate minimum column width (especially for GitHub display) -->
| Folder name<br/>__________ | Notebook<br/>__________ | Experiment name (Benchling)<br/>_________________________ | Description <br/>_____________________________________________________ | Sequencing<br/>________________ |
| :----------- | :-------- | :--------------------------- | :----------- | :--------------- |
| 20230717    | `20230717 HEK scBarcode.ipynb` | [2023-06-28 Split-Pool HEK Nuclei](https://benchling.com/s/etr-4FqDoTQiWpBfQlyOBYhz) | Sequencing of barcodes from serial and limiting dilutions of SPRITE Zero (concentration-doubling)-barcoded HEK nuclei | NextSeq 77x77 |
| 20230831 (or 20230831_barcodes) | `20230830 In vitro barcoding Nanopore.ipynb` | [2023-08-09 DPM ligation to split-pooled oligo](https://benchling.com/s/etr-XcdY7Za2GoVwcRjFEqBF) | Nanopore sequencing of gel-cut ~197 bp and ~250 bp bands of amplified, unblocked barcoded oligo (Oligo + Odd + ER + dA + DPM + Odd + Y) | Nanopore Flongle |
| 20230831_scBarcode | `20230831_scBarcode.ipynb` | [2023-08-10 HEK H3K4me3 scChIP-seq](https://benchling.com/s/etr-V1QGaZkTpBrchx0XtlvM) | Sequencing of barcodes from serial dilutions of FACS-sorted singlets and clumps of SPRITE Zero (concentration-doubling)-barcoded HEK nuclei | NextSeq 101x201 |
| 20230831_HEK_H3K4me3_scChIPseq | `20230831_HEK_H3K4me3_scChIPseq.ipynb` | [2023-08-10 HEK H3K4me3 scChIP-seq](https://benchling.com/s/etr-V1QGaZkTpBrchx0XtlvM) | Sequencing of genomic DNA and barcodes from H3K4me3 ChIP of 1500 flow-sorted singlet SPRITE Zero (concentration-doubling)-barcoded HEK nuclei | NextSeq 101x201 |
| 20231002 | `20231002.ipynb` | [2023-09-25 Single Cell Barcode Troubleshooting](https://benchling.com/s/etr-fV6EV9txrU2wIEb8MfxX) | Sequencing of barcodes from individual flow-sorted SPRITE Zero (concentration-doubling)-barcoded HEK nuclei, using old SPRITE Zero plates | NextSeq 151x151 |
| 20231017 | `20231017.ipynb` | [2023-10-09 Single Cell Barcode Troubleshooting, v2](https://benchling.com/s/etr-5fnlp2r3TvdgBXmWeyxQ) | Sequencing of barcodes from individual flow-sorted SPRITE Zero (concentration-doubling)-barcoded HEK nuclei, barcoded using fresh SPRITE Zero tag plates | NextSeq 51x51 |
| 20231107 | `20231107.ipynb` | [2023-11-02 Single Cell Barcode Troubleshooting, v3](https://benchling.com/s/etr-55XImqxPMe9Y2BOH1EEd) | Sequencing of barcodes from individual flow-sorted HEK nuclei, using high concentration of barcodes with wash steps in between rounds | AVITI 100x200 |
| 20231208 | `20231208.ipynb` | [2023-11-29 Single Cell Barcode Troubleshooting, v4 (terminal tag vs. EDTA quench)](https://benchling.com/s/etr-3bKIM8CScL814XrVEUBG) | Sequencing of barcodes from SPRITE Zero-barcoded HEK nuclei, comparing EDTA quench vs. terminal tag ligation between each round of tag ligation | AVITI 120x180 |
| 20231230 | `20231230.ipynb` | [2023-12-22 Tag Plate Contamination Test](https://benchling.com/s/etr-Umk5xAsODA8uMRHekEBi) | Test contamination of my SPRITE Zero and NYLigOdd tag plates | AVITI 150x150 |
| 20240112 | `20240112.ipynb` | [2024-01-05 Tag Plate Contamination Test, v2](https://benchling.com/s/etr-78hBQk3CEBX5dD6ojPmr) | Test contamination of Andrew Perez's SPRITE Zero and NYLigOdd tag plates | NextSeq 51x51 |
| 20240124 | `20240124.ipynb` | [2024-01-18 Tag Plate Contamination, v3](https://benchling.com/s/etr-0kEaLIIMFTjHJbNXzjpb) | Test contamination of stock, unannealed SPRTIE Zero R1-R4 tag plates. (Anneal new SPRITE Zero R1-R4 and NYLigOdd tag plates, then check for contamination.) | AVITI 120x180 |