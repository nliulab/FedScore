FedScore: A privacy-preserving framework for federated scoring system development
=========================

R and Python code for our proposed FedScore framework for generation of scoring systems across multiple sites in a privacy preserving way. See our new [Preprint](https://arxiv.org/abs/2303.00282) for the whole story.

## Table of contents
* [Introduction](#introduction)
* [Usage](#usage)
* [Citation](#citation)
* [Contact](#contact)

## Introduction

Cross-institutional collaboration has gained popularity in recent years as a way to accelerate medical research and facilitate quality improvement. Federated learning (FL), sometimes referred to as distributed learning or distributed algorithms, can avoid data sharing by collectively training algorithms without exchanging patient-level data. There exist many applications of FL for medical image data, most of which use black box models from computer vision. Interpretable models, on the contrary, have fewer instances of FL applications despite their popularity in clinical research.

As a type of interpretable risk scoring model, scoring systems have been employed in practically every diagnostic area of medicine. Regardless of diverse development strategies, scoring systems have usually been created using singlesource data, limiting application at other sites if the development data has insufficient sample size or is not representative. Although it is possible to develop scoring systems on pooled data the process of doing such pooling, as noted previously, is time consuming and difficult to achieve due to privacy restrictions. As a result, frameworks for building scoring systems in a
federated manner are needed to overcome such difficulties. 

To fill this gap, we propose FedScore, a first-of-its-kind framework for building federated scoring systems across multiple sites. 

<div class="figure" style="text-align: center">

<img src="figures/Figure1.pdf" width="70%"/>

</div>

## Usage

The structure of this repository is detailed as follows:

- `Benchmark_scripts/...` contains the scripts for benchmark dataset generation (master_data.csv).
- `Benchmark_scripts/...` contains the scripts for building the various task-specific benchmark models.
-  


## Citation

Li, S., Ning, Y., Ong, M. E. H., Chakraborty, B., Hong, C., Xie, F., ... & Liu, N. (2023). FedScore: A privacy-preserving framework for federated scoring system development. arXiv preprint arXiv:2303.00282.

## Contact

- Siqi Li(Email: <siqili@u.duke.nus.edu>)
- Nan Liu (Email: <liu.nan@duke-nus.edu.sg>)
