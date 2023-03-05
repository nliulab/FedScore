FedScore: A privacy-preserving framework for federated scoring system development
=========================

R and Python code for our proposed FedScore framework for generation of scoring systems across multiple sites in a privacy preserving way. See our new [Preprint](https://arxiv.org/abs/2303.00282) for the whole story.

## Table of contents
* [Introduction](#introduction)
* [Usage](#usage)
* [Citation](#citation)
* [Contact](#contact)

## Introduction

Clinical decisions in the emergency department play an important role in optimizing urgent patient care and scarce resources. And unsurprisingly, machine learning based clinical prediction models have been widely adopted in the field of emergency medicine.

In parellel to the rise of clinical prediction models, there has also been a rapid increase in adoption of Electronic Health records (EHR) for patient data. The Medical Information Mart for Intensive Care ([MIMIC)-IV]((https://physionet.org/content/mimiciv/1.0/)) and [MIMIC-IV-ED](https://physionet.org/content/mimic-iv-ed/1.0/) are examples of EHR databases that contain a vast amount of patient information.

There is therefore a need for publicly available benchmark datasets and models that allow researchers to produce comparable and reproducible results. 

For the previous iteration of the MIMIC database (MIMIC-III), several benchmark pipelines have published in [2019](https://github.com/YerevaNN/mimic3-benchmarks) and [2020](https://github.com/MLforHealth/MIMIC_Extract).

Here, we present a workflow that generates a benchmark dataset from the MIMIC-IV-ED database and constructs benchmark models for three ED-based prediction tasks.


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
