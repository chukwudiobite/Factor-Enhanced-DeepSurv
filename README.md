# Factor-Enhanced-DeepSurv

📌 Project Overview

This repository contains the code and materials for the research paper:

"Factor enhanced DeepSurv: A deep learning approach for predicting survival probabilities in cirrhosis data"Published in Computers in Biology and Medicine (2025)Authors: Chukwudi Paul Obite, Emmanuella Onyinyechi Chukwudi, Merit Uchechukwu, Ugochinyere Ihuoma Nwosu (https://www.sciencedirect.com/science/article/abs/pii/S0010482525003142)

🎯 Objective

To develop and validate a novel deep learning survival model called Factor Enhanced DeepSurv (FE-DeepSurv) that integrates factor analysis and deep neural networks to predict survival probabilities in cirrhosis patients.

📁 Dataset

Source: Kaggle Cirrhosis Dataset

Observations: 276 patients

Variables: 17 predictors + survival time + censoring indicator

🧠 Methods

Data Transformation: Split survival time into 10 intervals to handle censoring

Dimensionality Reduction: Factor Analysis applied to reduce predictors

Model Architecture: Deep neural network with 3 hidden layers (ReLU), output layer (Sigmoid)

Optimization: Adam optimizer to minimize cross-entropy loss

Benchmark Comparison:

Cox Proportional Hazards Model

DeepSurv

DeepHit

Random Survival Forest (RSF)

🔢 Performance Metrics

C-Index (Concordance Index)

Brier Score (BS)

Integrated Brier Score (IBS)

🗂 Repository Structure

fe-deepsurv/
├── data/                # Cirrhosis dataset (linked or included if license permits)
├── notebooks/           # Jupyter notebooks or R scripts for modeling and evaluation
├── src/                 # Source code for factor analysis and neural network
├── figures/             # Plots, survival curves, model comparison visuals
├── paper/               # Published paper or preprint version
├── README.md            # Project summary (this file)
└── requirements.txt     # Dependencies

📜 Citation

If you use this work, please cite:

@article{obite2025fe,
  title={Factor enhanced DeepSurv: A deep learning approach for predicting survival probabilities in cirrhosis data},
  author={Obite, Chukwudi Paul and others},
  journal={Computers in Biology and Medicine},
  volume={189},
  year={2025},
  doi={10.1016/j.compbiomed.2025.109963}
}
