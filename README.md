# democracyMarkovChain
Python code for the estimation of political regime transition matrices


This code accompanies the manuscript:
Quantifying the end of history through a Bayesian Markov-chain approach
by Florian Klimm

It allows the reproduction of all figures and results in the manuscript. In particular, this includes
- inference of Markov-chain transition probabilities from time-series data
- computation of expected regime changes under this Markov model
- estimation of the End of History as steady-state distribution
- estimation of Survival curves
- test the prediction of the Markov model

The associated Jupyter Notebooks are available under `/code`.

The numerical results and figures are included in `/results`.


### Prerequisites
- Python (tested for 3.7.10)
- some standard Python libraries (numpy, seaborn)
- [lifelines](https://github.com/CamDavidsonPilon/lifelines) for Kaplan-Meier estimates of survival functions
