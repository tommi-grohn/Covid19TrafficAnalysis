# Analysis Covid-19 with Google's traffic data in Finland, Norway and Sweden

The most important document is 'Report.pdf'. There is written all the essentials: abstract, introduction, methology, results and conclusion of this project.

In the file 'ModelFitting.ipynb' the stochastic SEAPIR-model is fitted for each region separately. 
The first manually defined variable is 'current_region'. By changing that value you can get 9 different fits overall. 
This notebook uses 'SEAPIR.stan'.

The notebook 'DataCleaning.ipynb' uses Data.ods (Data collected from Finnish, Norwish and Swedish health institutions) and creates dataframes 
- 'DfRegions.csv' 
- 'DfInvestigatedTimePeriod.csv' 

To show that the model works at least for a simulated data, there has been created a notebook 'SimulatedModelFitting.ipynb'

An interesting illustration how Google's traffic data also could be used is in 'AppendixTrafficAverages.pdf'
