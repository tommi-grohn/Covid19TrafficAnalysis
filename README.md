# Analyzing the spread of Covid-19 with Google's traffic data in Finland, Norway and Sweden

Open right away the most important document, 'Report.pdf'. There is written all the essentials: abstract, introduction, methology, results and conclusion of this project.

In the file 'ModelFitting.ipynb' the stochastic SEAPIR-model is fitted for each region separately with a lot of comments. 
There is clearly described where the manually defined variable can be changed 'current_region'. By changing that value you can get 9 different fits overall (one for each region). 
This notebook uses 'SEAPIR.stan'.

The notebook 'DataCleaning.ipynb' uses Data.ods (Data collected from Finnish, Norwish and Swedish health institutions) and creates dataframes:
- 'DfRegions.csv' 
- 'DfInvestigatedTimePeriod.csv' 

To show that the model works at least for a simulated data, there has been created a notebook 'SimulatedModelFitting.ipynb'.

An interesting illustration how Google's traffic data also could be used is in 'AppendixTrafficAverages.pdf'.
