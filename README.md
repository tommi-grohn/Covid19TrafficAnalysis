# Covid19TrafficAnalysis

The most important document is 'report.pdf'. There is written abstract, introduction, methology, results and conclusion of this project.

In the file 'Model_Fitting.ipynb' the stochastic SEAPIR-model is fitted for each region separately. 
The first manually defined variable is 'current_region'. By changing that value you can get 9 different fits overall. 
This notebook uses 'SEAPIR_Model.stan'.

The notebook 'Data_Cleaning.ipynb' uses Data.ods (Data collected from Finnish, Norwish and Swedish health institutions) and creates dataframes 
- df_regions.csv 
- df_investigated_time_period.csv 

To show that the model works at least for a simulated data, there has been created a notebook SimulatedModelFitting.ipynb

An interesting illustration how traffic data also could be used is in AppendixTrafficAverages.pdf
