# RichardsGrowthModel
R code for fitting a Richards growth model as multiple sub-models + model selection by predictive stacking

* DoubleRichardsModelSelection.r contains R code to evaluate a data example. It calls:
  + DoubleRichardsModelFunctions.r
  + grapedata.csv
* SimulationStudy_Richards_model_stacking.r contains R code to simulate the performance of the predictive stacking procedure for the Richards model. The import of the code in
  + RichardsModelFunctions.r is required
* stacking_simulation_results.Rda contains simulation results, which can be loaded within the simulation code.
