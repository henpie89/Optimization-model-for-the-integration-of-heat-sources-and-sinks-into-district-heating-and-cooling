# -DTU
The copyright owner is the Technical University of Denmark.

The model is called "Optimization model for the integration of district heating, district cooling, heat sources and heat sinks"

Input data and description of the model can be found under the following DOI: 10.11583/DTU.11363198

The optimization model was developed in the General Algebraic Modelling System (GAMS, version 24.8.3) using the CPLEX solver, version 12.7.0.0. The software code can be found at the provided link. 
The purpose of the model is to investigate which heat sources, heat sinks and/or combination of them are best suited for heat pumps and chillers to supply district heating and cooling. Mixed-integer linear programming was used to carry out an annual optimization on an hourly basis to minimize total annualized costs.
The model takes investment costs, operating and maintenance (O&M) costs, short-term hot and cold water storage, seasonal temperature variations of heat sources, heat sinks and the district heating network, capacity limitations of the heat sources/sinks as well as the distance from the heat sources/sinks to the district heating and cooling network into account.
The model was developed with the intention of using it during an early planning stage of a new district heating and/or district cooling area when details about required production capacity are unknown and no decision about the heat source and heat sink has been made.

The PhD project was part of the project “EnergyLab Nordhavn - New Urban Energy Infrastructures” and was funded by EUDP (Energy Technology Development and Demonstration), project number: 64014-0555.
