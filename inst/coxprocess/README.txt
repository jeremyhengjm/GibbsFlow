Log-Gaussian Cox process model
- Run demo_gibbsvelocity_dim.R, demo_gibbsvelocity_vi.R and demo_gibbsvelocity_ep.R followed by analyse_gibbsflow_integration_dim.R and analyse_gibbsflow_integration_initial.R to investigate stiffness of Gibbs flow
- Run demo_variational_inference.R to compute variational Bayes approximation of the posterior distribution
- Run demo_expectation_propagation.R to compute expectation-propagation approximation of the posterior distribution
- Run demo_ais, demo_gibbsflow_sis, demo_gibbsflow_ais files and demo_svgd.R to see how each method performs on this model
- Run repeat files using repeats.R followed by analyse_initial_results.R and analyse_scaling_result.R to compare performance
- To run on the coarser 10 x 10 grid, line 6 of src/gibbsflow_coxprocess.cpp should also be changed to const int ngrid = 10;
- See transportmaps.py for an implementation of the transport methodology developed by Youssef Marzouk and his research group. This builds on the TransportMaps Python package http://transportmaps.mit.edu/docs/ 
