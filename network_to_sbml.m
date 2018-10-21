function err = network_to_sbml(stoichiometry, reversiblity, model_name)
  addpath('../CellNetAnalyzer/');
  startcna ;
  startcna(1);
  network.stoichMat = stoichiometry;
  network.cnapMin = reversiblity;
  network.cnapMin(reversiblity == 1) = -inf;
  cnap= CNAgenerateMFNetwork(network)
  mkdir('./models',model_name);
  err = CNAMFNetwork2sbml(cnap,strcat('./models/', model_name, '/', 'sbml.xml'));
