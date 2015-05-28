clc;
clear;
fclose all;

CLImAT('./data', './results', './genotype', './configuration/config.txt');
CLImAT_plot('./data', './results', './plots');
