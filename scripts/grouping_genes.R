library(tidyverse)
library(ggplot2)

setwd("~/Desktop/jacob/mcycdb/MCycDB/")


metadata <- read.csv("metadata_filtered.csv")


metadata$cycle <- NA




cenmetpat <- c('mtrA', 'mtrB', 'mtrC', 'mtrD', 'mtrE', 'mtrF', 'mtrG', 'mtrH', 
               'mcrA', 'mcrB', 'mcrC', 'mcrD', 'mcrG', 'hdrA1', 'hdrB1', 'hdrC1',
               'hdrA2', 'hdrB2', 'hdrC2', 'hdrD', 'hdrE', 'mvhA', 'mvhD', 'mvhG',
               'ehaA', 'ehaB', 'ehaC', 'ehaD', 'ehaE', 'ehaF', 'ehaG', 'ehaH', 
               'ehaI', 'ehaJ', 'ehaK', 'ehaL', 'ehaM', 'ehaN', 'ehaO', 'ehaP', 
               'ehaQ', 'ehaR', 'ehbA', 'ehbB', 'ehbC', 'ehbD', 'ehbE', 'ehbF', 
               'ehbG', 'ehbH', 'ehbI', 'ehbJ', 'ehbK', 'ehbL', 'ehbM', 'ehbN', 
               'ehbO', 'ehbP', 'ehbQ', 'mbhL', 'mbhK', 'mbhJ', 'rnfA', 'rnfB', 
               'rnfC', 'rnfD', 'rnfE', 'rnfG', 'echA', 'echB', 'echC', 'echD', 
               'echE', 'echF', 'vhoA', 'vhoC', 'vhoG', 'vhtA', 'vhtC', 'vhtG', 
               'vhuA', 'vhuU', 'vhuD', 'vhuG', 'vhcA', 'vhcD', 'vhcG', 'fpoA', 
               'fpoB', 'fpoC', 'fpoD', 'fpoF', 'fpoH', 'fpoI', 'fpoJ', 'fpoK', 
               'fpoL', 'fpoM', 'fpoN', 'fpoO', 'fqoA', 'fqoD', 'fqoF', 'fqoH', 
               'fqoJ', 'fqoK', 'fqoL', 'fqoM', 'fqoN', 'frhA', 'frhB', 'frhD', 
               'frhG', 'fruA', 'fruB', 'fruD', 'fruG', 'frcA', 'frcB', 'frcD', 
               'frcG')


hydro_met <- c('fmdA', 'fmdB', 'fmdC', 'fmdD', 'fmdE', 'fmdF', 'fwdA', 'fwdB', 
               'fwdC', 'fwdD', 'fwdE', 'fwdF', 'fwdG', 'fwdH', 'ftr', 'mch', 
               'mtdA', 'mtdB', 'hmd', 'mer', 'metF')

aceti_met <- c('acs', 'acsA', 'acsC', 'acsD', 'ackA', 'pta', 'cdhC', 'cdhD', 
               'cdhE', 'cdhA', 'cdhB', 'acdA', 'acdB', 'cooS', 'cooF')

methyl_met <- c('mtaA', 'mtaB', 'mtaC', 'mtbA', 'mtbB', 'mtbC', 'mtsA', 'mtsB', 
                'torA', 'torC', 'torD', 'torY', 'torZ', 'mttB', 'mttC', 'mtmB', 'mtmC')

aom <- c('fmdA', 'fmdB', 'fmdC', 'fmdD', 'fmdE', 'fmdF', 'fwdA', 'fwdB', 'fwdC', 
         'fwdD', 'fwdE', 'fwdF', 'fwdG', 'fwdH', 'ftr', 'mch', 'mtdA', 'mtdB', 
         'hmd', 'mer', 'metF', 'mtrA', 'mtrB', 'mtrC', 'mtrD', 'mtrE', 'mtrF',
         'mtrG', 'mtrH', 'mcrA', 'mcrB', 'mcrC', 'mcrD', 'mcrG', 'hdrA1', 'hdrB1', 
         'hdrC1', 'hdrA2', 'hdrB2', 'hdrC2', 'hdrD', 'hdrE', 'rnfA', 'rnfB', 'rnfC', 
         'rnfD', 'rnfE', 'rnfG', 'frhA', 'frhB', 'frhD', 'frhG', 'fruA', 'fruB', 
         'fruD', 'fruG', 'frcA', 'frcB', 'frcD', 'frcG', 'fpoA', 'fpoB', 'fpoC', 
         'fpoD', 'fpoF', 'fpoH', 'fpoI', 'fpoJ', 'fpoK', 'fpoL', 'fpoM', 'fpoN', 
         'fpoO', 'fqoA', 'fqoD', 'fqoF', 'fqoH', 'fqoJ', 'fqoK', 'fqoL', 'fqoM', 
         'fqoN', 'narG', 'narZ', 'narH', 'narY', 'narB', 'napH', 'nrfH', 'nrfA', 
         'nxrA', 'nod', 'nirK', 'nirS', 'cytC')

oxid_met_c1 <- c('pqqA', 'pqqB', 'pqqC', 'pqqD', 'pqqE', 'pqqF', 'mmoX', 'mmoY', 
                 'mmoZ', 'mmoB', 'mmoC', 'mmoD', 'pmoA', 'pmoB', 'pmoC', 'amoA', 
                 'amoB', 'amoC', 'xoxF1', 'xoxF2', 'xoxF4', 'xoxF5', 'mxaF', 'mxaI',
                 'mxaJ', 'mxaG', 'mxaA', 'mxaC', 'mxaK', 'mxaL', 'mxaD', 'mauA', 'mauB',
                 'mauC', 'mauD', 'mauE', 'mauF', 'mgsA', 'mgsB', 'mgsC', 'mgdA', 
                 'mgdB', 'mgdD', 'qhpA', 'dcmR', 'dcmA', 'mdh', 'tmm')

oxid_formaldehyde <- c('gfa', 'frmA', 'frmB', 'fghA', 'fae', 'mtdB', 'mch', 'mdo',
                       'fdm', 'fdhA-K00148')

oxid_formate <- c('fdwA', 'fdwB', 'fdwE', 'fdsB', 'fdsD', 'fdsG', 'fdoG', 'fdoH', 
                  'fdoI', 'fdhF', 'fdhA', 'fdhB')

serine <- c('fchA', 'mtdA', 'dfrA1', 'dfrA12', 'dfrA10', 'dfrA19', 'folA', 'glyA',
            'hprA', 'gck', 'eno', 'ppc', 'mtkA', 'mtkB', 'mcl', 'gpmI', 'gpmB',
            'apgM', 'serA', 'serC', 'serB', 'thrH', 'psp', 'porA', 'porB', 'porD', 
            'porG', 'pps', 'mdh-K00024')

rump <- c('hxlB', 'hxlA', 'pfkA', 'pfkB', 'pfkC', 'pfp', 'fbp', 'fbaA', 'fbaB', 'glpX',
          'fbp3', 'fae-hps', 'fbp-SEBP')


for(value in cenmetpat){
  if(metadata$cycle[metadata$gene == value] == ){
    metadata$cycle[metadata$gene == value] <- "cenmetpat"
  }
  else{
    metadata$cycle[metadata$gene == value] <- paste0(metadata$cycle[metadata$gene == value], ",cenmetpat")
    
  }
}

for(value in hydro_met){
  metadata$cycle[metadata$gene == value] <- "hydro_met"
}

for(value in aceti_met){
  metadata$cycle[metadata$gene == value] <- "aceti_met"
}

for(value in methyl_met){
  metadata$cycle[metadata$gene == value] <- "methyl_met"
}

for(value in aom){
  metadata$cycle[metadata$gene == value] <- "aom"
}

for(value in oxid_met_c1){
  metadata$cycle[metadata$gene == value] <- "oxid_met_c1"
}

for(value in oxid_formaldehyde){
  metadata$cycle[metadata$gene == value] <- "oxid_formaldehyde"
}

for(value in oxid_formate){
  metadata$cycle[metadata$gene == value] <- "oxid_formate"
}

for(value in serine){
  metadata$cycle[metadata$gene == value] <- "serine"
}

for(value in rump){
  metadata$cycle[metadata$gene == value] <- "rump"
}


write.csv(metadata, "metadata_wgroups.csv", row.names = F)


check <- filter(metadata, is.na(cycle) == TRUE)
