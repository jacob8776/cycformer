
gene_list = ["mtrABCDEFGH", "mcrABCDG", "hdrA1B1C1A2B2C2DE", "mcrABCDG", "hdrA1B1C1A2B2C2DE", "mvhADG", 
"ehaABCDEFGHIJKLMNOPQR", "ehbABCDEFGHIJKLMNOPQ", "mbhLKJ", "rnfABCDEG", 
"echABCDEF", "vhoACG", "vhtACG", "vhuAUDG", "vhcADG", "fpoABCDFHIJKLMNO", 
"fqoADFHJKLMN", "frhABDG", "fruABDG", "frcABDG"]

cenmetpat = []

for gene in gene_list: 
  for i in gene[3:]:
    if i.isalpha(): 
      cenmetpat.append(gene[0:3] + i)
    else: 
      cenmetpat[-1] = cenmetpat[-1] + i
    




