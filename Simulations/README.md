<b>'simulate_de_novo_mutations.R' </b>

Script 'simulate_de_novo_mutations.R' generates synthetic dataset of de novo mutations based on the mutational spectrum. 

It takes the dataset with de novo mutations from public trios (wt dataset) and also the dataset of de novo mutations in the offspring of fathers with *POLD1* L474P variant from the current study (pold dataset). It calculates the 96-context mutational spectra for both datasets summarizing mutations by all available samples.

After this, it generates de novo mutations for the number of families equal to those in wt dataset. It samples mutations to this synthetic dataset according to their frequency in the wt de novo spectrum or pold de novo spectrum.The proportion of samples corresponding to the pold de novo spectrum can be changes using 'percent_of_pold_samples' variable (initially percent_of_pold_samples = 0).

Number of offspring generated per family can also be changed using 'n_simulations_per_family' variable (initially n_simulations_per_family=1)
