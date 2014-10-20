n=1000 #genes
m=1000 #experiments

#Find all possible gene combinations from 2 to 1000 out of 1000 genes
total_gene_combos=0
for (k in 2:n){
gene_combos=choose(n,k)
total_gene_combos=total_gene_combos+gene_combos
}
total_gene_combos

#Find all possible exp combinations from 10 to 1000 out of 1000 experiments
total_exp_combos=0
for (k in 10:m){
exp_combos=choose(m,k)
total_exp_combos=total_exp_combos+exp_combos
}
total_exp_combos

#The total number of subspaces is 
#the number of gene combinations times the number of experiment combinations
total_subspaces=total_gene_combos*total_exp_combos
total_subspaces
