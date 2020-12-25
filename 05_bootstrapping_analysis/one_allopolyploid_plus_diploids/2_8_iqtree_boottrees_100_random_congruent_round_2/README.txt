See README.md

## 2.12) Random selection of 100 re-labelled boottrees (additional round) for each directory/allele

Run in 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate to copy to 2_8_iqtree_boottrees_100_random_congruent_round_2 directory
```
for d in *.relabelled;
do for f in *.root.ph;
do shuf -n 100 -e $d/$f | xargs cp --parents -t ../2_8_iqtree_boottrees_100_random_congruent_round_2;
done; done
```
We check the amount of boottrees in 2_8_iqtree_boottrees_100_random_congruent_round_2 directory

