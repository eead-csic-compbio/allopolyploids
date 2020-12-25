See README.md

## 2.7) Random selection of 100 boottrees (root.ph) for each directory/allele in the folder 2_3_iqtree_boottrees_1000_congruent

We run the commands in the folder 2_3_iqtree_boottrees_1000_congruent
```
for d in *.relabelled;
do for f in *.root.ph;
do shuf -n 100 -e $d/$f | xargs cp --parents -t ../2_4_iqtree_boottrees_100_random_congruent;
done; done

