See README.md

Run in the folder 2_6_iqtree_boottrees_1000_congruent_round_2
```
for d in *.relabelled;
do for f in *.ph;
do (cd $d && cat $f | sed 's/;/;\n/g' > $d.boottrees);
done; done
```
We move the concatenated boottrees files to new directory 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate
```
for d in *.relabelled;
do for f in *.boottrees;
do (cd $d && mv $f ../../2_7_iqtree_boottrees_1000_congruent_round_2_concatenate);
done; done
```
We rename the boottrees files in 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate directory
```
for file in *.boottrees;
do mv "$file" "$(basename "$file" .boottrees.relabelled.boottrees).boottrees";
done
```


## 2.11) Re-label and compute the statistics of the boottrees (additional round) with congruent diploid skeleton for each directory/allele

We run the scripts/_bootstrap_label_stats.pl tool to re-label the boottrees and calculate the statistics of homeolog-types (A, B, C,...)
```
nohup scripts/_bootstrap_label_stats.pl \
2_7_iqtree_boottrees_1000_congruent_round_2_concatenate > stats_all_boottrees_congruent_round_2.tsv &
```


## 2.12) Random selection of 100 re-labelled boottrees (additional round) for each directory/allele

Run in 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate to copy to 2_8_iqtree_boottrees_100_random_congruent_round_2 directory
```
for d in *.relabelled;
do for f in *.root.ph;
do shuf -n 100 -e $d/$f | xargs cp --parents -t ../2_8_iqtree_boottrees_100_random_congruent_round_2;
done; done
```
We check the amount of boottrees in 2_8_iqtree_boottrees_100_random_congruent_round_2 directory

