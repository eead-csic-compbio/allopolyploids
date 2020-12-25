See README.md

## 2.9) Additional round of re-labelling congruent boottrees

There may be alleles where not all of your boottrees have been labelled. This can occur if one polyploid homeolog branch is inserted between the two outgroups (g. Osat and Hvul) and you have not been defined any ad-hoc labelling rules for this possibility because the genus under study is more recent than both outgroups. To get more consistent statistics in this cases, we re-label all boottrees with a congruent diploid skeleton (not just 100 boottrees) and conduct the random selection of 100 re-labelled boottrees.

Firstly, we create a list (e.g. list_alleles_round_2.txt in the folder 2_3_iqtree_boottrees_1000_congruent) of alleles (directories) that do not show 100 random boottrees re-labelled. We detect 79 alleles with some homeolog without labelling. We copy these to the filder 2_6_iqtree_boottrees_1000_congruent_round_2 and repeat the previous steps.
Run in the folder 2_3_iqtree_boottrees_1000_congruent
```
cat list_alleles_round_2.txt | xargs cp -r --parents -t ../2_6_iqtree_boottrees_1000_congruent_round_2/
```


## 2.10) We put together in a single boottrees file each of the 100 boottrees (additional round) of each allele.

Run in the folder 2_6_iqtree_boottrees_1000_congruent_round_2
```
for d in *.relabelled;
do for f in *.ph;
do (cd $d && cat $f | sed 's/;/;\n/g' > $d.boottrees);
done; done
```
We move the concatenated boottrees files to new directory 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate

