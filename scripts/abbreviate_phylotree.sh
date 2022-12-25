sed -E 's|\[&Classification=\"[a-zA-Z \/]*\"\]||g'  ../data/phylo/SupFigs1-3-FigTree.tre \
| sed -E 's|\,Classification=\"[a-zA-Z \/]*\"\]|\]|g' \
| sed -E 's|\[\"Full_name\"=\".*\"\]||g' >../data/phylo/SupFigs1-3-FigTree_abr.tre

sed -E 's|\[&Classification=\"[a-zA-Z \/]*\"\]||g'  ../data\smirarab-1kp-b3326da\speciestrees\SupFigs1-3-FigTree.tre \
| sed -E 's|\,Classification=\"[a-zA-Z \/]*\"\]|\]|g' \
| sed -E 's|\[&ASTRAL_Caryo1.*Full_name.*\"\]||g' > ../data/phylo/SupFigs1-3-FigTree_abr.tre
