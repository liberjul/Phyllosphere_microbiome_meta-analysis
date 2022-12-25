> orders_with_tips.txt

while IFS= read -r line
do
  ORDER=$(echo $line | grep -o "[[:alpha:]]*")
  echo $ORDER
  grep "$ORDER" ../data/smirarab-1kp-b3326da/speciestrees/SupFigs1-3-FigTree.tre | head -n1 | cut -f2 | cut -d[ -f1
  TIP=$(grep "$ORDER" ../data/smirarab-1kp-b3326da/speciestrees/SupFigs1-3-FigTree.tre | head -n1 | cut -f2 | cut -d[ -f1)
  echo $ORDER $TIP >> orders_with_tips.txt
done < ../orders_used.txt

cat orders_with_tips.txt
