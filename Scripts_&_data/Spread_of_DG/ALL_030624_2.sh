# mafft ALL_030624_1.fas > ALL_030624_1_NEW.fas
# mafft ALL_030624_2.fas > ALL_030624_2_NEW.fas

IQTREE_v_1_6_12/iqtree -s ALL_030624_2.fas -m MFP -mem 10Go -mset GTR -nt 4 # -b 100

BEAGLEPATH="/usr/local/lib/"
java -Djava.library.path="$BEAGLEPATH" -jar beast_1_10_5.jar ALL_030624_2.xml
