# installing "nodadi" version of scripts
cd
rm -rf AFS-analysis-with-moments
git clone https://github.com/z0on/AFS-analysis-with-moments.git
cp ~/AFS-analysis-with-moments/multimodel_inference/nodadi/* ~/bin

# jp's saf data
cp /work/05703/jprippe/jp_shared/mcav/* .

# angsd bootstrap (with sleep in front to make sure bootstraps use different random seeds)
export GENOME_REF=Mcavernosa_July2018_phased.fasta
>b100
for B in `seq 1 100`; do
echo "sleep $B && realSFS c2.saf.idx c4.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 6 -P 1 -resample_chr 1 >c24_$B
sleep $B && realSFS c1.saf.idx c4.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 6 -P 1 -resample_chr 1 >c14_$B
sleep $B && realSFS c1.saf.idx c3.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 6 -P 1 -resample_chr 1 >c13_$B
sleep $B && realSFS c1.saf.idx c2.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 6 -P 1 -resample_chr 1 >c12_$B
sleep $B && realSFS c2.saf.idx c3.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 6 -P 1 -resample_chr 1 >c23_$B
sleep $B && realSFS c3.saf.idx c4.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 6 -P 1 -resample_chr 1 >c34_$B" >>b100;
done
ls5_launcher_creator.py -j b100 -n b100 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch b100.slurm


# sizes of 1d frequency spectra
>sizes
for C in `seq 1 4`; do
realSFS c${C}.saf.idx >c${C}.sfs1;
cat c${C}.sfs1 | wc -w >>sizes;
done

cat sizes
75
63
47
17

# "bagging" (averaging) by 5 bootstrap reps (removing the first one since it is for the whole dataset, according to Nate Pope)
for B in `seq 1 100`; do
echo "63 17" >c24_${B}.sfs;
tail -5 c24_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> c24_${B}.sfs;
echo "75 17" >c14_${B}.sfs;
tail -5 c14_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> c14_${B}.sfs;
echo "47 17" >c34_${B}.sfs;
tail -5 c34_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> c34_${B}.sfs;
echo "75 63" >c12_${B}.sfs;
tail -5 c12_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> c12_${B}.sfs;
echo "63 47" >c23_${B}.sfs;
tail -5 c23_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> c23_${B}.sfs;
echo "75 47" >c13_${B}.sfs;
tail -5 c13_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> c13_${B}.sfs;
done

cp *.sfs /work/05703/jprippe/jp_shared/

for K in `seq 1 2`; do
echo $K;
2dSFSplot_angsd_p.py c12_${K}.sfs 74 62;
2dSFSplot_angsd_p.py c12_${K}.sfs 68 56;
2dSFSplot_angsd_p.py c12_${K}.sfs 60 50;
done

# will use 80% down-projection

# ----------- running multimodel analysis: all models galore, 10 bootstrap replicates (in 2 batches of 3 starts)

cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
NREPS=3
>mods
for i in `seq 1 $NREPS`;do 
	cat allmodels >>mods;
done

#for F in *sfs1; do
#F=`echo $F | perl -pe 's/c(.+)\..+/$1/'`;
#echo $F;
#done

CONTRAST=c12
ARGS="c1 c2 60 50 0.02 0.005"

>100runs1
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runs1;
done
ls5_launcher_creator.py -j 100runs1 -n 100runs1 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runs1.slurm


>100runs2
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runs2;
done
ls5_launcher_creator.py -j 100runs2 -n 100runs2 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runs2.slurm


>100runs5
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runs5;
done
ls5_launcher_creator.py -j 100runs5 -n 100runs5 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runs5.slurm

>100runs6
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runs6;
done
ls5_launcher_creator.py -j 100runs6 -n 100runs6 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runs6.slurm


grep RESULT 100runs[1256].o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c12.res

# best likelihood:
cat c12.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'
grep "sc3m " c12.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5,6 -d " " c12.res >c12.likes

grep "sc3ielsm " c12.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

#-------------------- c34 

cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
NREPS=3
>mods
for i in `seq 1 $NREPS`;do 
	cat allmodels >>mods;
done

CONTRAST=c34
ARGS="c3 c4 38 12 0.02 0.005"

>100runs3
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runs3;
done
ls5_launcher_creator.py -j 100runs3 -n 100runs3 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runs3.slurm


>100runs4
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runs4;
done
ls5_launcher_creator.py -j 100runs4 -n 100runs4 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runs4.slurm

grep RESULT 100runs[34].o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c34.res

# best likelihood:
cat c34.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5,6 -d " " c34.res >c34.likes

grep "sc2ielsm1 " c34.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'


#------- c13


cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
NREPS=3
>mods
for i in `seq 1 $NREPS`;do 
cat allmodels >>mods;
done

CONTRAST=c13
ARGS="c1 c3 60 38 0.02 0.005"

>100runsc13.1
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc13.1;
done
ls5_launcher_creator.py -j 100runsc13.1 -n 100runsc13.1 -t 1:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc13.1.slurm


>100runsc13.2
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc13.2;
done
ls5_launcher_creator.py -j 100runsc13.2 -n 100runsc13.2 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc13.2.slurm


grep RESULT 100runsc13.*.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c13.res

# best likelihood:
grep "sc2iel " c13.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5,6 -d " " c13.res >c13.likes

#------- c14

cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
NREPS=3
>mods
for i in `seq 1 $NREPS`;do 
cat allmodels >>mods;
done

CONTRAST=c14
ARGS="c1 c4 60 12 0.02 0.005"

>100runsc14.1
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc14.1;
done
ls5_launcher_creator.py -j 100runsc14.1 -n 100runsc14.1 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc14.1.slurm


>100runsc14.2
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc14.2;
done
ls5_launcher_creator.py -j 100runsc14.2 -n 100runsc14.2 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc14.2.slurm

grep RESULT 100runsc14.*.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c14.res

# best likelihood:
grep "sc2ielsm2 " c14.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5,6 -d " " c14.res >c14.likes


#------- c24

cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
NREPS=3
>mods
for i in `seq 1 $NREPS`;do 
cat allmodels >>mods;
done

CONTRAST=c24
ARGS="c2 c4 50 12 0.02 0.005"

>100runsc24.1
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc24.1;
done
ls5_launcher_creator.py -j 100runsc24.1 -n 100runsc24.1 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc24.1.slurm


>100runsc24.2
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc24.2;
done
ls5_launcher_creator.py -j 100runsc24.2 -n 100runsc24.2 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc24.2.slurm

grep RESULT 100runsc24.*.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c24.res

# best likelihood:
grep "sc2ielsm1 " c24.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5,6 -d " " c24.res >c24.likes


#------- c23

cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
NREPS=3
>mods
for i in `seq 1 $NREPS`;do 
cat allmodels >>mods;
done

CONTRAST=c23
ARGS="c2 c3 50 38 0.02 0.005"

>100runsc23.1
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc23.1;
done
ls5_launcher_creator.py -j 100runsc23.1 -n 100runsc23.1 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc23.1.slurm


>100runsc23.2
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsc23.2;
done
ls5_launcher_creator.py -j 100runsc23.2 -n 100runsc23.2 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsc23.2.slurm

grep RESULT 100runsc23.*.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c23.res

# best likelihood:
grep "sc2ielsm1 " c23.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5,6 -d " " c23.res >c23.likes


#----- 100 bootstraps for the winning model

module load python2
CONTRAST=c12
ARGS="c1 c2 60 50 0.02 0.005"
WINNER="SC3ielsm1.py"

NREPS=6
>mods
for i in `seq 1 $NREPS`;do 
echo $WINNER >>mods;
done

>100runsw1
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsw1;
done
ls5_launcher_creator.py -j 100runsw1 -n 100runsw1 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsw1.slurm

grep RESULT 100runsw1.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c12wins.res
# best fit
cat c12wins.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'




CONTRAST=c34
ARGS="c3 c4 38 12 0.02 0.005"
WINNER="SC2ielsm1.py"

NREPS=6
>mods
for i in `seq 1 $NREPS`;do 
echo $WINNER >>mods;
done

>100runsw2
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100runsw2;
done
ls5_launcher_creator.py -j 100runsw2 -n 100runsw2 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100runsw2.slurm

grep RESULT 100runsw2.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c34wins.res
# best fit
cat c34wins.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'


#----- 100 boots c13

module load python2
CONTRAST=c13
WINNER=`ls ${CONTRAST}.res.* | perl -pe 's/.+\.//'`
ARGS="c1 c3 60 38 0.02 0.005"

NREPS=6
>mods
for i in `seq 1 $NREPS`;do 
echo $WINNER >>mods;
done

>100b.c13
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100b.c13;
done

ls5_launcher_creator.py -j winboots -n winboots -t 1:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch winboots.slurm

ls5_launcher_creator.py -j 100b.c13 -n winboots -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100b.c13.slurm


grep RESULT 100b.c13.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c13wins.res
# best fit
cat c13wins.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

#----- 100 boots c23

module load python2
CONTRAST=c23
ARGS="c2 c3 50 38 0.02 0.005"
WINNER="SC2ielsm1.py"

NREPS=6
>mods
for i in `seq 1 $NREPS`;do 
echo $WINNER >>mods;
done

>100b.c23
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100b.c23;
done
ls5_launcher_creator.py -j 100b.c23 -n 100b.c23 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100b.c23.slurm

grep RESULT 100b.c23.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c23wins.res
# best fit
cat c23wins.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

#----- 100 boots c14

module load python2
CONTRAST=c14
ARGS="c1 c4 60 12 0.02 0.005"
WINNER="SC2ielsm2.py"

NREPS=6
>mods
for i in `seq 1 $NREPS`;do 
echo $WINNER >>mods;
done

>100b.c14
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100b.c14;
done
ls5_launcher_creator.py -j 100b.c14 -n 100b.c14 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100b.c14.slurm

grep RESULT 100b.c14.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c14wins.res
# best fit
cat c14wins.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'

#----- 100 boots c24

module load python2
CONTRAST=c24
ARGS="c2 c4 50 12 0.02 0.005"
WINNER="SC2ielsm1.py"

NREPS=6
>mods
for i in `seq 1 $NREPS`;do 
echo $WINNER >>mods;
done

>100b.c24
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>100b.c24;
done
ls5_launcher_creator.py -j 100b.c24 -n 100b.c24 -t 2:00:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch 100b.c24.slurm

grep RESULT 100b.c24.o* -A 4 | grep -v Launcher | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >c24wins.res
# best fit
cat c24wins.res | awk 'NR == 1 || $5 > max {line = $0; max = $5}END{print line}'



#=========================================================

#--------------- GADMA

#installing dadi
cd
git clone https://bitbucket.org/gutenkunstlab/dadi.git
cd dadi
python setup.py install --user

# installing Pillow
python -m pip install Pillow

# installing GADMA
cd
rm -rf GADMA
git clone https://github.com/ctlab/GADMA.git
cd GADMA
python setup.py install --user


# writing GADMA parameters file

echo "# param_file
Input file : c12_100.sfs_60_50
Output directory : gadma_12dp
Population labels : c1 , c2
Initial structure : 1,1" >gadma_params

rm -rf gadma_12*
echo "gadma -p params">gad
echo "gadma --resume gadma_12dp">gad
ls5_launcher_creator.py -j gad -n gad -t 48:00:00 -e matz@utexas.edu -w 1 -a tagmap -q normal
sbatch gad.slurm


