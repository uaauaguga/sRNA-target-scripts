pt=models.b10.c64.with.targets.da0.5/36.pt
outdir=output/sRNA.KP.performance
mkdir -p $outdir
#scripts/inference.py  -f dataset/KP.test.sRNA.fa -o $outdir/scores.bed -m $pt
#scripts/performance-evaluation.py -i $outdir/scores.bed -o $outdir/performance.txt
scripts/performance-evaluation.py -i $outdir/scores.GCF_000742755.1.bed -o $outdir/scores.GCF_000742755.1.sRNA.txt
