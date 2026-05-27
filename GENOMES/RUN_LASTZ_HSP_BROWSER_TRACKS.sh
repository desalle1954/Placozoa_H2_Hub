#!/bin/bash
set -e

REF="H2_GNM.fasta"
OUTDIR="../LASTZ_HSP_BROWSER_TRACKS"
THREADS=4

mkdir -p ${OUTDIR}/{GENERAL,BED,BIGBED,SUMMARY}

which lastz || { echo "lastz missing"; exit 1; }
which bedToBigBed || { echo "bedToBigBed missing"; exit 1; }

awk '
BEGIN{RS=">"}
NR>1{
    split($0,a,"\n")
    split(a[1],b," ")
    seq=""
    for(i=2;i<=length(a);i++){seq=seq a[i]}
    gsub(/[ \t\r]/,"",seq)
    print b[1] "\t" length(seq)
}
' ${REF} > ${OUTDIR}/H2.chrom.sizes

GENOMES=(
H0_GNM.fasta
H1_GNM.fasta
H13_GNM.fasta
H15_GNM.fasta
H17_GNM.fasta
H19_GNM.fasta
H24_GNM.fasta
H25_GNM.fasta
)

SUMMARY=${OUTDIR}/SUMMARY/lastz_hsp_summary.tsv
echo -e "genome\tbed_lines\tbigbed_size" > ${SUMMARY}

for GENOME in "${GENOMES[@]}"
do
    PREFIX=$(basename ${GENOME} .fasta)

    echo
    echo "ALIGNING ${PREFIX} vs H2 with unchained LASTZ local HSP mode"

    lastz \
        ${REF}[multiple] \
        ${GENOME}[multiple] \
        --step=10 \
        --seed=match12 \
        --gapped \
        --format=general:name1,start1,end1,name2,start2,end2,strand2,identity,score \
        > ${OUTDIR}/GENERAL/${PREFIX}_vs_H2.general

    awk '
    BEGIN{OFS="\t"}
    NR>1{
        chr=$1
        start=$2-1
        end=$3
        q=$4
        strand=$7
        ident=$8
        gsub(/%/,"",ident)
        score=int(ident*10)
        if(score>1000){score=1000}
        if(score<0){score=0}
        if(start<0){start=0}
        len=end-start
        if(len>=20){
            print chr,start,end,q,score,strand
        }
    }
    ' ${OUTDIR}/GENERAL/${PREFIX}_vs_H2.general \
    > ${OUTDIR}/BED/${PREFIX}_vs_H2.bed

    sort -k1,1 -k2,2n \
        ${OUTDIR}/BED/${PREFIX}_vs_H2.bed \
        > ${OUTDIR}/BED/${PREFIX}_vs_H2.sorted.bed

    bedToBigBed \
        ${OUTDIR}/BED/${PREFIX}_vs_H2.sorted.bed \
        ${OUTDIR}/H2.chrom.sizes \
        ${OUTDIR}/BIGBED/${PREFIX}_vs_H2.hspLastz.bb

    BEDLINES=$(wc -l < ${OUTDIR}/BED/${PREFIX}_vs_H2.bed)
    BBSIZE=$(ls -lh ${OUTDIR}/BIGBED/${PREFIX}_vs_H2.hspLastz.bb | awk '{print $5}')

    echo -e "${PREFIX}\t${BEDLINES}\t${BBSIZE}" >> ${SUMMARY}
    echo "BED lines: ${BEDLINES}; bigBed size: ${BBSIZE}"
done

echo
echo "DONE:"
cat ${SUMMARY}
