#!/bin/bash

set -e

# ============================================================
# SETTINGS
# ============================================================

REF="H2_GNM.fasta"

OUTDIR="../LASTZ_BROWSER_ALIGNMENTS"

THREADS=4

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/LAV
mkdir -p ${OUTDIR}/AXT
mkdir -p ${OUTDIR}/BED
mkdir -p ${OUTDIR}/BIGBED
mkdir -p ${OUTDIR}/SUMMARY

# ============================================================
# CHECK SOFTWARE
# ============================================================

echo
echo "====================================================="
echo "CHECKING SOFTWARE"
echo "====================================================="
echo

which lastz || { echo "lastz missing"; exit 1; }
which bedToBigBed || { echo "bedToBigBed missing"; exit 1; }

# ============================================================
# BUILD CHROM.SIZES
# ============================================================

echo
echo "====================================================="
echo "BUILDING H2 CHROM.SIZES"
echo "====================================================="
echo

awk '
BEGIN{
    RS=">"
}
NR>1{

    split($0,a,"\n")

    header=a[1]

    seq=""

    for(i=2;i<=length(a);i++){
        seq=seq a[i]
    }

    gsub(/[ \t\r]/,"",seq)

    split(header,b," ")

    print b[1] "\t" length(seq)
}
' ${REF} > ${OUTDIR}/H2.chrom.sizes

# ============================================================
# GENOMES
# ============================================================

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

# ============================================================
# SUMMARY
# ============================================================

SUMMARY=${OUTDIR}/SUMMARY/lastz_summary.tsv

echo -e "genome\tbed_lines" > ${SUMMARY}

# ============================================================
# RUN
# ============================================================

echo
echo "====================================================="
echo "RUNNING LASTZ LOCAL ALIGNMENTS"
echo "====================================================="
echo

for GENOME in "${GENOMES[@]}"
do

    PREFIX=$(basename ${GENOME} .fasta)

    echo
    echo "-----------------------------------------------------"
    echo "ALIGNING ${PREFIX} vs H2"
    echo "-----------------------------------------------------"
    echo

    # --------------------------------------------------------
    # LASTZ
    # --------------------------------------------------------

    lastz \
        ${REF}[multiple] \
        ${GENOME}[multiple] \
        --step=20 \
        --seed=match12 \
        --notransition \
        --gfextend \
        --chain \
        --gapped \
        --format=general:name1,start1,end1,name2,start2,end2,strand2,identity,length1,length2 \
        > ${OUTDIR}/LAV/${PREFIX}_vs_H2.general

    # --------------------------------------------------------
    # GENERAL → BED
    # --------------------------------------------------------

    awk '
    BEGIN{
        OFS="\t"
    }

    NR>1{

        chr=$1
        start=$2
        end=$3

        qname=$4

        strand=$7

        ident=$8

        gsub(/%/,"",ident)

        score=int(ident*10)

        if(score > 1000){
            score=1000
        }

        alnlen=end-start

        if(alnlen >= 50){

            print chr,
                  start,
                  end,
                  qname,
                  score,
                  strand
        }
    }
    ' ${OUTDIR}/LAV/${PREFIX}_vs_H2.general \
    > ${OUTDIR}/BED/${PREFIX}_vs_H2.bed

    BEDLINES=$(wc -l < ${OUTDIR}/BED/${PREFIX}_vs_H2.bed)

    echo "BED LINES: ${BEDLINES}"

    # --------------------------------------------------------
    # SORT
    # --------------------------------------------------------

    sort -k1,1 -k2,2n \
        ${OUTDIR}/BED/${PREFIX}_vs_H2.bed \
        > ${OUTDIR}/BED/${PREFIX}_vs_H2.sorted.bed

    # --------------------------------------------------------
    # BIGBED
    # --------------------------------------------------------

    bedToBigBed \
        ${OUTDIR}/BED/${PREFIX}_vs_H2.sorted.bed \
        ${OUTDIR}/H2.chrom.sizes \
        ${OUTDIR}/BIGBED/${PREFIX}_vs_H2.bb

    echo
    echo "BUILT:"
    echo "${OUTDIR}/BIGBED/${PREFIX}_vs_H2.bb"

    echo -e "${PREFIX}\t${BEDLINES}" \
        >> ${SUMMARY}

done

# ============================================================
# DONE
# ============================================================

echo
echo "====================================================="
echo "DONE"
echo "====================================================="
echo

echo "BIGBED FILES:"
ls -lh ${OUTDIR}/BIGBED/*.bb

echo
echo "SUMMARY:"
echo ${SUMMARY}

echo

