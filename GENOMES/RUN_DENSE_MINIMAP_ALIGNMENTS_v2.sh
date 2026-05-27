#!/bin/bash

set -e

# ============================================================
# SETTINGS
# ============================================================

REF="H2_GNM.fasta"

OUTDIR="../DENSE_MINIMAP_ALIGNMENTS_v2"

THREADS=4

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/PAF
mkdir -p ${OUTDIR}/BED
mkdir -p ${OUTDIR}/BIGBED
mkdir -p ${OUTDIR}/SUMMARY

# ============================================================
# CHECKS
# ============================================================

echo
echo "====================================================="
echo "CHECKING SOFTWARE"
echo "====================================================="
echo

which minimap2 || { echo "minimap2 missing"; exit 1; }
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

echo
head ${OUTDIR}/H2.chrom.sizes

echo
echo "TOTAL SCAFFOLDS:"
wc -l ${OUTDIR}/H2.chrom.sizes

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

SUMMARY=${OUTDIR}/SUMMARY/dense_alignment_summary.tsv

echo -e "genome\tpaf_lines\tbed_lines" > ${SUMMARY}

# ============================================================
# RUN
# ============================================================

echo
echo "====================================================="
echo "RUNNING HIGH-SENSITIVITY DENSE ALIGNMENTS"
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
    # HIGH-SENSITIVITY MINIMAP2
    # --------------------------------------------------------

    minimap2 \
        -x asm10 \
        -c \
        --cs \
        -t ${THREADS} \
        ${REF} \
        ${GENOME} \
        > ${OUTDIR}/PAF/${PREFIX}_vs_H2.paf

    RAW=$(wc -l < ${OUTDIR}/PAF/${PREFIX}_vs_H2.paf)

    echo
    echo "RAW PAF LINES: ${RAW}"

    # --------------------------------------------------------
    # DENSE BED CONVERSION
    # --------------------------------------------------------

    awk '
    BEGIN{
        OFS="\t"
    }

    {

        qname=$1
        qstart=$3
        qend=$4
        strand=$5

        tname=$6
        tstart=$8
        tend=$9

        matches=$10
        alnlen=$11

        identity=(matches/alnlen)

        score=int(identity*1000)

        if(score < 0){score=0}
        if(score > 1000){score=1000}

        # ----------------------------------------------------
        # LOWER FILTERS FOR TRUE DENSE TRACKS
        # ----------------------------------------------------

        if(alnlen >= 100){

            print tname,
                  tstart,
                  tend,
                  qname,
                  score,
                  strand
        }
    }
    ' ${OUTDIR}/PAF/${PREFIX}_vs_H2.paf \
    > ${OUTDIR}/BED/${PREFIX}_vs_H2.bed

    BEDLINES=$(wc -l < ${OUTDIR}/BED/${PREFIX}_vs_H2.bed)

    echo "BED LINES: ${BEDLINES}"

    # --------------------------------------------------------
    # SORT BED
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

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    echo -e "${PREFIX}\t${RAW}\t${BEDLINES}" \
        >> ${SUMMARY}

done

# ============================================================
# FINISHED
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

