#!/bin/zsh

# sequence transformations
cdna() {tr "AaCcGgTt" "TtGgCcAa" <&0 >&1}
rdna() {rev <&0 >&1}
rcdna() {rev | tr "AaCcGgTt" "TtGgCcAa" <&0 >&1}

########################################
#
# File Manipulation
#
########################################

# remove spaces & tabs preceding the end of line
# convert each series of spaces & tabs to a single tab
# Not clear why this doesn't work: sed -E 's/[ \t]+$//' -E 's/[ \t]+/\t/g' $1
cleantsv() {
    sed -E 's/[[:space:]]+$//' $1 | sed -E 's/[[:blank:]]+/\t/g'
}

# generate shell variables for the sequences in a primer file
# file should have 3 columns: 1=label, 2=forward primer sequence, 3=reverse primer sequence
setseqvars() {$(awk '{print "export f"$1"="$2} {print "export r"$1"="$3}' $1)}

# Write one line for every pair of lines in file1, file2
allpairs() {
    while IFS= read -r line1; do
        while IFS= read -r line2; do
            echo $line1 '\t' $line2
        done < "$2"
    done < "$1"
}

# for paired primer file with fields: tag1, seq1, tag2, seq2,
# write a list of pairs as: tag1-tag2, seq1, seq2
# e.g., allpairs f_primers.tsv r_primers.tsv | allpairs
pairprimers() {
    allpairs $1 $2 | cleantsv | awk 'BEGIN{OFS="\t";} {print $1"-"$3,$2,$4}'
}

# extract sequences for forward & reverse primers from primers.txt or similarly formatted stream
# e.g., rprimers < primers.txt
rprimers() {awk '$4 ~ /R/ {print $2}' <&0 >&1}
fprimers() {awk '$4 ~ /F/ {print $2}' <&0 >&1}
rcrprimers() {rprimers <&0 | rcdna}
rcfprimers() {fprimers <&0 | rcdna}

# extract four sets of sequences from the primers.txt file, for forward, reverse, and reverse compliment of each
# first argument should be the name of the primers file `primers.txt`
# e.g. split_primers primers.txt
split_primers() {
    cat $1 | rprimers > reverse.txt
    cat $1 | fprimers > forward.txt
    cat $1 | rcrprimers > rcreverse.txt
    cat $1 | rcfprimers > rcforward.txt
}

########################################
#
# Highlighting & Viewing
#
########################################

export GREP_OPTIONS="--colour=always"

# pairgrep() {
#   pr --sep-string=' | ' -m -T -W 80 <(grep --color=always -e '$' -e $1 $2) <(grep --color=always -e '$' -e $1 $3) |less
# }


alias compare='pr -t -m -w 160'

# Display corresponding reads side by side, with search sequence highlighted.
# e.g. pairgrep ACGT R1.fastq R2.fastq
# pairgrep() {
#     pr --sep-string=' | ' -m -T -W $COLUMNS $2 $3 |grep --color=always -e '$' -e $1 |less
# }

pairgrep() {
    pr -m -t -w $COLUMNS $2 $3 |grep --color=always -e '$' -e $1 |less
}

prettytsv() {
    if [[ "$#" -eq 0 ]]; then
        csvtk pretty -t | less
    else
        csvtk pretty -t $1 | less
    fi
}

# highlight primers in different colors
# assumes the presence of four files containing the forward, reverse, reverse compliment forward & reverse primer sequences,
# as produced by `split_primers`
# e.g. phighlight seqeunces_R1.fastq
phighlight () {
    seqkit -w 0 fq2fa $1 | \
    GREP_COLOR='31;40' grep --colour=always -e '$' -f forward.txt | \
    GREP_COLOR='30;41' grep --colour=always -e '$' -f rcforward.txt | \
    GREP_COLOR='32;40' grep --colour=always -e '$' -f reverse.txt | \
    GREP_COLOR='30;42' grep --colour=always -e '$' -f rcreverse.txt | \
    less
}

phighlight_barcodes() {
    REF_PRE="TTCTTGACGAGTTCTTCTGA"
    REF_POST="ACGCGTCTGGAACAATCAAC"
    seqkit -w 0 fq2fa $1 | \
    GREP_COLOR='31;40' grep --colour=always -e '$' -f forward.txt | \
    GREP_COLOR='30;41' grep --colour=always -e '$' -f rcforward.txt | \
    GREP_COLOR='32;40' grep --colour=always -e '$' -f reverse.txt | \
    GREP_COLOR='30;42' grep --colour=always -e '$' -f rcreverse.txt | \
    GREP_COLOR='33;40' grep --colour=always -e '$' -e $REF_PRE | \
    GREP_COLOR='30;43' grep --colour=always -e '$' -e $REF_POST | \
    GREP_COLOR='30;106' grep --colour=always -E '([GC][AT]){9,}[GC]?|$' | \
    less
}

grep_sequence () {
    seqkit -w 0 fq2fa $2 | \
        GREP_COLOR='31;40;1' grep -e '$' -e $1 --colour=always | \
        GREP_COLOR='32;40' grep -e '$' -e $REF_PRE --colour=always | \
        GREP_COLOR='34;40' grep -e '$' -e $REF_POST --colour=always | \
        GREP_COLOR='36;40' grep -e '$' -e $SEQ_FWD --colour=always | \
        GREP_COLOR='35;40' grep -e '$' -e $SEQ_REV --colour=always | \
        less -S
}

function highlight2 () {
    # https://github.com/kepkin/dev-shell-essentials/blob/master/highlight.sh
    declare -A fg_color_map
    fg_color_map[black]=30
    fg_color_map[red]=31
    fg_color_map[green]=32
    fg_color_map[yellow]=33
    fg_color_map[blue]=34
    fg_color_map[magenta]=35
    fg_color_map[cyan]=36

    fg_c=$(echo -e "\e[1;${fg_color_map[$1]}m")
    c_rs=$'\e[0m'
    sed -u -E s"/($2)/$fg_c\1$c_rs/g"
}

# highlight regex matches
# e.g:
# > cat file.fasta | highlight red "ACTGACTG"
#
function highlight() {
    declare -A color_map
    color_map[black]="$(tput setaf 0)"
    color_map[red]="$(tput setaf 1)"
    color_map[green]="$(tput setaf 2)"
    color_map[yellow]="$(tput setaf 3)"
    color_map[blue]="$(tput setaf 4)"
    color_map[magenta]="$(tput setaf 5)"
    color_map[cyan]="$(tput setaf 6)"
    color_map[bold_red]="$(tput setaf 1)""$(tput setbf)"
    color_map[bold_green]="$(tput setaf 2)""$(tput setbf)"
    color_map[bold_blue]="$(tput setaf 4)""$(tput setbf)"
    color_map[bold_yellow]="$(tput setaf 3)""$(tput setbf)"

    fg_c=$color_map[$1]
    c_rs=`tput sgr0`

    sed -u -E "s/$2/$fg_c&$c_rs/g"
}

function dnaHighlight() {
    highlight bold_blue 'G'|highlight bold_red 'A'|highlight bold_green 'C'|highlight bold_yellow 'T'
}


# Given two TSV files with fields sample, frequency and barcode, display them side-by-side.
# Showing only rows with the matching sample name.
# Requires csvtk.

function compcounts() {
    pr -m -t -w $COLUMNS \
       <(csvtk grep -t -f sample -p $1 $2|csvtk cut -t -f frequency,barcode|csvtk pretty -t) \
       <(csvtk grep -t -f sample -p $1 $3|csvtk cut -t -f frequency,barcode|csvtk pretty -t)
}

# List the column names (from the first row)
colnames () {
    csvtk head $1 | csvtk transpose| csvtk cut -f 1
}

########################################
#
# Searching, etc
#
########################################

function count() {
    pattern_file='/Users/johann/bio/nijhawanlab/source/sequences.sh'
    #SEQ_REF_PRE="TTCTTGACGAGTTCTTCTGA"
    #SEQ_REF_POST="ACGCGTCTGGAACAATCAAC"
    #SEQ_FWD="GC.......GATATTGCTGAAGAGCTTG"
    #SEQ_REV="TT........CAGAGGTTGATTGTTCCAGA"
    #SEQ_FWD_RC="CAAGCTCTTCAGCAATATC.......GC"
    #SEQ_REV_RC="TCTGGAACAATCAACCTCTG........AA"
    #SEQ_TEST="GTCTGTCTCAGTCACACACAGTGT"

    declare -A pattern_lookup
    pattern_table=`awk -F= '{print $1,$2}' $pattern_file`
    for name value in `awk -F= '{print $1,$2}' $pattern_file`
    do
        #echo "$name -> $value"
        pattern_lookup[$name]=$value
    done

    # echo "pattern_lookup"
    #for key val in ${(@fkv)pattern_lookup}; do
    #    echo "$key -> $val"
    #done

    for key val in ${(@fkv)pattern_lookup}
    do
        cmd="grep --count -E $val $1"
        # echo $cmd
        echo "$key\t$val\t$1\t"$(grep --count -E ${val} ${1})
        #$($cmd)
        #echo "$key\t$val\t$1\t"`$cmd`
    done
}

# check to make sure sequence IDs are the same in two files
checkids () {
    diff <(seqkit seq -n -i $1) <(seqkit seq -n -i $2)
}

# A more concise output that counts the lines.
# The count is not accurate except that 0 means the IDs are the same (0 differences).
checkids2 () {
    diff -U 0 <(seqkit seq -n -i $1) <(seqkit seq -n -i $2) | grep -v ^@ | wc -l
}

# Pull sequence records from the forward, reverse & merged read FASTQ files,
# with IDs matching the argument passed to pull_seqs.  Write them to a file and
# run MUSCLE to align them.
# e.g.,
# pull_seqs VK003 $ID

function pull_seqs () {
    FWD_FILE=$1"_R1_001.fastq"
    REV_FILE=$1"_R2_001.fastq"
    MERGE_FILE=$1"_flash2.fastq"
    OUT_FILE="temp.ala"
    echo "Forward: $FWD_FILE Reverse: $REV_FILE Merged: $MERGE_FILE" >&2

    seqkit grep -p $2 $FWD_FILE | seqkit fq2fa -w0 > temp.fasta
    seqkit grep -p $2 $REV_FILE | seqkit seq -r -p -w 0 | seqkit fq2fa -w0 >> temp.fasta
    seqkit grep -p $2 $MERGE_FILE | seqkit fq2fa -w0 >> temp.fasta
    muscle -align temp.fasta -output $OUT_FILE
    echo "Alignment written to $OUT_FILE." >&2
    echo ">$2"
    seqkit seq -s $OUT_FILE
}

# testing shell glob expansion
function expand() {
    echo $~1
    file_list=$~1
    echo $file_list
}

########################################
#
# Amplicon-related
#
########################################

# Takes a single command-line argument specifying the prefix
# for the FASTQ files.

function flash_merge () {
    PREFIX=$1
    FWD_FILE=$1"_R1_001.fastq.gz"
    REV_FILE=$1"_R2_001.fastq.gz"
    LOG_FILE=$1"_flash.log"
    FRAGS_FILE=$1".extendedFrags.fastq"
    MERGED_FILE=$1"_merged.fastq"

    flash2 -M 230 -o $1 $FWD_FILE $REV_FILE 2>&1 | tee $LOG_FILE
    mv $FRAGS_FILE $MERGED_FILE
}

# Takes a single argument with the file name of a FASTQ
# file, e.g. the file containing the merged R1 & R2 reads.

guess_amplicon () {
    HCT116_SEQ="GTTCTTCTGASWSWSWSWSWSWSWSWSWSWSWSWACGCGTCTG"
    SUDHL5_SEQ="GGAGGGCTGASWSWSWSWSWSWSWSWSWSWTCTGGAACAA"
    N=$(seqkit grep -d -C -s -p "$HCT116_SEQ" $1)
    echo "HCT116: $N"
    N=$(seqkit grep -d -C -s -p "$SUDHL5_SEQ" $1)
    echo "SUDHL5: $N"
}

# Process a samples.tsv file for use in the Colab amplicons workflow.

procsamples () {
    cleantsv $1 | csvtk -t add-header --names "sample,fwd_primer,rev_primer"
}

########################################

echo "biocommands.sh sourced"
