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
    sed -E 's/[ \t]+$//' $1 | sed -E 's/[ \t]+/\t/g'
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
# Searching & Highlighting
#
########################################

# pairgrep() {
#   pr --sep-string=' | ' -m -T -W 80 <(grep --color=always -e '$' -e $1 $2) <(grep --color=always -e '$' -e $1 $3) |less
# }

# Display corresponding reads side by side, with search sequence highlighted.
# e.g. pairgrep ACGT R1.fastq R2.fastq
pairgrep() {
    pr --sep-string=' | ' -m -T -W $COLUMNS $2 $3 |grep --color=always -e '$' -e $1 |less
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

function highlight2() {
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

function highlight() {
    declare -A color_map
    color_map[black]="$(tput setaf 0)"
    color_map[red]="$(tput setaf 1)"
    color_map[green]="$(tput setaf 2)"
    color_map[yellow]="$(tput setaf 3)"
    color_map[blue]="$(tput setaf 4)"
    color_map[magenta]="$(tput setaf 5)"
    color_map[cyan]="$(tput setaf 6)"
	
    fg_c=$color_map[$1]
    c_rs=`tput sgr0`

    sed -u -E "s/$2/$fg_c&$c_rs/g"
}

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

function expand() {
    echo $~1
    file_list=$~1
    echo $file_list
}
