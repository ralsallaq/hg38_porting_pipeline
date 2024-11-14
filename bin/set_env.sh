module load ucsc
module load vcftools
module load blat
cbload () {
    local path=/research/rgs01/resgen/system/sjcbinit/cbload.sh;
    if [ ! -f $path ]; then
        echo "cbload.sh not found at $path" 1>&2;
    else
        source $path "$@";
    fi
}

setcbenv ()
{
    local path=/research/rgs01/resgen/system/sjcbinit/setcbenv.sh;
    if [ ! -f $path ]; then
        echo "setcbenv.sh not found at $path" 1>&2;
    else
        source $path "$@";
    fi
}

setcbenv dev branches/hg38
cbload --set official
cbload configs
cbload common-scripts-internal
cbload seq_anls_ops
cbload snv-annovar
