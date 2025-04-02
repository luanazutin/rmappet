process rmats_run {
    container { params.containers.rmats }
    publishDir { "$params.outputdir/rmats/run/" }, mode: 'copy'
    label 'md'
    
    input:
    tuple val(conditions), path(condition_a_bams), path(condition_b_bams)
    path annotation

    output:
    tuple val(comparison), path("$comparison"), emit: data

    script:
    comparison = "${conditions.a}_vs_${conditions.b}"
    a_bams_joined = condition_a_bams.join(",")
    b_bams_joined = condition_b_bams.join(",")
    """
    # Preparação dos arquivos de entrada
    echo $a_bams_joined > ${conditions.a}.txt
    echo $b_bams_joined > ${conditions.b}.txt
    
    # Execução do rMATS
    mkdir -p $comparison
    rmats.py \\
        --b1 ${conditions.a}.txt \\
        --b2 ${conditions.b}.txt \\
        --novelSS \\
        -t $params.libtype \\
        --readLength $params.readlen \\
        --variable-read-length \\
        --gtf $annotation \\
        --nthread $task.cpus \\
        --od $comparison \\
        --tmp ${comparison}/tmp \\
        --task post \\
        --tstat

    # Validação crítica das colunas
    for file in $comparison/*.txt; do
        if ! head -n1 \$file | grep -q "transcript_ID"; then
            echo "ERRO CRÍTICO: transcript_ID não encontrado em \$file" >&2
            echo "Colunas disponíveis: \$(head -n1 \$file | tr '\t' '\n')" >&2
            exit 1
        fi
    done

    # Limpeza
    rm -r ${comparison}/tmp
    """
}

process rmats_parse_coords {
    container { params.containers.python }
    publishDir { "${params.outputdir}/rmats/results" }, mode: 'copy'
    label 'sm'

    input:
    tuple val(comparison), path(data)

    output:
    tuple val(comparison), path('*.jcec.tsv'), path('*.jc.tsv'), emit: data

    script:
    """
    # Processamento seguro com verificações explícitas
    for file in ${data}/*.txt; do
        awk -F'\t' '
        BEGIN {
            OFS = "\t";
            error = 0;
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if (\$i == "transcript_ID") trans_col = i;
                if (\$i == "chr") chr_col = i;
                if (\$i == "strand") strand_col = i;
                if (\$i == "exonStart_0base") start_col = i;
                if (\$i == "exonEnd") end_col = i;
            }
            if (trans_col == 0 || chr_col == 0 || strand_col == 0 || start_col == 0 || end_col == 0) {
                print "ERRO: Colunas obrigatórias não encontradas no cabeçalho" > "/dev/stderr";
                print "Cabeçalho atual: " \$0 > "/dev/stderr";
                error = 1;
                exit 1;
            }
            print "transcript_ID", "coordinates";
        }
        NR > 1 && !error {
            print \$trans_col, \$chr_col ":" \$start_col "-" \$end_col ":" \$strand_col;
        }' \$file > \$(basename \$file .txt).jcec.tsv || exit 1
    done
    """
}

process rmats_filter {
    container { params.containers.python }
    publishDir { "${params.outputdir}/rmats/results" }, mode: 'copy'
    label 'sm'

    input:
    tuple val(comparison), path(jcec), path(jc)

    output:
    tuple val(comparison), path('*.significant.jcec.tsv'), path('*.significant.jc.tsv')

    script:
    """
    # Filtro com verificações de colunas
    awk -F'\t' -v mindiff="$params.rmats.mindiff" -v maxfdr="$params.rmats.maxfdr" '
    BEGIN {
        OFS = "\t";
        error = 0;
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if (\$i == "transcript_ID") trans_col = i;
            if (\$i == "FDR") fdr_col = i;
            if (\$i == "IncLevelDifference") diff_col = i;
        }
        if (trans_col == 0 || fdr_col == 0 || diff_col == 0) {
            print "ERRO: Colunas de filtragem não encontradas" > "/dev/stderr";
            print "Colunas disponíveis: " \$0 > "/dev/stderr";
            error = 1;
            exit 1;
        }
        print \$0;
    }
    NR > 1 && !error {
        if (\$fdr_col <= maxfdr && \$diff_col >= mindiff) {
            print \$0;
        }
    }' ${jcec} > ${comparison}.significant.jcec.tsv || exit 1

    # Processamento opcional para arquivos JC
    if [ -e "${jc}" ]; then
        awk 'NR == 1 || (\$NF <= $params.rmats.maxfdr && \$(NF-1) >= $params.rmats.mindiff)' ${jc} \
            > ${comparison}.significant.jc.tsv
    else
        touch ${comparison}.significant.jc.tsv
    fi
    """
}
}
