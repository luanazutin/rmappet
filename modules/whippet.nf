process whippet_index {
  container { params.containers.whippet }
  label 'lg'

  input:
  path genome
  path annotation

  output:
  path 'index', emit: index

  script:
  """
  # Validação prévia dos arquivos de entrada
  if [ ! -s "$genome" ]; then
    echo "ERRO: Arquivo do genoma vazio ou não encontrado: $genome" >&2
    exit 1
  fi
  if [ ! -s "$annotation" ]; then
    echo "ERRO: Arquivo de anotação vazio ou não encontrado: $annotation" >&2
    exit 1
  fi

  # Configuração do ambiente Julia (se usando Singularity)
  ${ workflow.profile.contains('singularity') ? "julia --project=/code/whippet/ -e 'using Pkg; Pkg.instantiate()'" : "" }

  # Geração do índice
  mkdir -p index
  whippet-index.jl \\
    --fasta $genome \\
    --gtf $annotation \\
    -x index/graph \\
    --suppress-low-tsl

  # Validação pós-execução
  if [ ! -s "index/graph.jls" ]; then
    echo "ERRO: Falha na geração do índice Whippet" >&2
    exit 1
  fi
  """

  stub:
  """
  mkdir -p index
  echo "whippet-index.jl --fasta $genome --gtf $annotation -x index/graph --suppress-low-tsl" > index/command.txt
  touch index/graph.jls
  """
}

process whippet_quant {
  publishDir { "${params.outputdir}/whippet/quant" }, mode: 'copy'
  container { params.containers.whippet }
  label 'md'

  input:
  path index
  tuple val(sampleID), path(reads), val(condition)

  output:
  tuple val(sampleID), path('*.psi.gz'), val(condition), emit: psi
  path '*.gene.tpm.gz',                                  emit: gene_tpm
  path '*.isoform.tpm.gz',                               emit: isoform_tpm
  path '*.jnc.gz',                                       emit: junctions
  path '*.map.gz',                                       emit: map

  script:
  """
  # Validação do índice
  if [ ! -s "index/graph.jls" ]; then
    echo "ERRO: Índice Whippet não encontrado ou vazio" >&2
    exit 1
  fi

  # Quantificação
  whippet-quant.jl \\
    $reads \\
    -o $sampleID \\
    -x index/graph.jls

  # Validação da saída (garantindo que transcript_id está presente)
  for file in ${sampleID}*.psi.gz; do
    if ! zcat \$file | head -n1 | grep -q "transcript_id"; then
      echo "ERRO: transcript_id não encontrado em \$file" >&2
      zcat \$file | head -n1 | tr '\t' '\n' >&2
      exit 1
    fi
  done
  """

  stub:
  """
  echo "whippet-quant.jl $reads -o $sampleID -x index/graph.jls" > ${sampleID}.log
  touch ${sampleID}.psi.gz
  touch ${sampleID}.gene.tpm.gz
  touch ${sampleID}.isoform.tpm.gz
  touch ${sampleID}.jnc.gz
  touch ${sampleID}.map.gz
  """
}

process whippet_delta {
  publishDir { "${params.outputdir}/whippet/delta" }, mode: 'copy'
  container { params.containers.whippet }
  label 'md'

  input:
  tuple val(conditions), path(condition_a_quants), path(condition_b_quants)

  output:
  tuple val(comparison), path('*.diff.gz'), emit: delta

  script:
  comparison = "${conditions.a}_vs_${conditions.b}"
  """
  # Validação dos arquivos de entrada
  for file in ${condition_a_quants} ${condition_b_quants}; do
    if [ ! -s "\$file" ]; then
      echo "ERRO: Arquivo de quantificação vazio: \$file" >&2
      exit 1
    fi
  done

  # Análise diferencial
  whippet-delta.jl \\
    -a ${condition_a_quants.join(',')} \\
    -b ${condition_b_quants.join(',')} \\
    -o $comparison \\
    -s 2

  # Validação da saída
  if [ ! -s "${comparison}.diff.gz" ]; then
    echo "ERRO: Falha na análise diferencial" >&2
    exit 1
  fi
  """

  stub:
  comparison = "${conditions.a}_vs_${conditions.b}"
  """
  echo "whippet-delta.jl -a ${condition_a_quants.join(',')} -b ${condition_b_quants.join(',')} -o $comparison -s 2" > ${comparison}.log
  touch ${comparison}.diff.gz
  """
}

process whippet_filter {
  container { params.containers.python }
  publishDir { "${params.outputdir}/whippet/delta" }, mode: 'copy'
  label 'sm'

  input:
  tuple val(comparison), path(delta)

  output:
  tuple val(comparison), path('*.significant.tsv'), emit: data

  script:
  """
  # Validação prévia
  if [ ! -s "$delta" ]; then
    echo "ERRO: Arquivo .diff.gz vazio ou não encontrado" >&2
    exit 1
  fi

  # Filtro com garantia de transcript_id
  zcat $delta | head -n1 | grep -q "transcript_id" || {
    echo "ERRO: transcript_id não encontrado no arquivo de entrada" >&2
    exit 1
  }

  # Processamento
  whippet_filter.py \\
    $comparison \\
    $delta \\
    $params.whippet.mindiff \\
    $params.whippet.minprob \\
    --id-type transcript_id  # Garante uso explícito de transcript_id

  # Validação pós-processamento
  if [ ! -s "${comparison}.significant.tsv" ]; then
    echo "ERRO: Nenhum resultado significativo encontrado" >&2
    exit 1
  fi
  """

  stub:
  """
  echo "whippet_filter.py $comparison $delta $params.whippet.mindiff $params.whippet.minprob --id-type transcript_id" > ${comparison}.log
  touch ${comparison}.significant.tsv
  """
}
