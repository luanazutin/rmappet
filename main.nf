include { rmats_run; rmats_parse_coords; rmats_filter }                    from './modules/rmats.nf'
include { whippet_index; whippet_quant; whippet_delta; whippet_filter }    from './modules/whippet.nf'
include { overlap; overlap_filter }                                        from './modules/rmappet.nf'

// Funcao  para combinar amostras (mantida do original)
def combinations(channel) {
    channel
        .groupTuple(by: 0)
        .toList()
        .map {
            [it, it]
            .combinations()
            .findAll{ a, b -> a[0] < b[0] }
        }
        .flatMap()
        .map {
            it -> [
                [ a: it[0][0], b: it[1][0] ],
                it[0][1],
                it[1][1]
            ]
        }
}

// Canais de entrada (modificados para BAMs)
Channel.fromPath( params.samplesheet )
    .set { samplesheet_ch }
Channel.fromPath( params.genome )
    .set{ genome_ch }
Channel.fromPath( params.annotation )
    .set{ annotation_ch }

// Ler o samplesheet e espera uma coluna "bam" (nao mais "read1"/"read2")
samplesheet_ch.splitCsv( header: true )
    .take( params.dev ? 1: -1 )
    .map{ row -> [
            row.sample_id, file( row.bam ), row.condition  // Agora usa BAM!
    ]}
    .set { bam_ch }

// Workflow modificado (sem fastp/STAR/samtools)
workflow {
    // rMATS (usa BAMs diretamente)
    bam_ch
        .map { it -> [ it[2], it[1], it[0] ] }  // [condition, bam, sample_id]
        .set { rmats_ch }

    rmats_run( combinations( rmats_ch ), annotation_ch.collect() )
    rmats_parse_coords( rmats_run.out.data )
    rmats_filter( rmats_parse_coords.out.data )

    // Whippet (ainda precisa do indice, mas ignora alinhamento)
    whippet_index( genome_ch, annotation_ch )
    whippet_quant( whippet_index.out.index.collect(), bam_ch )  // Usa BAMs
    whippet_quant.out.psi
        .map { it -> [ it[2], it[1], it[0] ] }
        .set { whippet_ch }
    whippet_delta( combinations( whippet_ch ) )
    whippet_filter( whippet_delta.out.delta )

    // Overlap (mesmo do original)
    rmats_parse_coords.out.data
        .map { it -> [ it[0], it[1] ] }
        .join( whippet_delta.out.delta, by: 0 )
        .set { overlap_ch }
    
    overlap( overlap_ch )
    overlap_filter( overlap.out.data )
}
