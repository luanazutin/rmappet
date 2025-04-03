include { rmats_run; rmats_parse_coords; rmats_filter }                    from './modules/rmats.nf'
include { whippet_index; whippet_quant; whippet_delta; whippet_filter }    from './modules/whippet.nf'
include { overlap; overlap_filter }                                        from './modules/rmappet.nf'

// Functions
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

// Channels
Channel.fromPath( params.samplesheet )
    .splitCsv( header: true )
    .take( params.dev ? 1: -1 )
    .map{ row -> [ row.sample_id, file( row.bam ), row.condition ] }
    .set { bam_ch }

Channel.fromPath( params.annotation )
    .set { annotation_ch }

// Workflow 
workflow {
    
    // rMATS (análise de splicing alternativo)
    rmats_run( combinations(samtools_sort.out.bam.map{ it -> [ it[2], it[1][0], it[0] ] }), annotation_ch.collect() )
    rmats_parse_coords( rmats_run.out.data )
    rmats_filter( rmats_parse_coords.out.data )
    
    // Whippet (análise de eventos de splicing)
    whippet_index( annotation_ch )
    whippet_quant( whippet_index.out.index.collect(), samtools_sort.out.bam.map{ it -> [ it[0], it[1][0] ] } )
    whippet_delta( combinations(whippet_quant.out.psi.map{ it -> [ it[2], it[1], it[0] ] }) 
    whippet_filter( whippet_delta.out.delta )
    
    // Overlap (integração rMATS + Whippet)
    rmats_parse_coords.out.data
        .map { it -> [ it[0], it[1] ] }
        .join( whippet_delta.out.delta, by: 0 )
        .set { overlap_ch }
    
    overlap( overlap_ch )
    overlap_filter( overlap.out.data )
}
