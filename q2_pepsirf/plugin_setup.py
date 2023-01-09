#!/usr/bin/env python
import importlib
import q2_pepsirf

from qiime2.plugin import (Plugin,
                        SemanticType,
                        model,
                        Int,
                        Range,
                        MetadataColumn,
                        Categorical,
                        Str,
                        List,
                        Visualization,
                        Metadata,
                        TypeMap,
                        Choices,
                        Float,
                        Bool)

from q2_pepsirf.format_types import (
    EnrichedPeptideDirFmt, Normed, NormedDifference, 
    NormedDiffRatio, NormedRatio, NormedSized, 
    PairwiseEnrichment, PeptideIDListFmt, ProteinAlignment, Zscore, RawCounts,
    PepsirfContingencyTSVDirFmt, PepsirfContingencyTSVFormat, 
    ZscoreNan, ZscoreNanDirFmt, ZscoreNanFormat, PeptideBinFormat,
    PeptideBinDirFmt, PeptideBins, PepsirfInfoSumOfProbesDirFmt,
    PepsirfInfoSumOfProbesFmt, InfoSumOfProbes, PepsirfInfoSNPNFormat,
    PepsirfInfoSNPNDirFmt, InfoSNPN, EnrichThresh, EnrichThreshFileDirFmt,
    EnrichThreshFileFormat, SubjoinMultiFileFmt, SubjoinMultiFileDirFmt, MultiFile,
    ProteinFasta, ProteinFastaFmt, ProteinFastaDirFmt, PeptideFasta, PeptideFastaFmt,
    PeptideFastaDirFmt, Link, PepsirfLinkTSVFormat, PepsirfLinkTSVDirFmt,
    PeptideIDListDirFmt, PeptideIDList, PepsirfDMP, PepsirfDMPDirFmt, PepsirfDMPFormat,
    DeconvSingluar, PepsirfDeconvSingularFormat, PepsirfDeconvSingularDirFmt, ScorePerRound,
    ScorePerRoundFmt, ScorePerRoundDirFmt, PepsirfDeconvBatchDirFmt, DeconvBatch,
    PeptideAssignMapFormat, PeptideAssignMapDirFmt, PeptideAssignmentMap, DemuxFif,
    PepsirfDemuxFifDirFmt, PepsirfDemuxFifFmt, DemuxSampleList, PepsirfDemuxSampleListFmt,
    PepsirfDemuxSampleListDirFmt, DemuxIndex, PepsirfDemuxIndexDirFmt, PepsirfDemuxIndexFmt,
    DemuxLibrary, PepsirfDemuxLibraryDirFmt, PepsirfDemuxLibraryFmt, DemuxFastq, PepsirfDemuxFastqDirFmt,
    PepsirfDemuxFastqFmt, DemuxDiagnostic, PepsirfDemuxDiagnosticDirFmt, PepsirfDemuxDiagnosticFormat,
    ProteinAlignmentManifestFormat, ProteinAlignment, PeptideToProteinAlignmentFormat,
    ProteinAlignmentDirFormat, MutantReference, MutantReferenceFileFmt, MutantReferenceDirFmt, ProteinAlignmentFmt
    )
import q2_pepsirf.actions as actions
import q2_pepsirf.actions.zscore as zscore
import q2_pepsirf.actions.enrich as enrich
import q2_pepsirf.actions.info as info
import q2_pepsirf.actions.norm as norm
import q2_pepsirf.actions.bin as bin
import q2_pepsirf.actions.subjoin as subjoin
import q2_pepsirf.actions.link as link
import q2_pepsirf.actions.deconv as deconv
import q2_pepsirf.actions.demux as demux

from q2_types.feature_table import FeatureTable, BIOMV210DirFmt

# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.
plugin = Plugin("pepsirf", version=q2_pepsirf.__version__,
                website="https://github.com/LadnerLab/q2-pepsirf",
                description="Qiime2 Plug-in for the use of pepsirf within qiime.")

# Register the sematic types as formats
plugin.register_formats(PepsirfContingencyTSVFormat,
                        PepsirfContingencyTSVDirFmt,
                        PeptideIDListFmt,
                        EnrichedPeptideDirFmt,
                        ZscoreNanFormat,
                        ZscoreNanDirFmt, 
                        PeptideBinFormat, 
                        PeptideBinDirFmt,
                        PepsirfInfoSNPNDirFmt,
                        PepsirfInfoSNPNFormat,
                        PepsirfInfoSumOfProbesDirFmt,
                        PepsirfInfoSumOfProbesFmt,
                        EnrichThreshFileFormat,
                        EnrichThreshFileDirFmt,
                        SubjoinMultiFileFmt,
                        SubjoinMultiFileDirFmt,
                        ProteinFastaFmt,
                        ProteinFastaDirFmt,
                        PeptideFastaFmt,
                        PeptideFastaDirFmt,
                        PepsirfLinkTSVFormat,
                        PepsirfLinkTSVDirFmt,
                        PeptideIDListDirFmt,
                        PepsirfDMPFormat,
                        PepsirfDMPDirFmt,
                        PepsirfDeconvSingularFormat,
                        PepsirfDeconvSingularDirFmt,
                        ScorePerRoundFmt,
                        ScorePerRoundDirFmt,
                        PepsirfDeconvBatchDirFmt,
                        PeptideAssignMapFormat,
                        PeptideAssignMapDirFmt,
                        PepsirfDemuxFifFmt,
                        PepsirfDemuxFifDirFmt,
                        PepsirfDemuxSampleListFmt,
                        PepsirfDemuxSampleListDirFmt,
                        PepsirfDemuxIndexFmt,
                        PepsirfDemuxIndexDirFmt,
                        PepsirfDemuxLibraryFmt,
                        PepsirfDemuxLibraryDirFmt,
                        PepsirfDemuxFastqFmt,
                        PepsirfDemuxFastqDirFmt,
                        PepsirfDemuxDiagnosticFormat,
                        PepsirfDemuxDiagnosticDirFmt,
                        ProteinAlignmentManifestFormat,
                        PeptideToProteinAlignmentFormat,
                        ProteinAlignmentFmt,
                        ProteinAlignmentDirFormat,
                        MutantReferenceFileFmt,
                        MutantReferenceDirFmt
                        )

# register all semantic types
plugin.register_semantic_types(
        Normed, NormedDifference, NormedDiffRatio,
        NormedRatio, NormedSized, Zscore, RawCounts,
        PairwiseEnrichment)
plugin.register_semantic_type_to_format(
        FeatureTable[
                Normed | NormedDifference |
                NormedDiffRatio | NormedRatio 
                | NormedSized | Zscore | RawCounts],
        PepsirfContingencyTSVDirFmt)
plugin.register_semantic_type_to_format(
        PairwiseEnrichment,
        EnrichedPeptideDirFmt
)
plugin.register_semantic_type_to_format(
        ZscoreNan,
        ZscoreNanDirFmt
)
plugin.register_semantic_type_to_format(
        PeptideBins,
        PeptideBinDirFmt
)
plugin.register_semantic_type_to_format(
        InfoSNPN,
        PepsirfInfoSNPNDirFmt
)
plugin.register_semantic_type_to_format(
        InfoSumOfProbes,
        PepsirfInfoSumOfProbesDirFmt
)
plugin.register_semantic_type_to_format(
        EnrichThresh,
        EnrichedPeptideDirFmt
)
plugin.register_semantic_type_to_format(
        MultiFile,
        SubjoinMultiFileDirFmt
)
plugin.register_semantic_type_to_format(
        ProteinFasta,
        ProteinFastaDirFmt
)
plugin.register_semantic_type_to_format(
        PeptideFasta,
        PeptideFastaDirFmt
)
plugin.register_semantic_type_to_format(
        Link,
        PepsirfLinkTSVDirFmt
)
plugin.register_semantic_type_to_format(
        PeptideIDList,
        PeptideIDListDirFmt
)
plugin.register_semantic_type_to_format(
        PepsirfDMP,
        PepsirfDMPDirFmt
)
plugin.register_semantic_type_to_format(
        DeconvSingluar,
        PepsirfDeconvSingularDirFmt
)
plugin.register_semantic_type_to_format(
        ScorePerRound,
        ScorePerRoundDirFmt
)
plugin.register_semantic_type_to_format(
        DeconvBatch,
        PepsirfDeconvBatchDirFmt
)
plugin.register_semantic_type_to_format(
        PeptideAssignmentMap,
        PeptideAssignMapDirFmt
)
plugin.register_semantic_type_to_format(
        DemuxFif,
        PepsirfDemuxFifDirFmt
)
plugin.register_semantic_type_to_format(
        DemuxSampleList,
        PepsirfDemuxSampleListDirFmt
)
plugin.register_semantic_type_to_format(
        DemuxIndex,
        PepsirfDemuxIndexDirFmt
)
plugin.register_semantic_type_to_format(
        DemuxLibrary,
        PepsirfDemuxLibraryDirFmt
)
plugin.register_semantic_type_to_format(
        DemuxFastq,
        PepsirfDemuxFastqDirFmt
)
plugin.register_semantic_type_to_format(
        DemuxDiagnostic,
        PepsirfDemuxDiagnosticDirFmt
)
plugin.register_semantic_type_to_format(
        ProteinAlignment,
        ProteinAlignmentDirFormat
)
plugin.register_semantic_type_to_format(
        MutantReference,
        MutantReferenceDirFmt
)

# create a type map to change outputs dependent on str choice
T_approach, T_out = TypeMap ({
        Str%Choices("col_sum"): Normed,
        Str%Choices("diff"): NormedDifference,
        Str%Choices("diff_ratio"): NormedDiffRatio,
        Str%Choices("ratio"): NormedRatio,
        Str%Choices("size_factors"): NormedSized
})

# action set up for norm module
plugin.methods.register_function(
        function=norm.norm,
        inputs={
                'peptide_scores': FeatureTable[RawCounts | Normed],
                'negative_control': FeatureTable[RawCounts | Normed]
        },
        parameters={
                'normalize_approach': T_approach,
                'negative_id': Str,
                'negative_names': List[Str],
                'precision': Int % Range(0, None),
                'pepsirf_binary': Str,
                'outfile': Str
        },
        outputs=[
                ('qza_output', FeatureTable[T_out])
        ],
        input_descriptions={
                'peptide_scores': "Name of FeatureTable matrix file containing peptide scores. This file should be in the same format "
                                "as the output from the demux module.",
                'negative_control': "Name of FeatureTable matrix file containing data for sb samples."
        },
        parameter_descriptions={
                'normalize_approach': "'col_sum': Normalize the scores using a column-sum method. Output per peptide is the score per "
                                "million for the sample (i.e., summed across all peptides). "
                                "'size_factors': Normalize the scores using the size factors method (Anders and Huber 2010). "
                                "'diff': Normalize the scores using the difference method. For each peptide and sample, the difference "
                                "between the score and the respective peptide's mean score in the negative controls is determined. "
                                "'ratio': Normalize the scores using the ratio method. For each peptide and sample, the ratio of score "
                                "to the respective peptide's mean score in the negative controls is determined. "
                                "'diff_ratio': Normalize the scores using the difference-ratio method. For each peptide and sample, "
                                "the difference between the score and the respective peptide's mean score in the negative controls is "
                                "first determined. This difference is then divided by the respective peptide's mean score in the "
                                "negative controls.",
                'negative_id': "Optional approach for identifying negative controls. Provide a unique string at the start of all "
                                "negative control samples.",
                'negative_names': "Optional approach for identifying negative controls. "
                                "Space-separated list of negative control sample names.",
                'precision': "Output score precision. The scores written to the output will be output to this many decimal places.",
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF."
        },
        output_descriptions={
                'qza_output': "the FeatureTable (.qza) output based on the normalized approach given by user"
        },
        name='pepsirf norm module',
        description="Normalize raw count data with pepsirf's norm module"
)

# action set up for zscore module
plugin.methods.register_function(
        function=zscore.zscore,
        inputs={
                'scores': FeatureTable[Normed | RawCounts | NormedDifference | NormedDiffRatio],
                'bins': PeptideBins
        },
        parameters={
                'trim': Float % Range(0.0, 100.0),
                'hdi': Float % Range(0.0, None),
                'num_threads': Int % Range(1, None),
                'pepsirf_binary': Str,
                'outfile': Str
        },
        outputs=[
                ('zscore_output', FeatureTable[Zscore]),
                ('nan_report', ZscoreNan)
        ],
        input_descriptions={
                'scores': "Name of the file to use as input. Should be a score matrix in the "
                        "format as output by the demux and subjoin modules. Raw or normalized read counts can be used.",
                'bins': "Name of the file containing bins, one bin per line, as output by the bin module. Each bin contains a "
                        "tab-delimited list of peptide names."
        },
        parameter_descriptions={
                'trim': "Percentile of lowest and highest counts within a bin to ignore when calculating "
                        "the mean and standard deviation.",
                'hdi': "Alternative approach for discarding outliers prior to calculating mean and stdev. If provided, this "
                        "argument will override --trim, which trims evenly from both sides of the distribution. For --hdi, the "
                        "user should provide the high density interval to be used for calculation of mean and stdev. For "
                        "example, '--hdi 0.95' would instruct the program to utilize the 95% highest density interval (from each "
                        "bin) for these calculations.",
                'num_threads': "The number of threads to use for analyses.",
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF."
        },
        output_descriptions={
                'zscore_output': "Name for the output Z scores file. This file will be a FeatureTable[Zscore] with the same "
                               "dimensions as the input score file. Each peptide will be written with its z-score within each sample.",
                'nan_report': "Name of the file to write out information regarding peptides that are given a zscore of 'nan'. This "
                        "occurs when the mean score of a bin and the score of the focal peptide are both zero. This will be a "
                        "ZscoreNan file, with three columns per line. The first column will contain the name of the peptide, "
                        "the second will be the name of the sample, and the third will be the bin number of the probe. This bin "
                        "number corresponds to the line number in the bins file, within which the probe was found."
        },
        name='pepsirf zscore module',
        description="Calculate Z scores for each peptide in each sample with pepsirf's zscore module"
)

# action set up for enrich module
plugin.methods.register_function(
        function=enrich.enrich,
        inputs={
                'zscores': FeatureTable[Zscore],
                'raw_scores': FeatureTable[RawCounts],
                'thresh_file': EnrichThresh,
                'col_sum': FeatureTable[Normed]

        },
        parameters={
                'exact_z_thresh': Str,
                'raw_constraint': Int % Range(0, None),
                'enrichment_failure': Bool,
                'truncate': Bool,
                'pepsirf_binary': Str,
                'source': MetadataColumn[Categorical],
                'exact_cs_thresh': Str,
                'low_raw_reads': Bool,
                'outfile': Str
        },
        outputs=[
                ('dir_fmt_output', PairwiseEnrichment)
        ],
        input_descriptions={
                'zscores': "FeatureTable containing z scores of the normalized read counts. "
                        "Fist column header must be 'Sequence Name' as produced by pepsirf.",
                'col_sum': "FeatureTable containing the normalized read counts. "
                        "Fist column header must be 'Sequence Name' as produced by pepsirf.",
                'raw_scores': "This matrix must contain the raw counts for each Peptide for every sample of "
                        "interest. If included, '--raw_score_constraint' must also be specified.",
                'thresh_file': " The name of a tab-delimited file containing one tab-delimited matrix filename "
                        "and threshold(s), one per line. If providing more than z score matrix."
        },
        parameter_descriptions={
                'exact_z_thresh': "Exact z score thresholds.",
                'exact_cs_thresh': "Exact col-sum threholds.",
                'raw_constraint': "The minimum total raw count across all peptides for a sample to be "
                                "included in the analysis.This provides a way to impose a minimum read "
                                "count for a sample to be evaluated.",
                'enrichment_failure': "For each sample set that does not result in the generation of an enriched peptide file, a row of two "
                                "tab-delimited columns is provided: the first column contains the replicate names (comma-delimited) and "
                                "the second column provides the reason why the sample did not result in an enriched peptide file. "
                                "This file is output to the same directory as the enriched peptide files. The 'Reason' column will contain "
                                "one of the following: 'Raw read count threshold' or 'No enriched peptides'.",
                'truncate': "By default each filename in the output directory will include every replicate name joined by the "
                        "'join_on' value. Alternatively, if more than two replicates are being evaluated, then you may include "
                        "this flag to stop the filenames from including more than 3 samplenames in the output. When this flag is "
                        "used, the output names will be of the form 'A~B~C~1more', for example.",
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'source': "Metadata file containing all sample names and their source groups. "
                        "Used to create pairs tsv to run pepsirf enrich module.",
                'low_raw_reads': "By default samples with any replicates below the raw read threshold will be dropped when this flag is included, replicates with reads above the threshold will be kept.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF."
        },
        output_descriptions={
                'dir_fmt_output': "Directory formatted qza containing lists of enriched peptides"
        },
        name='pepsirf enrich module',
        description="Determines which peptides are enriched in samples with pepsirf's enrich module"
)

# action set up for infoSNPN module
plugin.methods.register_function(
        function=info.infoSNPN,
        inputs={
                'input': FeatureTable[
                        Normed | NormedDifference | NormedDiffRatio
                        | RawCounts | NormedDiffRatio | NormedSized],
        },
        parameters={
                'get': Str%Choices("samples", "probes"),
                'pepsirf_binary': Str,
                'outfile': Str
        },
        outputs=[
                ('snpn_output', InfoSNPN)
        ],
        input_descriptions={
                'input': "An input score matrix to gather information from."
        },
        parameter_descriptions={
                'get': "Specify weather you want to collect sample names or probe/peptide names",
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF."
        },
        output_descriptions={
                'snpn_output':"InfoSNPN file in the form of a file with no header, one sample name per line."
        },
        name='pepsirf info (samples/peptide names) module',
        description="Gathers the number of samples and peptides in the matrix using pepsirf's info modules"
)

# action set up for infoSumOfProbes module
plugin.methods.register_function(
        function=info.infoSumOfProbes,
        inputs={
                'input': FeatureTable[
                        Normed | NormedDifference | NormedDiffRatio
                        | RawCounts | NormedDiffRatio | NormedSized],
        },
        parameters={
                'pepsirf_binary': Str,
                'outfile': Str
        },
        outputs=[
                ('sum_of_probes_output', InfoSumOfProbes)
        ],
        input_descriptions={
                'input': "An input score matrix to gather information from."
        },
        parameter_descriptions={
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF."
        },
        output_descriptions={
                'sum_of_probes_output':"InfoSumOfProbes file, The first entry in each column will be the name of the "
                                "sample, and the second will be the sum of the peptide/probe scores for the sample."
        },
        name='pepsirf info (sum of probes) module',
        description="Gathers the sum of probes per sample in the matrix using pepsirf's info module"
)

# create a typemap to allow other input types when bool is true
T_norm, T_flag, T_out = TypeMap ({
        (Normed, Bool%Choices(False)): PeptideBins,
        (NormedDifference | NormedDiffRatio | NormedRatio | NormedSized, Bool%Choices(True)): PeptideBins
})

# action set up for bin module
plugin.methods.register_function(
        function=bin.bin,
        inputs={
                'scores': FeatureTable[T_norm],
        },
        parameters={
                'pepsirf_binary': Str,
                'bin_size': Int % Range(1, None),
                'round_to': Int % Range(0, None),
                'allow_other_normalization': T_flag,
                'outfile': Str
        },
        outputs=[
                ('bin_output', T_out)
        ],
        input_descriptions={
                'scores': "Input tab-delimited normalized score matrix file to use for peptide binning. "
                        "This matrix should only contain info for the negative control samples that "
                        "should be used to generate bins (see subjoin module for help generatinig input "
                        "matrix). Peptides with similar scores, summed across the negative controls, "
                        "will be binned together."
        },
        parameter_descriptions={
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'bin_size': "The minimum number of peptides that a bin must contain. If a bin would be "
                        "created with fewer than bin_size peptides, it will be combined with the next "
                        "bin until at least bin_size peptides are found.",
                'round_to': "The 'rounding factor' for the scores parsed from the score matrix prior to "
                        "binning. Scores found in the matrix will be rounded to the nearest 1/10^x for a "
                        "rounding factor x. For example, a rounding factor of 0 will result in rounding "
                        "to the nearest integer, while a rounding factor of 1 will result in rounding to "
                        "the nearest tenth.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF."
        },
        output_descriptions={
                'bin_output':"PeptideBins file that contains  one bin per line and each "
                        "bin will be a tab-delimited list of the names of the peptides in the bin."
        },
        name='pepsirf bin module',
        description="Creates groups of peptides with similar starting abundances with pepsirf's bin module"
)

# create a typemap to choose output type based on str choice
s_approach, s_out = TypeMap ({
        Str%Choices("raw"): RawCounts,
        Str%Choices("col_sum"): Normed,
        Str%Choices("diff"): NormedDifference,
        Str%Choices("diff_ratio"): NormedDiffRatio,
        Str%Choices("ratio"): NormedRatio,
        Str%Choices("size_factors"): NormedSized
})

# action set up for subjoin module
plugin.methods.register_function(
        function=subjoin.subjoin,
        inputs={
                'multi_file': MultiFile
        },
        parameters={
                'multi_file_input': List[Str],
                'pepsirf_binary': Str,
                'input_type': s_approach,
                'outfile': Str,
                'subjoin_input': Str, 
                'filter_peptide_names': Bool,
                'duplicate_evaluation': Str%Choices("include", "combine", "ignore")
        },
        outputs=[
                ('subjoin_output', FeatureTable[s_out])
        ],
        input_descriptions={
                'multi_file': "The name of a tab-delimited file containing score matrix and sample name "
                        "list filename pairs, one per line. Each of these pairs must be a score "
                        "matrix and a file containing the names of samples (or peptides, if "
                        "specified) to keep from the input the score matrix. The score matrix "
                        "should be of the format output by the demux module, with sample names on "
                        "the columns and peptide names on the rows. The namelist must have one "
                        "name per line, but can optionally have 2, if renaming samples in the "
                        "subjoin output. Optionally, a name list can be omitted if all samples "
                        "from the input matrix should be included in the output."
        },
        parameter_descriptions={
                'multi_file_input': "To use multiple name lists with multiple score matrices, include this argument "
                                "with a list of file names. Optionally, a name list can be omitted if all samples "
                                "from the input matrix should be included in the output. "
                                "Ex:(file1,file1_samples file2,file2_samples) or Ex:(file1 file2 file3).",
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'input_type': "Specify the type of file being inputted into subjoin, in order to produce the correct "
                        "output file format type. raw = FeatureTable[RawCounts], col_sum = FeatureTable[Normed], etc.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF.",
                'subjoin_input': " Comma-separated filenames (For example: "
                                "score_matrix.tsv,sample_names.txt ). Each of these pairs must be a score "
                                "matrix and a file containing the names of samples (or peptides, if "
                                "specified) to keep in the score matrix. The score matrix should be of "
                                "the format output by the demux module, with sample names on the columns "
                                "and peptide names on the rows. The namelist must have one name per line, "
                                "but can optionally have 2. If 2 tab-delimited names are included on one "
                                "line, the name in the first column should match the name in the input "
                                "matrix file, while the name in the second column will be output. "
                                "Therefore, this allows for the renaming of samples in the output.",
                'filter_peptide_names':"Flag to include if the name lists input to the input or multi_file "
                                        "options should be treated as peptide (i.e. row) names instead of sample "
                                        "(i.e. column) names. With the inclusion of this flag, the input files "
                                        "will be filtered on peptide names (rows) instead of sample names (column).",
                'duplicate_evaluation': "Defines what should be done when sample or peptide names are not unique "
                                        "across files being joined. Currently, three different duplicate "
                                        "evaluation strategies are available: - combine: Combine (with addition) the values "
                                        "associated with identical sample/peptide names from different files. "
                                        "- include: Include each duplicate, adding a suffix to the duplicate "
                                        "name detailing the file from which the sample came. - ignore: Ignore the possibility "
                                        "of duplicates. Behavior is undefined when duplicates are  encountered in this "
                                        "mode Therefore, this mode is not recommended."
        },
        output_descriptions={
                'subjoin_output':"Name for the output score matrix file. The output will be in the form of "
                                "the input, but with only the specified values (samplenames or peptides) "
                                "found in the namelists."
        },
        name='pepsirf subjoin module',
        description="Manipulate matrix files with pepsirf's subjoin module"
)

# action set up for link module
plugin.methods.register_function(
        function=link.link,
        inputs={
                'protein_file': ProteinFasta,
                'peptide_file': PeptideFasta
        },
        parameters={
                'pepsirf_binary': Str,
                'outfile': Str,
                'meta': Str,
                'kmer_size': Int,
                'kmer_redundancy_control': Bool,
        },
        outputs=[
                ('link_output', Link)
        ],
        input_descriptions={
                'protein_file': "Name of fasta file containing protein sequences of interest.",
                'peptide_file': " Name of fasta file containing aa peptides of interest. These will generally be "
                                "peptides that are contained in a particular assay."
        },
        parameter_descriptions={
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF.",
                'meta': "Optional method for providing taxonomic information for each protein contained "
                        "in '--protein_file'. Three comma-separated strings should be provided: 1) name "
                        "of tab-delimited metadata file, 2) header for column containing protein "
                        "sequence name and 3) header for column containing ID to be used in creating"
                        "the linkage map.",
                'kmer_size':"Kmer size to use when creating the linkage map.",
                'kmer_redundancy_control': "Optional modification to the way scores are calculated. If this flag is used, "
                                        "then instead of a peptide receiving one point for each kmer it shares with "
                                        "proteins of a given taxonomic group, it receives 1 / ( the number of times the "
                                        "kmer appears in the provided peptides ) points."
        },
        output_descriptions={
                'link_output':"Name of the file to which output is written. Output will be in the form of a "
                        "tab-delimited file with a header. Each entry/row will be of the form: "
                        "peptide_name TAB id:score,id:score, and so on. By default, 'score' is defined "
                        "as the number of shared kmers."
        },
        name='pepsirf link module',
        description="Defines linkages between taxonomic groups and peptides based on shared kmers with pepsirf's link module"
)

# shared inputs for deoncv singular mode and batch mode
deconv_shared_inputs = {
                'linked': Link,
                'id_name_map': PepsirfDMP
        }

# shared parameters for deoncv singular mode and batch mode
deconv_shared_parameters = {
                'pepsirf_binary': Str,
                'outfile': Str,
                'threshold': Int,
                'scoring_strategy': Str%Choices("summation", "integer", "fraction"),
                'score_filtering': Bool,
                'score_tie_threshold': Float,
                'score_overlap_threshold': Float,
                'single_threaded': Bool
        }

# shared input descriptions for deoncv singular mode and batch mode
deconv_shared_input_descript = {
                'linked': "Name of linkage map to be used for deconvolution. It should be in the "
                        "format output by the 'link' module.",
                'id_name_map': "Optional file containing mappings from taxonomic id to taxon name. This file "
                        "should be formatted like the file 'rankedlineage.dmp' from NCBI. It is "
                        "recommended to either use this file or a subset of this file that contains all "
                        "of the taxon ids linked to peptides of interest. If included, the output will "
                        "contain a column denoting the name of the species as well as the id."
        }
# shared parameter descriptions for deoncv singular mode and batch mode
deconv_shared_parameters_descript = {
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF.",
                'threshold': "Minimum score that a taxon must obtain in order to be included in the "
                        "deconvolution report.",
                'scoring_strategy':"Scoring strategies 'summation', 'integer', or 'fraction' can be specified. "
                                "By not including this flag, summation scoring will be used by default. "
                                "The --linked file passed must be of the form created by the link module. "
                                "This means a file of tab-delimited values, one per line. "
                                "Each line is of the form peptide_name TAB id:score,id:score, and so on. "
                                "An error will occurif input is not in this format. For summation scoring, the score "
                                "assigned to each peptide/ID pair is determined by the ':score' portion of the --linked file. "
                                "For example, assume a line in the --linked file looks like the following: peptide_1 TAB 123:4,543:8 "
                                "The IDs '123' and '543' will receive scores of 4 and 8 respectively. "
                                "For integer scoring, each ID receives a score of 1 for every enriched peptide to which "
                                "it is linked (':score' is ignored). For fractional scoring, the score is assigned to each "
                                "peptide/ID pair is defined by 1/n for each peptide, where n is the number of IDs to which a peptide is linked. "
                                "In this method of scoring peptides, a peptide with fewer linked IDs is worth more points..",
                'score_filtering': "Include this option if you want filtering to be done by the score of each taxon, rather than the count "
                                "of linked peptides. If used, any taxon with a score below '--threshold' will be removed from consideration, "
                                "even if it is the highest scoring taxon. Note that for integer "
                                "scoring, both score filtering and count filtering (default) "
                                "are the same. If this flag is not included, then any species whose count falls below '--threshold' will be "
                                "removed from consideration. Score filtering is best suited for the summation scoring method.",
                'score_tie_threshold': "Threshold for two species to be evaluated as a tie. Note that this value can be either an integer "
                                "or a ratio that is in (0,1). When provided as an integer this value dictates the difference in score that "
                                "is allowed for two taxa to be considered as potentially tied. For example, if this flag is provided with "
                                "the value 0, then two or more taxa must have the exact same score to be tied. If this flag is provided with "
                                "the value 4, then the difference between the scores of two taxa must be no "
                                "greater than 4 to be considered tied. For example, if taxon 1 has a score of 5, and taxon 2 has a score "
                                "anywhere between the integer values in [1,9], then these species will be considered tied, and their tie "
                                "will be evaluated as dictated by the specified '--score_overlap_threshold'. If the argument provided to "
                                "this flag is in (0, 1), then the score for a taxon must be at least this proportion of the score for the "
                                "highest scoring taxon, to trigger a tie. So if species 1 has the highest score with a score of 9, and "
                                "species 2 has a score of 5, then this flag must be provided with value >= 4/5 = 0.8 for the species 1 and "
                                "2 to be considered tied. Note that any values provided to this flag that are in the set { x: x >= 1 } - Z, "
                                "where Z is the set of integers, will result in an error. So 4.45 is not a valid value, "
                                "but both 4 and 0.45 are.",
                'score_overlap_threshold': "Once two species have been determined to be tied, according to '--score_tie_threshold', "
                                        "they are then evaluated as a tie. To use integer tie evaluation, where species must share an "
                                        "integer number of peptides, not a ratio of their total peptides, provide this argument with a "
                                        "value in the interval [1, inf). For ratio tie evaluation, which is used when this argument is "
                                        "provided with a value in the interval (0,1), two taxon must reciprocally share at least the "
                                        "specified proportion of peptides to be reported together. For example, suppose species 1 shares "
                                        "half (0.5) of its peptides with species 2, but species 2 only shares a tenth (0.1) of its "
                                        "peptides with species 1. These two will only be reported together if score_overlap_threshold' <= 0.1.",
                'single_threaded': "By default this module uses two threads. Include this option with no arguments "
                                "if you only want only one thread to be used."
        }

# action set up for deconv_singular module
plugin.methods.register_function(
        function=deconv.deconv_singular,
        inputs={
                'enriched': PeptideIDList,
                **deconv_shared_inputs
        },
        parameters=deconv_shared_parameters,
        outputs=[
                ('deconv_output', DeconvSingluar),
                ('score_per_round', ScorePerRound)
        ],
        input_descriptions={
                'enriched': "single file containing the names of enriched peptides, one per line. "
                        "Each peptide contained in this file (or files) should have a corresponding "
                        "entry in the '--linked' input file.",
                **deconv_shared_input_descript
        },
        parameter_descriptions=deconv_shared_parameters_descript,
        output_descriptions={
                'deconv_output':"Name of the file to which output is written. Output will be in the form of a "
                        "tab-delimited file with a header.",
                'score_per_round':"Name of directory to write counts/scores to after every round. If "
                                "included, the counts and scores for all remaining taxa will be recorded after "
                                "every round. Filenames will be written in the format '$dir/round_x', where x "
                                "is the round number. The original scores will be written to "
                                "'$dir/round_0'. A new file will be written to the directory after each "
                                "subsequent round. If this flag is included and the specified directory "
                                "exists, the program will exit with an error."
        },
        name='pepsirf deconv singular mode module',
        description="converts a list of enriched peptides into a parsimony-based list of likely taxa to which the assayed individual has likely "
                "been exposed with pepsirf's deconv singular mode module"
)

# action set up for deconv_batch module
plugin.methods.register_function(
        function=deconv.deconv_batch,
        inputs={
                'enriched_dir': PairwiseEnrichment,
                **deconv_shared_inputs
        },
        parameters={
                'outfile_suffix': Str,
                'mapfile_suffix': Str,
                'remove_file_types': Bool,
                **deconv_shared_parameters
        },
        outputs=[
                ('deconv_output', DeconvBatch),
                ('score_per_round', ScorePerRound),
                ('peptide_assignment_map', PeptideAssignmentMap)
        ],
        input_descriptions={
                'enriched_dir': "Name of a directory containing files, that contain the names of enriched "
                               "peptides, one per line. Each Peptide contained within these files should have "
                               "a corresponding entry in the '--linked' input file.", 
                **deconv_shared_input_descript
        },
        parameter_descriptions={
                'outfile_suffix': "Used for batch mode only. When specified, the name of each file "
                                "written to the output  directory will have this suffix.",
                'mapfile_suffix': "Used for batch mode only. When specified, the name of each "
                                "'--peptide_assignment_map' will have this suffix.",
                'remove_file_types': "Use this flag to exclude input file ('--enrich') extensions from the names "
                                "of output files. Not used in singular mode.",
                **deconv_shared_parameters_descript
        },
        output_descriptions={
                'deconv_output':"Name of the file to which output is written. Output will be in the form of a "
                        "tab-delimited file with a header.",
                'score_per_round':"Name of directory to write counts/scores to after every round. If "
                                "included, the counts and scores for all remaining taxa will be recorded after "
                                "every round. Filenames will be written in the format '$dir/round_x', where x "
                                "is the round number. The original scores will be written to "
                                "'$dir/round_0'. A new file will be written to the directory after each "
                                "subsequent round. If this flag is included and the specified directory "
                                "exists, the program will exit with an error.",
                'peptide_assignment_map': "Optional output. If specified, a map detailing which peptides were assigned "
                                       "to which taxa will be written. If this module is run in batch mode, this will "
                                       "be used as a directory name for the peptide maps to be stored. Maps will be "
                                        "tab-delimited files with the first column being peptide names; the second "
                                       "column containing a comma-separated list of taxa to which the peptide was "
                                       "assigned; the third column will be a list of the taxa with which the peptide "
                                       "originally shared a kmer. Note that the second column will only contain "
                                       "multiple values in the event of a tie."
        },
        name='pepsirf deconv batch mode module',
        description="converts a list of enriched peptides into a parsimony-based list of likely taxa to which the assayed individual has likely "
                "been exposed with pepsirf's deconv batch mode module"
)

# action set up for demux module
plugin.methods.register_function(
        function=demux.demux,
        inputs={
                'input_r1': DemuxFastq,
                'input_r2': DemuxFastq,
                'index': DemuxIndex,
                'samplelist': DemuxSampleList,
                'fif': DemuxFif,
                'library': DemuxLibrary
        },
        parameters={
                'seq': Str,
                'read_per_loop': Int,
                'num_threads': Int,
                'phred_base': Int,
                'phred_min_score': Int,
                'sindex': Str,
                'translate_aggregates': Bool,
                'concatemer': Bool,
                'sname': Str,
                'index1': Str,
                'index2': Str,
                'pepsirf_binary': Str,
                'outfile': Str
        },
        outputs=[
                ('raw_counts_output', FeatureTable[RawCounts]),
                ('diagnostic_output', DemuxDiagnostic)
        ],
        input_descriptions={
                'input_r1': " Fastq-formatted file containing reads with DNA tags. If PepSIRF was NOT "
                        "compiled with Zlib support, this file must be uncompressed. If PepSIRF was "
                        "compiled with Zlib support, then this file can be uncompressed or compressed "
                        "using gzip. In this case, the file format will be automatically "
                        "determined.",
                'input_r2': "Optional index-only fastq file. If PepSIRF was NOT compiled with Zlib "
                        "support, this file must be uncompressed. If PepSIRF was compiled "
                        "with Zlib support, then this file can be uncompressed or compressed using "
                        "gzip. In this case, the file format will be automatically determined. "
                        "Note that if this argument is not supplied, only 'index1' will be used to "
                        "identify samples.",
                'index': "Name of fasta-formatted file containing forward and (potentially) reverse index "
                        "sequences. Sequence names must match exactly with those supplied in the "
                        "'samplelist'.",
                'samplelist': "A tab-delimited list of samples with a header row and one sample per line. "
                        "This file must contain at least one index column and one sample name "
                        "column. Multiple index columns may be included. This file can also include "
                        "additional columns that will not be used for the demultiplexing. Specify "
                        "which columns to use with the '--sname', '--sindex1', and '--sindex2' "
                        "flags. If '-fif' is used, then only '-sname' will be used.",
                'fif': "The flexible index file can be provided as an alternative to the '--index1' and "
                "'--index2' options. The file must use the following format: a tab-delimited "
                "file with 5 ordered columns: 1) index name, which should correspond to a "
                "header name in the sample sheet, 2) read name, which should be either 'r1' "
                "or 'r2' (not case-sensitive) to specify whether the index is in '--input_r1' or "
                "'--input_r2', 3) index start location (0-based, inclusive), 4) index length "
                "and 5) number of mismatched to allow. '--index1', '--index2', '--sname', "
                "'--sindex1', and 'sindex2' will be ignored if this option is provided.",
                'library': "Fasta-formatted file containing reference DNA tags. If this flag is not "
                        "included, reference-independent demultiplexing will be performed. In "
                        "reference-independent mode, each sequence in the region specified by "
                        "'--seq' will be considered its own reference, and the observed sequences "
                        "will be used as the row names in the output count matrix."
        },
        parameter_descriptions={
                'seq': "Positional information for the DNA tags. This argument must be passed in "
                "the same format specified for 'index1'.",
                'read_per_loop': "The number of fastq records read a time. A higher value will result in "
                                "more memory usage by the program, but will also result in fewer disk "
                                "accesses, increasing performance of the program.",
                'num_threads': "Number of threads to use for analyses.",
                'phred_base': "Phred base to use when parsing fastq quality scores. Valid options include "
                        "33 or 64.",
                'phred_min_score': "The minimum average phred-scaled quality score for the DNA tag portion "
                                "of a read for it to be considered for matching. This means that if the "
                                "average phred33/64 score for a read at the expected locations of the DNA tag "
                                "is not at least this then the read will be discarded.",
                'sindex': "Used to specify the header for the index 1 and optional index 2 column in "
                        "the samplelist. This is an alternative to using the '--fif'' option.",
                'translate_aggregates': "Include this flag to use translation-based aggregation. In this "
                                "mode, counts for nt sequences will be combined if they translate into the "
                                "same aa sequence. Note: When this mode is used, the name of the aggregate "
                                "sequence will be the sequence that was a result of the translation. Therefore, "
                                "this mode is most appropriate for use with reference-independent "
                                "demultiplexing.",
                'concatemer': "Concatenated adapter/primer sequences (optional). The presence of this "
                        "sequence within a read indicates that the expected DNA tag is not present. If "
                        "supplied, the number of times this concatemer is recorded in the input "
                        "file is reported.",
                'sname': "Used to specify the header for the sample name column in the samplelist. "
                        "By default 'SampleName' is set as the column header name.",
                'index1': "Positional information for index1 (i.e barcode 1). This argument must be "
                        "passed as 3 comma-separated values. The first item represents the (0-based) "
                        "expected start position of the first index; the second represents the length "
                        "of the first index; and the third represents the number of mismatches "
                        "that are tolerated for this index. An example is '--index1 12,12,1'. This "
                        "says that the index starts at (0-based) position 12, the index is 12 "
                        "nucleotides long, and if a perfect match is not found, then up to one "
                        "mismatch will be tolerated.",
                'index2': "Positional information for index2, optional. This argument must be passed "
                        "in the same format specified for '--index1'. If '--input2' is provided, "
                        "this positional information is assummed to refer to the reads contained in this "
                        "second, index-only fastq file. If '--input_r2' is NOT provided, this "
                        "positional information is assumed to refer to the reads contained in the "
                        "'--input_r1' fastq file.",
                'pepsirf_binary': "The binary to call pepsirf on your system.",
                'outfile': "The outfile that will produce a list of inputs to PepSIRF."
        },
        output_descriptions=None,
        name='pepsirf demux module',
        description="takes the following parameters and outputs counts for each reference"
                "sequence (i.e. probe/peptide) for each sample with pepsirf's demux module "
                "(MUST precompile pepsirf's develop branch to run this module)"
)

# import all sematic type transformers
importlib.import_module("q2_pepsirf.transformers")
