#' Testthat helper utility to locate files used
#' for package tests
#' @param filename Name of the test file.
#' @param location Directory of the test file.
#' @return Returns the file to be tested.
#' @noRd
.testfile <- function(filename, location="extdata") {
    if (file.exists(filename)) return(filename)
    f <- system.file(location, filename, package="StructuralVariantAnnotation")
    if (!file.exists(f)) {
        f <- file.path(getwd(), "inst", location, filename)
    }
    assertthat::assert_that(file.exists(f))
    return(f)
}
#' Loading a VCF containing the given records
#' @param record string vector of record to write
#' @return A VCF object.
#' @noRd
.testrecord <- function(record) {
    filename=tempfile(fileext=".vcf")
    write(paste0(c(
        "##fileformat=VCFv4.2",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=INV,Description=\"Inversion\">",
        "##ALT=<ID=TRA,Description=\"Translocation\">",
        "##ALT=<ID=INS,Description=\"Insertion\">",
        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">",
        "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">",
        "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">",
        "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">",
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=PARID,Number=1,Type=String,Description=\"ID of partner breakend\">",
		"##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">",
    	"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    	"##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">",
    	"##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">",
    	"##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">",
        "##contig=<ID=chr1,length=249250621>",
		"##contig=<ID=chr2,length=1000000000>",
		"##contig=<ID=chr3,length=1000000000>",
		"##contig=<ID=chr4,length=1000000000>",
		"##contig=<ID=chr5,length=1000000000>",
		"##contig=<ID=chr6,length=1000000000>",
		"##contig=<ID=chr7,length=1000000000>",
		"##contig=<ID=chr8,length=1000000000>",
		"##contig=<ID=chr9,length=1000000000>",
		"##contig=<ID=chr10,length=1000000000>",
		"##contig=<ID=chr11,length=1000000000>",
    	"##contig=<ID=chr12,length=1000000000>",
    	"##contig=<ID=chrM,length=16571>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        record), collapse="\n"),
        file=filename)
    vcf <- readVcf(.testfile(filename), "")
	# readVcf holds a file handle open so we won't be able to delete the
	# file even if we tried
    #file.remove(filename)
    return(vcf)
}
#' Loading a VCF containing the given records
#' @param record string vector of record to write
#' @return A VCF object.
#' @noRd
.testrecordv44 <- function(record) {
	filename=tempfile(fileext=".vcf")
	write(paste0(c(
		"##fileformat=VCFv4.4",
		"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">",
		"##INFO=<ID=NOVEL,Number=0,Type=Flag,Description=\"Indicates a novel structural variation\">",
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
		"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of the structural variant\">",
		"##INFO=<ID=CIPOS,Number=.,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">",
		"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">",
		"##INFO=<ID=HOMLEN,Number=A,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">",
		"##INFO=<ID=HOMSEQ,Number=A,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">",
		"##INFO=<ID=BKPTID,Number=A,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">",
		"##INFO=<ID=MEINFO,Number=.,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">",
		"##INFO=<ID=METRANS,Number=.,Type=String,Description=\"Mobile element transduction info of the form CHR,START,END,POLARITY\">",
		"##INFO=<ID=DGVID,Number=A,Type=String,Description=\"ID of this element in Database of Genomic Variation\">",
		"##INFO=<ID=DBVARID,Number=A,Type=String,Description=\"ID of this element in DBVAR\">",
		"##INFO=<ID=DBRIPID,Number=A,Type=String,Description=\"ID of this element in DBRIP\">",
		"##INFO=<ID=MATEID,Number=A,Type=String,Description=\"ID of mate breakends\">",
		"##INFO=<ID=PARID,Number=A,Type=String,Description=\"ID of partner breakend\">",
		"##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of associated event\">",
		"##INFO=<ID=EVENTTYPE,Number=1,Type=String,Description=\"Type of associated event\">",
		"##INFO=<ID=CILEN,Number=.,Type=Integer,Description=\"Confidence interval for the SVLEN field\">",
		"##INFO=<ID=DP,Number=A,Type=Integer,Description=\"Read Depth of segment containing breakend\">",
		"##INFO=<ID=DPADJ,Number=A,Type=Integer,Description=\"Read Depth of adjacency\">",
		"##INFO=<ID=CN,Number=A,Type=Integer,Description=\"Copy number of segment containing breakend\">",
		"##INFO=<ID=CNADJ,Number=A,Type=Integer,Description=\"Copy number of adjacency\">",
		"##INFO=<ID=CICN,Number=.,Type=Integer,Description=\"Confidence interval around copy number for the segment\">",
		"##INFO=<ID=CICNADJ,Number=.,Type=Integer,Description=\"Confidence interval around copy number for the adjacency\">",
		"##INFO=<ID=SVCLAIM,Number=A,Type=String,Description=\"Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively.\">",
		"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">",
		"##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">",
		"##FORMAT=<ID=CNL,Number=G,Type=Float,Description=\"Copy number genotype likelihood for imprecise events\">",
		"##FORMAT=<ID=CNP,Number=G,Type=Float,Description=\"Copy number posterior probabilities\">",
		"##FORMAT=<ID=NQ,Number=1,Type=Integer,Description=\"Phred style probability score that the variant is novel\">",
		"##FORMAT=<ID=HAP,Number=1,Type=Integer,Description=\"Unique haplotype identifier\">",
		"##FORMAT=<ID=AHAP,Number=1,Type=Integer,Description=\"Unique identifier of ancestral haplotype\">",
		"##ALT=<ID=DEL,Description=\"Deletion\">",
		"##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">",
		"##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">",
		"##ALT=<ID=DUP,Description=\"Duplication\">",
		"##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">",
		"##ALT=<ID=INS,Description=\"Insertion of novel sequence\">",
		"##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">",
		"##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">",
		"##ALT=<ID=INV,Description=\"Inversion\">",
		"##ALT=<ID=CNV,Description=\"Copy number variable region\">",
		"##contig=<ID=chrA,length=10000>",
		"##contig=<ID=chrB,length=10000>",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
		record), collapse="\n"),
		file=filename)
	vcf <- readVcf(.testfile(filename), "")
	# readVcf holds a file handle open so we won't be able to delete the
	# file even if we tried
	#file.remove(filename)
	return(vcf)
}