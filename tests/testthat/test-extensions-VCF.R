context('main functions')
example <- readVcf(.testfile("vcf4.2.example.sv.vcf"), "")
simple <- readVcf(.testfile("simple.vcf"), "")
breakend <- readVcf(.testfile("breakend.vcf"), "")
multipleAlleles <- readVcf(.testfile("multipleAltSVs.vcf"), "")
representations <- readVcf(.testfile("representations.vcf"))


breakdancer <- readVcf(.testfile("breakdancer-1.4.5.vcf"), "")
#cortex <- readVcf(.testfile("cortex-1.0.5.14.vcf"), "")
#crest <- readVcf(.testfile("crest.vcf"), "")
delly <- readVcf(.testfile("delly-0.6.8.vcf"), "")
#gasv <- readVcf(.testfile("gasv-20140228.vcf"), "")
gridss <- readVcf(.testfile("gridss.vcf"), "")
gridss_missing_partner <- readVcf(.testfile("gridss-missingPartner.vcf"), "")
#lumpy <- readVcf(.testfile("lumpy-0.2.11.vcf"), "")
pindel <- readVcf(.testfile("pindel-0.2.5b6.vcf"), "")
#socrates <- readVcf(.testfile("socrates-1.13.vcf"), "")
tigra <- readVcf(.testfile("tigra-0.3.7.vcf"), "")
manta <- readVcf(.testfile("manta-0.29.6.vcf"), "")
manta111 <- readVcf(.testfile("manta-1.1.1.vcf"), "")
longranger <- readVcf(.testfile("COLO829.10X.longranger.largeSVs.vcf"), "")

test_that("INFO column import", {
	gr <- breakpointRanges(simple, info_columns=c("SVTYPE", "MATEID"))
	expect_equal("BNDBF", as.character(gr["BNDFB"]$MATEID))

	gr <- breakpointRanges(.testrecord(c("chr10	2991435	INV	N	<INV>	.	LowQual	SVTYPE=INV;CHR2=chr1;END=19357517;CT=3to5")))
})
test_that("longranger UNK", {
	gr <- breakpointRanges(longranger)[c("call_2416_1", "call_2416_2")]
	expect_equal(as.character(strand(gr)), c("*", "*"))
})
test_that("longranger IMPRECISE_DIR", {
	gr <- breakpointRanges(longranger)[c("call_534_bp1", "call_534_bp2")]
	expect_equal(c(2632546-10, 2632545-10+57120), start(gr))
	expect_equal(as.character(strand(gr)), c("*", "*"))
})
test_that("Delly TRA", {
	# https://groups.google.com/forum/#!msg/delly-users/6Mq2juBraRY/BjmMrBh3GAAJ
	# Sorry, I forgot to post this to the delly-users list:
	# For a translocation, you have 2 double strand breaks, one on chrA and one on chrB.
	# This creates 4 "dangling" ends, chrA_left, chrA_right, chrB_left, chrB_right.
	# For a translocation you can join chrA_left with chrB_left (3to3), chrA_left with chrB_right (3to5),
	# chrA_right with chrB_left (5to3) and chrA_right with chrB_right (5to5).
	# In fact for a typical reciprocal translocation in prostate cancer (where two chromosomes exchange their end)
	# Delly calls 2 translocations at the breakpoint, one 3to5 and one 5to3. But obviously not all translocations are reciprocal.
	# -Tobias
    gr <- breakpointRanges(.testrecord(c("chr10	2991435	TRA00000001	N	<TRA>	.	LowQual	CIEND=0,100;CIPOS=0,50;SVTYPE=TRA;CHR2=chr1;END=19357517;CT=3to5")))
    expect_equal(2, length(gr))
    expect_equal(c(2991435, 19357517), start(gr))
    expect_equal(c(2991485, 19357617), end(gr))
    expect_equal(c("+", "-"), as.character(strand(gr)))
    expect_equal(c("chr10", "chr1"), as.character(seqnames(gr)))
    gr <- breakpointRanges(delly)
})
test_that("tigra CTX", {
	gr <- breakpointRanges(tigra[3,])
	expect_equal(c(102520604 - 10, 70284 - 10), start(gr))
	expect_equal(c(102520604 + 10, 70284 + 10), end(gr))
	expect_equal(c("*", "*"), as.character(strand(gr)))
    expect_equal(c("chr12", "chrUn_gl000223"), as.character(seqnames(gr)))
})
test_that("pindel RPL", {
	gr <- breakpointRanges(pindel[5,])
	expect_equal(c(1029142, 1029183), start(gr))
	expect_equal(c("+", "-"), as.character(strand(gr)))
	expect_equal(c(51, 51), gr$insLen)
})
test_that("empty VCF", {
	expect_equal(0, length(breakpointRanges(.testrecord(c()))))
})
test_that(".hasSingleAllelePerRecord", {
	expect_false(.hasSingleAllelePerRecord(multipleAlleles))
	expect_true(.hasSingleAllelePerRecord(expand(multipleAlleles)))
})
test_that("isSymbolic", {
	expect_equal(
	    c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	    isSymbolic(simple))
})
test_that("isStructural", {
	expect_equal(
	    c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	    isStructural(simple))
	expect_true(isStructural(.testrecord("chr1	1	.	ATT	AGGA	.	.	")))
	expect_false(isStructural(.testrecord("chr1	1	.	ATT	NNN	.	.	")))
})
test_that(".svLen", {
	expect_equal(c(0, 0, 1, -1, 1, -2, NA, NA, NA, NA), .svLen(simple))
	expect_equal(.svLen(.testrecord("chr1	100	.	A	<DEL>	.	.	SVLEN=-1")), c(-1))
	expect_equal(.svLen(.testrecord("chr1	100	.	A	<DUP>	.	.	END=101")), c(1))
})
# unpaired breakend
test_that("partner fails if missing mate", {
	expect_error(partner(breakpointRanges(simple)[1,]))
})
test_that("breakpointRanges convert to breakend pairs", {
	gr <- breakpointRanges(simple)
	pairId <- c("INS", "DEL", "SYMINS", "SYMDEL")
	expect_true(all(paste0(pairId, "_bp1") %in% names(gr)))
	expect_true(all(paste0(pairId, "_bp2") %in% names(gr)))
	expect_equal(names(partner(gr))[names(gr) %in% paste0(pairId, "_bp1")], paste0(pairId, "_bp2"))
	expect_equal(names(partner(gr))[names(gr) %in% c("BNDFB", "BNDBF")], c("BNDBF", "BNDFB"))
})
test_that("breakpointRanges creates placeholder names", {
	expect_warning(expect_named(breakpointRanges(.testrecord(c(
		"chr1	100	.	A	<DEL>	.	.	SVLEN=-1",
		"chr1	100	.	A	<DEL>	.	.	SVLEN=-1")))),
		regex="Found 1 duplicate row names")
})
test_that("breakpointRanges non-symbolic alleles", {
	gr <- breakpointRanges(simple[c("INS", "DEL"),])
	expect_equal(4, length(gr))

	gr <- breakpointRanges(.testrecord("chr1	1	.	ATT	AGGA	.	.	"))
	expect_equal(start(gr), c(1, 4))
	expect_equal(gr$insSeq, c("GGA", "GGA"))
	expect_equal(gr$svLen, c(1, 1))

	gr <- breakpointRanges(.testrecord("chr1	2	.	TTT	AGGA	.	.	"))
	expect_equal(start(gr), c(1, 5))
	expect_equal(gr$insSeq, c("AGGA", "AGGA"))
	expect_equal(gr$svLen, c(1, 1))

	gr <- breakpointRanges(.testrecord("chr1	2	.	AGT	AGGA	.	.	"))
	expect_equal(start(gr), c(3, 5))
	expect_equal(gr$insSeq, c("GA", "GA"))
	expect_equal(gr$insLen, c(2, 2))
	expect_equal(gr$svLen, c(1, 1))

	gr <- breakpointRanges(.testrecord("chr1	2	.	AGGA	AG	.	.	"))
	expect_equal(start(gr), c(3, 6))
	expect_equal(as.character(strand(gr)), c("+", "-"))
	expect_equal(gr$insLen, c(0, 0))
	expect_equal(gr$svLen, c(-2, -2))
})
test_that("breakpointRanges intervals", {
	# Position assumed to the left aligned
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=0"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(100, 101))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=1"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(101, 102))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=2"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(102, 103))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMSEQ=AAAAAAAAAA"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(110, 111))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	CIPOS=-5,10"))
	expect_equal(start(gr), c(95, 96))
	expect_equal(end(gr), c(110, 111))

	# CIPOS over homology
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	CIPOS=-5,10;HOMLEN=50;HOMESEQ=A"))
	expect_equal(start(gr), c(95, 96))
	expect_equal(end(gr), c(110, 111))

	# HOMLEN over HOMSEQ
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=50;HOMESEQ=A"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(150, 151))
})
test_that("breakpointRanges DEL", {
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	SVLEN=-1"))
	expect_equal(start(gr), c(100, 102))
	expect_equal(gr$insLen, c(0, 0))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	END=101"))
	expect_equal(start(gr), c(100, 102))

	breakpointRanges(.testrecord("chr1	100	.	A	<DEL:WITH:SUBTYPE>	.	.	END=101"))
	breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	SVTYPE=DEL:WITH:SUBTYPE;END=101"))

	# warning about incompatable SVLEN and END fields
	#expect_warning(breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	END=101;SVLEN=-10")), "SVLEN")
})
test_that("breakpointRanges should fix positive DEL event size", {
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	SVLEN=10"))
	expect_equal(start(gr), c(100, 111))
	expect_equal(gr$insLen, c(0, 0))
})
test_that("breakpointRanges breakend", {
	expect_warning(gr <- breakpointRanges(breakend))
	expect_equal("parid_b", gr["parid_a",]$partner)
	expect_equal("mateid_b", gr["mateid_a",]$partner)
	expect_equal(partner(gr[c("parid_a", "parid_b"),]), gr[c("parid_b", "parid_a"),])
	expect_warning(breakpointRanges(breakend[c("mateid_a", "mateid_b", "multi_mateid")]), "Ignoring additional mate breakends")
	expect_warning(breakpointRanges(breakend[c("unpaired")]), "Removing [0-9]+ unpaired breakend variants")
	expect_equal(breakpointRanges(.testrecord(c(
		"chr1	100	a	N	N[chr1:105[	.	.	SVTYPE=BND;CIPOS=0,1;PARID=b",
		"chr1	105	b	N	]chr1:100]N	.	.	SVTYPE=BND;CIPOS=0,1;PARID=a"
	)))$svLen, c(-4, -4))
	expect_equal(breakpointRanges(.testrecord(c(
		"chr1	100	a	N	NAAAA[chr1:101[	.	.	SVTYPE=BND;CIPOS=0,1;PARID=b",
		"chr1	101	b	N	]chr1:100]TTTTN	.	.	SVTYPE=BND;CIPOS=0,1;PARID=a"
	)))$svLen, c(4, 4))
})
test_that("breakpointRanges INV", {
	# VCF example
	gr <- breakpointRanges(.testrecord("chr1	321682	INV0	T	<INV>	6	PASS	SVTYPE=INV;END=421681"))
	expect_equal(4, length(gr))
	expect_equal(start(gr), c(321682+1, 421681+1, 321682, 421681))
	expect_equal(as.character(strand(gr)), c("-", "-", "+", "+"))
	expect_equal(names(gr), c("INV0_bp1", "INV0_bp2", "INV0_bp3", "INV0_bp4"))

	gr <- breakpointRanges(.testrecord("chr1	321682	INV0	T	<INV>	6	PASS	SVTYPE=INV;END=421681;CIPOS=-2,1;CIEND=-3,4"))
	expect_equal(4, length(gr))
	expect_equal(start(gr), c(321682+1, 421681+1, 321682, 421681) + c(-2, -3, -2 ,-3))
	expect_equal(  end(gr), c(321682+1, 421681+1, 321682, 421681) + c( 1,  4,  1,  4))

	expect_error(breakpointRanges(.testrecord("chr1	321682	INV0	T	<INV>	6	PASS	SVTYPE=INV")))
})
test_that("breakpointRanges DUP", {
	# VCF example
	gr <- breakpointRanges(.testrecord(c(
		"chr1	12665100	.	A	<DUP>	14	PASS	SVTYPE=DUP;END=12686200;SVLEN=21100",
		"chr1	18665128	.	T	<DUP:TANDEM>	11	PASS	SVTYPE=DUP;END=18665204;SVLEN=76")))
	expect_equal(4, length(gr))
	expect_equal(start(gr), c(12665100+1, 18665128+1,
							  12686200, 18665204))
	expect_equal(as.character(strand(gr)), c("-", "-", "+", "+"))
	expect_equal(c(0,0,0,0), gr$insLen)

	gr <- breakpointRanges(.testrecord("chr1	12665100	.	A	<DUP>	14	PASS	SVTYPE=DUP;END=12686200;SVLEN=21100;CIPOS=-2,1;CIEND=-3,4"))
	expect_equal(2, length(gr))
	expect_equal(start(gr), c(12665100+1, 12686200) + c(-2,-3))
	expect_equal(  end(gr), c(12665100+1, 12686200) + c( 1, 4))

	expect_error(breakpointRanges(.testrecord("chr1	321682	.	T	<DUP>	.	.	SVTYPE=DUP")))

	gr <- breakpointRanges(.testrecord("chr12	5616362	chr12.5616362.DUP65536	A	<DUP>	.	.	SVLEN=65536;SVTYPE=DUP"))
	expect_equal(2, length(gr))
	expect_equal(start(gr), c(5616362+1, 5616362+65536))
	expect_equal(c("-", "+"), as.character(strand(gr)))
})

# test_that("manta merge should retain only unique events", {
# 	# VCF example
# 	gr <- breakpointRanges(manta)
# 	expect_equal(4, length(gr))
# })
test_that("manta 1.1.1", {
	gr <- breakpointRanges(manta111)
	expect_equal(2*10, length(gr))
})
test_that("manta INV3 should only have 1 breakpoint", {
	gr <- breakpointRanges(manta111[info(manta111)$INV3])
	expect_equal(c("+", "+"), as.character(strand(gr)))
})

test_that("nominalPosition should ignore confidence intervals", {
	# VCF example
	vcfExact <- .testrecord(c("chr1	100	.	A	<DEL>	14	PASS	SVTYPE=DEL;END=200"))
	vcfCI <- .testrecord(c("chr1	100	.	A	<DEL>	14	PASS	SVTYPE=DEL;END=200;CIEND=-10,10;CIPOS=-5,5"))
	expect_equal(breakpointRanges(vcfCI, nominalPosition=TRUE), breakpointRanges(vcfExact, nominalPosition=TRUE))
	expect_equal(ranges(breakpointRanges(vcfCI, nominalPosition=TRUE)), ranges(breakpointRanges(vcfExact, nominalPosition=TRUE)))
})
test_that("nominalPosition should ignore micro-homology", {
	# VCF example
	vcfExact <- .testrecord(c("chr1	100	.	A	<DEL>	14	PASS	SVTYPE=DEL;END=200"))
	vcfHOM <- .testrecord(c("chr1	100	.	A	<DEL>	14	PASS	SVTYPE=DEL;END=200;HOMLEN=15"))
	expect_equal(ranges(breakpointRanges(vcfHOM, nominalPosition=TRUE)), ranges(breakpointRanges(vcfExact, nominalPosition=TRUE)))
})
test_that("breakpointRanges should not include breakends", {
	expect_true(isSymbolic(.testrecord(c("chr1	100	.	A	AAA.	14	PASS	SVTYPE=BND"))))
	expect_equal(0, length(breakpointRanges(.testrecord(c("chr1	100	.	A	AAA.	14	PASS	SVTYPE=BND")))))
	expect_equal(0, length(breakpointRanges(.testrecord(c("chr1	1541062	gridss0_2065b	G	.GGTGGG	262.69	.	SVTYPE=BND")))))
})
test_that("breakendRanges should include breakends", {
	gr <- breakendRanges(.testrecord(c("chr1	100	.	A	TGC.	14	PASS	SVTYPE=BND")))
	expect_equal(1, length(gr))
	expect_equal("GC", gr$insSeq)
})
test_that("align_breakpoint should handle all orientations", {
	noname = function(x) { names(x) = NULL; return(x)}
	vcf = .testrecord(c(
		"chr1	1000	a	N	A[chr1:2000[	.	.	SVTYPE=BND;CIPOS=0,4;PARID=b",
		"chr1	2000	b	N	]chr1:1000]G	.	.	SVTYPE=BND;CIPOS=0,4;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(1002, 2002), start(rowRanges(vcf)))
	expect_equal(c(-2,2, -2,2), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("N[chr1:2002[", "]chr1:1002]N"), unlist(rowRanges(vcf)$ALT))
	expect_equal(c("+-", "-+"), .vcfAltToStrandPair(rowRanges(vcf)$ALT))

	vcf = .testrecord(c(
		"chr1	1000	a	N	]chr1:2000]A	.	.	SVTYPE=BND;CIPOS=0,4;PARID=b",
		"chr1	2000	b	N	G[chr1:1000[	.	.	SVTYPE=BND;CIPOS=0,4;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(1002, 2002), start(rowRanges(vcf)))
	expect_equal(c(-2,2, -2,2), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("]chr1:2002]N", "N[chr1:1002["), unlist(rowRanges(vcf)$ALT))
	expect_equal(c("-+", "+-"), .vcfAltToStrandPair(rowRanges(vcf)$ALT))

	vcf = .testrecord(c(
		"chr1	1000	a	N	A]chr1:2000]	.	.	SVTYPE=BND;CIPOS=-4,0;PARID=b",
		"chr1	2000	b	N	G]chr1:1000]	.	.	SVTYPE=BND;CIPOS=0,4;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(998, 2002), start(rowRanges(vcf)))
	expect_equal(c(-2,2, -2,2), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("N]chr1:2002]", "N]chr1:998]"), unlist(rowRanges(vcf)$ALT))
	expect_equal(c("++", "++"), .vcfAltToStrandPair(rowRanges(vcf)$ALT))

	vcf = .testrecord(c(
		"chr1	1000	a	N	[chr1:2000[A	.	.	SVTYPE=BND;CIPOS=-4,0;PARID=b",
		"chr1	2000	b	N	[chr1:1000[G	.	.	SVTYPE=BND;CIPOS=0,4;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(998, 2002), start(rowRanges(vcf)))
	expect_equal(c(-2,2, -2,2), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("[chr1:2002[N", "[chr1:998[N"), unlist(rowRanges(vcf)$ALT))
	expect_equal(c("--", "--"), .vcfAltToStrandPair(rowRanges(vcf)$ALT))
})
test_that("align_breakpoint centre should ensure odd length homology is consistent on both sides", {
	noname = function(x) { names(x) = NULL; return(x)}
	vcf = .testrecord(c(
		"chr1	1000	a	N	A[chr1:2000[	.	.	SVTYPE=BND;CIPOS=0,5;PARID=b",
		"chr1	2000	b	N	]chr1:1000]G	.	.	SVTYPE=BND;CIPOS=0,5;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(1002, 2002), start(rowRanges(vcf)))
	expect_equal(c(-2,3, -2,3), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("N[chr1:2002[", "]chr1:1002]N"), unlist(rowRanges(vcf)$ALT))

	vcf = .testrecord(c(
		"chr1	1000	a	N	]chr1:2000]A	.	.	SVTYPE=BND;CIPOS=0,5;PARID=b",
		"chr1	2000	b	N	G[chr1:1000[	.	.	SVTYPE=BND;CIPOS=0,5;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(1002, 2002), start(rowRanges(vcf)))
	expect_equal(c(-2,3, -2,3), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("]chr1:2002]N", "N[chr1:1002["), unlist(rowRanges(vcf)$ALT))

	vcf = .testrecord(c(
		"chr1	1000	a	N	A]chr1:2000]	.	.	SVTYPE=BND;CIPOS=-5,0;PARID=b",
		"chr1	2000	b	N	G]chr1:1000]	.	.	SVTYPE=BND;CIPOS=0,5;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(998, 2002), start(rowRanges(vcf)))
	expect_equal(c(-3,2, -2,3), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("N]chr1:2002]", "N]chr1:998]"), unlist(rowRanges(vcf)$ALT))

	vcf = .testrecord(c(
		"chr1	1000	a	N	[chr1:2000[A	.	.	SVTYPE=BND;CIPOS=-5,0;PARID=b",
		"chr1	2000	b	N	[chr1:1000[G	.	.	SVTYPE=BND;CIPOS=0,5;PARID=a"))
	vcf = align_breakpoints(vcf)
	expect_equal(c(998, 2002), start(rowRanges(vcf)))
	expect_equal(c(-3,2, -2,3), noname(unlist(info(vcf)$CIPOS)))
	expect_equal(c("[chr1:2002[N", "[chr1:998[N"), unlist(rowRanges(vcf)$ALT))
})
test_that("align_breakpoint should not touch other variants", {
	vcf = .testrecord(c(
		"chr1	1000	be1	N	AGT.	.	.	SVTYPE=BND;CIPOS=-5,0",
		"chr1	1000	b21	N	.AGT	.	.	SVTYPE=BND;CIPOS=-5,0",
		"chr1	1000	b21	N	<DEL>	.	.	SVTYPE=DEL;CIPOS=-5,0"))
	vcf = align_breakpoints(vcf)
	expect_equal(c("AGT.", ".AGT", "<DEL>"), unlist(rowRanges(vcf)$ALT))
})
test_that("breakpointRanges() should default to drop unpaired records.", {
	gr = expect_warning(breakpointRanges(gridss_missing_partner))
	expect_equal(2, length(gr))
})
test_that("breakpointRanges(inferMissingBreakends=TRUE) should add missing breakends.", {
	gr = breakpointRanges(gridss_missing_partner, inferMissingBreakends=TRUE)
	expect_equal(4, length(gr))
	expect_equal(c(18992158, 84963533, 84350, 4886681), start(gr))
	expect_equal(c("+", "-"), as.character(strand(gr))[1:2])
})

# CGTGTtgtagtaCCGTAA
#      -------       7bp del
# 0        1
# 123456789012345678

# Important symbolic allele info from the VCF specifications:
# 
# 1.6.1.4 If any of the ALT alleles is a symbolic allele (an angle-bracketed ID String “<ID>”) then the padding base is required and POS denotes the coordinate of the base preceding the polymorphism. 
# 1.6.1.8 End reference position (1-based), indicating the variant spans positions POS–END on reference/contig CHROM. Normally this is the position of the last base in the REF allele
test_that("representations are equivalent", {
    bpgr = breakpointRanges(representations)
    expect_equal(rep(-7, 6), bpgr$svLen)
    expect_equal(rep(c(5, 13), 3), start(bpgr))
    expect_equal(rep(c(5, 13), 3), end(bpgr))
})

# TODO: VCFv4.4 support
# - filter DEL/DUP on SVCLAIM
# - Basic support for CNV?
# - Deprecate SVTYPE

# VCFv4.4 PR: CIEND deprecate - use CILEN
    # need to define how CILEN interacts with CIPOS - there are multiple possible interpretation
    # [CIPOS_start, CIPOS_end], [CIEND_start, CIEND_end] implies CILEN=[CIPOS_end-CIEND_start, CIEND_end, CIPOS_start]
    # but the actual CILEN could be less. Homology
    # bonus: usually doesn't need to be written!
# VCFv4.4 PR: HOMLEN deprecate and replace with HOMPOS
    # Exact calls shouldn't have a CIPOS - they should use HOMPOS instead.
# VCFv4.4 PR: Event/Type should be type=A, not type=1

# VCFv4.4 PR: line 196: ALT=<BND> is not valid

# VCFv4.4 PR: update the \subsection{Encoding Structural Variants} example
	# The SV example should not include non-reserved subtypes

# VCFv4.4 PR: line 208: reserve all IUPAC ambiguity codes
# VCFv4.4 PR: line 208: can symoblic alleles be mixed with non-symbolic ones?
# VCFv4.4 TODO: why is the FORMAT CN field TYPE=1 ? Is it just the CN? Why Integer, not Float?

test_that("VCFv4.4 use ALT instead of SVTYPE", {
	vcf44 = .testrecordv44(c(
		# we'll ignore that the mate orientations are actually incorrect for now - we're just testing missing SVTYPE
		"chrA	1000	bp1	N	A]chrA:2000]	.	.	PARID=bp2",
		"chrA	2000	bp2	N	]chrA:1000]G	.	.	PARID=bp1",
		"chrA	1000	bp4	N	[chrA:2000[A	.	.	PARID=bp4",
		"chrA	2000	bp3	N	G[chrA:1000[	.	.	PARID=bp3",
		"chrA	2000	be1	N	.G	.	.	",
		"chrA	2000	be2	N	.G	.	.	",
		"chrA	1	del	A	<DEL>	0	.	SVLEN=10;SVCLAIM=J",
		"chrA	1	dup	A	<DUP>	0	.	SVLEN=10;SVCLAIM=J",
		"chrA	1	ins	A	<INS>	0	.	SVLEN=10",
		"chrA	1	inv	A	<INV>	0	.	SVLEN=10"))
	bpgr = breakpointRanges(vcf44)
	begr = breakendRanges(vcf44)
	expect_equal(sort(c(
		"bp1", "bp2", "bp3", "bp4",
		"del_bp1", "dup_bp1", "ins_bp1", "inv_bp1",
		"del_bp2", "dup_bp2", "ins_bp2", "inv_bp2",
		"inv_bp3", "inv_bp4")), sort(names(bpgr)))
	expect_equal(sort(c("be1", "be2")), sort(names(begr)))
})
test_that("VCFv4.4 CNV & SVCLAIM", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	1	cnv_del	A	<DEL>	0	.	SVLEN=10;SVCLAIM=D",
		"chrA	1	cnv_dup	A	<DUP>	0	.	SVLEN=10;SVCLAIM=D",
		"chrA	1	cnv	A	<CNV>	0	.	SVLEN=10",
		"chrA	1	del_actual	A	<DEL>	0	.	SVLEN=10;SVCLAIM=DJ",
		"chrA	1	dup_actual	A	<DUP>	0	.	SVLEN=10;SVCLAIM=DJ",
		"chrA	1	del	A	<DEL>	0	.	SVLEN=10;SVCLAIM=J",
		"chrA	1	dup	A	<DUP>	0	.	SVLEN=10;SVCLAIM=J",
		"chrA	1	ins	A	<INS>	0	.	SVLEN=10")))
	expected_names = c("del_actual", "dup_actual", "del", "dup", "ins")
	expect_equal(c(paste0(expected_names, "_bp1"), paste0(expected_names, "_bp2")), names(bpgr))
})
test_that("VCFv4.4 SVLEN>END", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	1	end	A	<DEL>	0	.	SVCLAIM=J;END=3",
		"chrA	1	svlen	A	<DEL>	0	.	SVCLAIM=J;SVLEN=5",
		"chrA	1	both	A	<DEL>	0	.	SVCLAIM=J;END=10;SVLEN=5")))
	expect_equal(4, start(bpgr["end_bp2"]))
	expect_equal(7, start(bpgr["svlen_bp2"]))
	expect_equal(7, start(bpgr["both_bp2"]))
})
test_that("VCFv4.4 CIEND>CILEN DEL", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	10	both1	A	<DEL>	0	.	SVLEN=5;END=15;CIPOS=-2,2;CILEN=-2,2;CIEND=0,0",
		"chrA	10	both2	A	<DEL>	0	.	SVLEN=5;END=15;CIPOS=-2,2;CILEN=0,0;CIEND=-2,2",
		"chrA	10	cilen	A	<DEL>	0	.	SVLEN=7;CIPOS=-2,2;CILEN=-3,3", #END=17
		"chrA	10	ciend	A	<DEL>	0	.	END=14;CIPOS=-2,2;CIEND=-1,3",
		"chrA	10	inferred_end_bounds	A	<DEL>	0	.	SVLEN=7;CIPOS=-2,2")))
	expect_equal(c(10,16) + c(-2,  0), start(bpgr[c("both1_bp1", "both1_bp2")]))
	expect_equal(c(10,16) + c( 2,  0),   end(bpgr[c("both1_bp1", "both1_bp2")]))
	
	expect_equal(c(10,16) + c(-2, -2), start(bpgr[c("both2_bp1", "both2_bp2")]))
	expect_equal(c(10,16) + c( 2,  2),   end(bpgr[c("both2_bp1", "both2_bp2")]))
	
	# widest end bounds is [leftmost start & shortest, rightmost start + longest]
	expect_equal(c(10,18) + c(-2, -2-3), start(bpgr[c("cilen_bp1", "cilen_bp2")]))
	expect_equal(c(10,18) + c( 2, +2+3),   end(bpgr[c("cilen_bp1", "cilen_bp2")]))
	
	expect_equal(c(10,15) + c(-2, -1), start(bpgr[c("ciend_bp1", "ciend_bp2")]))
	expect_equal(c(10,15) + c( 2,  3),   end(bpgr[c("ciend_bp1", "ciend_bp2")]))
		
	# inferred end bounds should match start since SVLEN is known and fixed
	expect_equal(c(10,18) + c(-2, -2), start(bpgr[c("inferred_end_bounds_bp1", "inferred_end_bounds_bp2")]))
	expect_equal(c(10,18) + c( 2,  2),   end(bpgr[c("inferred_end_bounds_bp1", "inferred_end_bounds_bp2")]))
})
test_that("VCFv4.4 CIEND CILEN INV", {
	# inversion has a different codepath
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	10	inv1	A	<INV>	0	.	SVLEN=5;CIPOS=-1,1;CILEN=-2,2;CIEND=0,0",
		"chrA	10	inv2	A	<INV>	0	.	SVLEN=5;CIPOS=-1,1;CILEN=-2,2")))
	expect_equal(c(16, 15) + 0,   start(bpgr[c("inv1_bp2", "inv1_bp4")]))
	expect_equal(c(16, 15) + 0,     end(bpgr[c("inv1_bp2", "inv1_bp4")]))
	expect_equal(c(16, 15) + -1-2,start(bpgr[c("inv2_bp2", "inv2_bp4")]))
	expect_equal(c(16, 15) + 1+2,   end(bpgr[c("inv2_bp2", "inv2_bp4")]))
})
test_that("VCFv4.4 CIEND CILEN INS", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	10	inv1	A	<INS>	0	.	SVLEN=5;CIPOS=-1,1;CILEN=-2,2;CIEND=0,0",
		"chrA	10	inv2	A	<INS>	0	.	SVLEN=5;CIPOS=-1,1;CILEN=-2,2")))
	# CILEN doesn't impact INS bounds
	expect_equal(11 + c(0,-1),start(bpgr[c("inv1_bp2", "inv2_bp2")]))
	expect_equal(11 + c(0, 1),  end(bpgr[c("inv1_bp2", "inv2_bp2")]))
})
test_that("VCFv4.4 ALT symbolic subtypes", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	1	ins	A	<INS:ME:ALU>	0	.	SVLEN=5",
		"chrA	1	dup	A	<DUP:TANDEM>	0	.	SVLEN=3"
	)))
	expect_equal(1, start(bpgr["ins_bp1"]))
	expect_equal(1,   end(bpgr["ins_bp1"]))
	expect_equal(2, start(bpgr["ins_bp2"]))
	expect_equal(2,   end(bpgr["ins_bp2"]))
	expect_equal(2, start(bpgr["dup_bp1"]))
	expect_equal(4, start(bpgr["dup_bp2"]))
	expect_equal("-", as.character(strand(bpgr["dup_bp1"])))
	expect_equal("+", as.character(strand(bpgr["dup_bp2"])))
})
test_that("VCFv4.4 EVENT", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	1	e1	A	<DEL>	0	.	SVLEN=5;SVCLAIM=J;EVENT=complex;EVENTTYPE=chromothripsis",
		"chrA	10	e2	A	<DEL>	0	.	SVLEN=5;SVCLAIM=J;EVENT=complex;EVENTTYPE=chromothripsis",
		"chrA	20	untagged	A	<DEL>	0	.	SVLEN=5;SVCLAIM=J")))
	expect_equal(c("complex", "complex", NA_character_), bpgr[c("e1_bp1", "e2_bp1", "untagged_bp1")]$event)
})
test_that("VCFv4.4 non-SV symbolic alleles", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	1	sym1	A	<*>	0	.	END=10",
		"chrA	10	sym2	A	<NON_REF>	0	.	")))
	expect_equal(0, length(bpgr))
})
test_that("SVLEN", {
	bpgr = breakpointRanges(.testrecordv44(c(
		"chrA	2	ins	A	<INS>	0	.	SVLEN=3",
		"chrA	2	dup	A	<DUP>	0	.	SVLEN=3",
		"chrA	2	del	A	<DEL>	0	.	SVLEN=3",
		"chrA	2	inv	A	<INV>	0	.	SVLEN=3",
		"chrA	2	cnv	A	<CNV>	0	.	SVLEN=3")))
	expect_equal(end(bpgr)-start(bpgr), rep(0, length(bpgr)))
	expect_equal(c(2, 3), start(bpgr[c("ins_bp1", "ins_bp2")]))
	expect_equal(c(3, 5), start(bpgr[c("dup_bp1", "dup_bp2")]))
	expect_equal(c(2, 6), start(bpgr[c("del_bp1", "del_bp2")]))
	expect_equal(c(3, 6, 2, 5), start(bpgr[c("inv_bp1", "inv_bp2", "inv_bp3", "inv_bp4")]))
	
	expect_equal(c("+", "-"), as.character(strand(bpgr[c("ins_bp1", "ins_bp2")])))
	expect_equal(c("-", "+"), as.character(strand(bpgr[c("dup_bp1", "dup_bp2")])))
	expect_equal(c("+", "-"), as.character(strand(bpgr[c("del_bp1", "del_bp2")])))
	expect_equal(c("-", "-", "+", "+"), as.character(strand(bpgr[c("inv_bp1", "inv_bp2", "inv_bp3", "inv_bp4")])))
})
# test_that("VCFv4.4 support multiple SV ALT alleles", {
# 	bpgr = breakpointRanges(.testrecordv44(c(
# 		"chrA	10	multi	A	<DEL>,<DUP>	0	.	SVLEN=5,10;CIPOS=-1,1,-2,2;SVCLAIM=J")))
# 	expect_equal(4, length(bpgr))
# 	expect_equal(c("multi_bp1", "multi_bp2", "multi_alt2_bp1", "multi_alt2_bp2"), names(bpgr))
# 	expect_equal(c(10-1, 10+5-1), start(bpgr[c("multi_bp1", "multi_bp2")]))
# 	expect_equal(c(10+1, 10+5+1),   end(bpgr[c("multi_bp1", "multi_bp2")]))
# 	expect_equal(c(10-2, 10+10-2), start(bpgr[c("multi_alt2_bp1", "multi_alt2_bp2")]))
# 	expect_equal(c(10+2, 10+10+2),   end(bpgr[c("multi_alt2_bp1", "multi_alt2_bp2")]))
# })






