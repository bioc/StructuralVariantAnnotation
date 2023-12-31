context('util functions')
expect_equal(elementExtract(CharacterList(c("a"), c("b", "c"), c(), c("d", "e", "f"), "g")), c("a", "b", NA_character_, "d", "g"))
expect_equal(elementExtract(CharacterList(c("b", "c"))), c("b"))
expect_equal(elementExtract(CharacterList()), character(0))
expect_equal(elementExtract(CharacterList(c())), c(NA_character_))
expect_equal(elementExtract(c("aa", "bb", "aa")), c("aa", "bb", "aa"))
expect_equal(elementExtract(c("aa", "bb", "aa"), 3), c(NA_character_, NA_character_, NA_character_))

expect_equal(elementExtract(DNAStringSetList(c("A"), c("C", "G"), c())), c("A", "C", NA_character_))
expect_equal(elementExtract(DNAStringSet(c("A", "C", "N"))), c("A", "C", "N"))

expect_equal(elementExtract(IntegerList(c(), c(1, NA), c(2, 3), c(4, 5))), c(NA, 1, 2, 4))
expect_equal(elementExtract(IntegerList(c(), c(1, NA), c(2, 3), c(4, 5)), 2), c(NA, NA, 3, 5))
expect_equal(elementExtract(IntegerList()), integer(0))

expect_equal(elementExtract(NULL), NULL)

expect_equal(.replaceNa(NULL, c(1,2,3)), c(1,2,3))
expect_equal(.replaceNa(c(1,NA,NA), c(3,2,NA)), c(1,2,NA))
expect_equal(.replaceNa(c(1,2,3), NULL), c(1,2,3))
expect_equal(.replaceNa(c(1,NA), c(1,2,3)), c(1,2))
expect_equal(.replaceNa(c(1,NA), integer(0)), c(1,NA))
# want to be able to perform operations on VCF columns that may not exist
expect_equal(.replaceNa(numeric(0), c(1,2)), c(1,2))
expect_equal(.replaceNa((NULL - c(2, NA)), c(1,2)), c(1,2))
expect_equal(c(1,2,7,3), c(1,2,NA,3) |> .replaceNa(c()) |> .replaceNa(c(5,6,7,8)))

expect_equal(.pairwiseLCPrefix(DNAStringSet(c("Aa", "CCCt", "N")), c("aA", "ccct", ""), ignore.case=TRUE), c(2, 4, 0))
expect_equal(.pairwiseLCPrefix(DNAStringSet(c("ACGT")), c("CGT"), ignore.case=TRUE), c(0))

# extensions-DataFrame.R
# expect_equal(.as.matrix.list(CharacterList(c("a"), c("b", "c"), c(), c("d", "e", "f"), "g")),
#              matrix(c("a", NA, NA,
#                       "b", "c", NA,
#                       NA, NA, NA,
#                       "d", "e", "f",
#                       "g", NA, NA), ncol=3, byrow=TRUE))
# expect_equal(.as.matrix.list(IntegerList(c(1), c(2), c(3))),
#              matrix(c(1,2,3), ncol=1, byrow=TRUE))
# expect_equal(.as.matrix.list(IntegerList(c(1,4), c(2,5), c(3,6))),
#              matrix(c(1,4,2,5,3,6), ncol=2, byrow=TRUE))
# expect_equal(.as.matrix.list(IntegerList(c(), c(2,5), c())),
#              matrix(c(NA, NA, 2, 5, NA, NA), ncol=2, byrow=TRUE))
# expect_equal(.as.matrix.list(IntegerList(c(), c(), c())),
#              matrix(integer(0), ncol=0, nrow=3))
# expect_equal(.as.matrix.list(CharacterList(c(), c(), c())),
#              matrix(character(0), ncol=0, nrow=3))


