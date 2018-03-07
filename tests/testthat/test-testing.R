library(RSeqLib)

library(testthat)

context("test-testing.R")

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("BWA works", {
    bwa <- BWA(seq = "CACTAGCTAGCTACGCGGGGGCGCGCGCGCGCGAAAAACACTTTCACAG")
    expect_message(show(bwa), "RSeqLib BWA object with params mc.cores = 1, hardclip = 0, keep_sec_with_frac_of_primary_score = 0.9, max_secondary = 10")
})

test_that("BWA_method works", {
    expect_equal(BWA_method("object"), "BWA__object")    
})

test_that("query works", {
    bwa <- BWA(seq = "CACTAGCTAGCTACGCGGGGGCGCGCGCGCGCGAAAAACACTTTCACAG")
    expect_equal(typeof(query(bwa, c("CACTAGCTAGCTACGCGGGGGCGCG", "CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG"))), "S4")
})
