library(dplyr)
context("Test tag.search()")
test_that("Test input dataset structure and search results feedbacks", {
  ## error feedback when no "group" column is in input dataset
  d = tibble(groupname = c("control","control", "lab1", "lab1"),
             mass = c(100, 200, 150, 250))
  expect_error(tag.search(d, delta = 50),
               regexp = "A variable/column named `group` is required")

  ## error feedback when no "mass" column is in input dataset
  d = tibble(group = c("control","control", "lab1", "lab1"),
             `m/z` = c(100, 200, 150, 250))
  expect_error(tag.search(d, delta = 50),
               regexp = "A variable/column named `mass` is required")

  ## error feedback when no "control" level is in "group" column
  d = tibble(group = c("ctrl","ctrl", "lab1", "lab1"),
             mass = c(100, 200, 150, 250))
  expect_error(tag.search(d, delta = 50),
               regexp = "`control` is required in the `group` column.")

  ## error feedback when no input in tag.search delta argument
  d = tibble(group = c("control","control", "lab1", "lab1"),
             mass = c(100, 200, 150, 250))
  expect_error(tag.search(d),
               regexp = "Mass shift input required.")


  ## check if successful remove duplicated rows (for group and mass only)
  d = tibble(group = c("control","control", "lab1", "lab1"),
             mass = c(100, 100, 150, 250),
             intensity = c(10, 20, 30, 40),
             peptide = c("a", "b", "c", "d"))
  expected.input = d %>% select(group, mass) %>% distinct() %>% left_join(d, by = c("group", "mass"))
  expected.output = tag.search(expected.input, delta = 50)[[1]]

  expect_equal(tag.search(d, delta = 50)[[1]], expected.output)


  ## check if augment input with "intensity" if it is lacking
  d = tibble(group = c("control","control", "lab1", "lab1"),
             mass = c(100, 100, 150, 250))
  tag.search(d, delta = 50)
  expect_output(tag.search(d, delta = 50)[[1]]$intensity)
  expect_equal(tag.search(d, delta = 50)[[1]]$intensity %>% class(), "numeric")

})




test_that("Test cases when no search found", {

  ## check feedback when no pair is found
  d = tibble(group = c("control","control", "lab1", "lab1"),
             mass = c(100, 180, 150, 250),
             intensity = c(500, 600, 200, 300))
  tag.search(d, delta = 40)
  expect_message(tag.search(d, delta = 40))


  # check feedback when no match is found
  d = tibble(group = c("control","control", "lab1", "lab1"),
             mass = c(100, 180, 150, 250))
  expect_equal(tag.search(d, delta = 50)[[2]],
               "Found paired peaks (mass differentiate by expected delta) , but not matched peaks (of the same mass).")


  ## check feedback when both pair and  match is found
  d = tibble(group = c("control","control", "lab1", "lab1"),
             mass = c(100, 150, 150, 250))
  expect_equal(tag.search(d, delta = 50)[[2]],
               "Found both paired peaks (mass differentiate by expected delta) and matched peaks (of the same mass).")


  ## check run when no Mismatch is found
  d = tibble(group = c("control","control",
                       "lab1", "lab1",
                       "lab2", "lab2", "lab2"),
             mass = c(100, 150,
                      150, 250,
                      150, 200, 250))
  d %>% tag.search(delta = c(50, 100))
  expect_output(d %>% tag.search(delta = c(50, 100)))

})



