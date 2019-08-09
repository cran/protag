context("Test mass spectra drawing functions: listplot and butterfly plot")
test_that("Test various cases of matching results",{

  ## test with 3 plots, NO MATCH FOUND
  d = tibble(group = c("control","control",
                       "lab1", "lab1",
                       "lab2", "lab2", "lab2"),
             mass = c(100, 150,
                      128, 365,
                      245, 345, 111))

  # manual check output
  tag.search(d, delta = 28) %>% tag.spectra.listplot()
  tag.search(d, delta = 28) %>% tag.spectra.butterflyplot() # note that "warnings" was cat output
  # autocheck
  expect_output(tag.search(d, delta = 28) %>% tag.spectra.listplot())
  expect_output(tag.search(d, delta = 28) %>% tag.spectra.butterflyplot())


  ## test with 2 plots, NO MATCH FOUND
  d = tibble(group = c("control","control",
                       "lab2", "lab2", "lab2"),
             mass = c(100, 150,
                      128, 245, 345)); d

  # manual check
  tag.search(d, delta = 28) %>% tag.spectra.listplot()
  tag.search(d, delta = 28) %>% tag.spectra.butterflyplot() # note that "warnings" was cat output
  # autocheck
  expect_output(tag.search(d, delta = 28) %>% tag.spectra.listplot())
  expect_output(tag.search(d, delta = 28) %>% tag.spectra.butterflyplot())


  ## test with 3 plots, NO MIIIIIISMATCH FOUND!
  d = tibble(group = c("control","control",
                       "lab1", "lab1",
                       "lab2", "lab2", "lab2"),
             mass = c(100, 300,
                      128, 300,
                      156, 328, 356)); d
  # manual check
  tag.search(d, delta = c(28, 56)) %>% tag.spectra.butterflyplot() # note that "warning" was cat output
  tag.search(d, delta = c(28, 56)) %>% tag.spectra.listplot()
  # autocheck
  expect_output( tag.search(d, delta = c(28, 56)) %>% tag.spectra.butterflyplot() )
  expect_output(tag.search(d, delta = c(28, 56)) %>% tag.spectra.listplot())



  ## test with 2 plots, NO MATCH, NO MISMATCH FOUND
  d = tibble(group = c("control","control",
                       "lab1", "lab1",
                       "lab2", "lab2", "lab2"),
             mass = c(100, 300,
                      128, 356,
                      128, 156, 356)); d
  # manual check
  tag.search(d, delta = c(28, 56)) %>% tag.spectra.butterflyplot()
  tag.search(d, delta = c(28, 56)) %>% tag.spectra.listplot()
  tag.search(d %>% filter(group %in% c("control", "lab1")),
             delta = c(28, 56)) %>% tag.spectra.butterflyplot()
  # autocheck
  expect_output(tag.search(d, delta = c(28, 56)) %>% tag.spectra.butterflyplot() )
  expect_output(tag.search(d, delta = c(28, 56)) %>% tag.spectra.listplot())
  expect_output(tag.search(d %>% filter(group %in% c("control", "lab1")),
                           delta = c(28, 56)) %>% tag.spectra.listplot())
})


