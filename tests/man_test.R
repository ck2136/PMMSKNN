# shinytest install
devtools::install_github("rstudio/shinytest")
library(shinytest)
shinytest::installDependencies()
# require phantomjs... now run recordTest() to record interaction with the application
recordTest(".")

library(testthat)
test_that("Applicationworks", {
  expect_pass(testApp("/srv/shiny-server/prc/",
                      testnames = "mytest",
                      compareImages = FALSE))
})

testApp("/srv/shiny-server/prc/", 
        testnames = "mytest",
        compareImages = FALSE)
