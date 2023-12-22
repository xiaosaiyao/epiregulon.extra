library(vdiffr)
AR <- data.frame(ID = c("ANDROGEN RESPONSE","PROLIFERATION","MAPK"),
p.adjust = c(0.001, 0.01, 0.04))
GATA6 <- data.frame(ID = c("STK33","PROLIFERATION","MAPK"), p.adjust = c(0.001, 0.01, 0.04))
enrichresults <- list(AR = AR, GATA6 = GATA6)


test_that("plotGseaNetwork works correctly",{
  expect_doppelganger(
    title = "Network plot",
    fig = plotGseaNetwork(tf = names(enrichresults), enrichresults = enrichresults)

  )
})
