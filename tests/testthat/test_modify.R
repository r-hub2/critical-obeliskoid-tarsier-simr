context("modify")

mod1 <- lmer(y ~ x + (1|g), data=simdata)
mod2 <- lmer(y ~ (1|g) + (1|x), data=simdata) # note that x has more groups, listed first by lme4

test_that("errors are thrown", {

    expect_error(fixef(fm1)["z"] <- 3, " is not the name of a fixed effect.")

    expect_error(sigma(fm2) <- 8, "sigma is not applicable for this model.")
    expect_error(sigma(fm2) <- 1, NA)

    expect_error(scale(fm2) <- 5, "scale is not applicable for this model.")

    expect_error(sigma(fglm) <- 8, "sigma is not applicable for this model.")
    expect_error(sigma(fglm) <- NULL, NA)
})

test_that("scale<- modifies VarCorr", {

    scale(fm1) <- 2

    expect_equal(attr(VarCorr(fm1), "sc"), 2)
})

test_that("named VarCorr assigned to correct random effects", {

    VarCorr(mod2) <- list(g=1, x=2)

    #expect_equal(c(VarCorr(mod2)$x), 2)
})

