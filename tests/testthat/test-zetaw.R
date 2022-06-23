test_that("zeta values", {
  expect_equal(
    zetaw(1+1i, g = c(5+3i, 5+3i)),
    0.802084165492408 - 0.381791358666872i
  )
})
