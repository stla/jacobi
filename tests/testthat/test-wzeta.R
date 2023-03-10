test_that("zeta values", {
  omega <- halfPeriods(c(5+3i, 5+3i))
  expect_equal(
    wzeta(1+1i, omega = omega),
    0.802084165492408 - 0.381791358666872i
  )
  omega <- halfPeriods(c(1, 0))
  expect_equal(
    wzeta(1+1i/2, omega = omega),
    0.796091125108 - 0.422887864713i
  )
  omega <- halfPeriods(c(0, 1))
  expect_equal(
    wzeta(1+1i/2, omega = omega),
    0.80847063824 - 0.409123683392i
  )
  # expect_equal(
  #   wzeta(1+1i, g = c(5+3i, 5+3i)),
  #   0.802084165492408 - 0.381791358666872i
  # )
  # expect_equal(
  #   wzeta(1+1i/2, g = c(1, 0)),
  #   0.796091125108 - 0.422887864713i
  # )
  # expect_equal(
  #   wzeta(1+1i/2, g = c(0, 1)),
  #   0.80847063824 - 0.409123683392i
  # )
})
