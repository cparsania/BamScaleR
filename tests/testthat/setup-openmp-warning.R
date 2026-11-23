# On a build without OpenMP (possible on macOS), the first call with
# threads > 1 warns once per session. Mark that warning as already given, so
# it does not show up in unrelated tests (the test of the warning resets the
# flag itself), and restore the session's state when the tests finish.
local({
  state <- BamScaleR:::.bamscale_state
  saved <- get0("warned_no_openmp", envir = state)
  if (identical(BamScaleR:::openmp_max_threads_cpp(), 0L)) {
    assign("warned_no_openmp", TRUE, envir = state)
  }
  withr::defer(
    assign("warned_no_openmp", saved, envir = state),
    testthat::teardown_env()
  )
})
