new_reporter <- function(verbose) {
  success <- function(message) if (verbose) cli::cli_alert_success(message)

  list(
    start = function() if (verbose) cli::cli_h1("ComBat-refQL"),
    success = success,
    info = function(message) if (verbose) cli::cli_alert_info(message),
    warning = function(message) cli::cli_warn(c("!" = message)),
    outcome = function(adjusted, unchanged, failed) {
      message <- sprintf("Adjusted %s | Unchanged %s | Failed %s",
        adjusted, unchanged, failed)
      if (verbose) {
        if (failed) cli::cli_alert_danger(message) else success(message)
      }
    },
    finish = function(seconds) {
      if (verbose) {
        cli::cli_text("")
        cli::cli_alert_info("Runtime: {format(round(seconds, 1), nsmall = 1)} s")
      }
    }
  )
}
