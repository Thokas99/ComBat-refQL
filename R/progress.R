new_reporter <- function(verbosity) {
  quiet <- verbosity == "quiet"
  completed <- 0L
  total <- 0L

  text <- function(message) if (!quiet) cli::cli_text(message)
  success <- function(message) if (!quiet) cli::cli_alert_success(message)

  list(
    success = success,
    info = function(message) if (!quiet) cli::cli_alert_info(message),
    warning = function(message) if (!quiet) cli::cli_alert_warning(message),
    reference = function(candidate) {
      if (verbosity == "verbose" && !interactive())
        text(sprintf("Scored reference candidate: %s", candidate))
    },
    start_mapping = function(steps) {
      total <<- steps
    },
    mapping = function(source, genes) {
      completed <<- completed + 1L
      if (verbosity == "verbose" &&
        (completed == 1L || completed == total ||
          completed %% max(1L, ceiling(total / 5)) == 0L)) {
        text(sprintf("Mapped batch %s: %s/%s gene chunks", source, completed, total))
      }
    },
    outcome = function(adjusted, unchanged, failed) {
      message <- sprintf("Adjusted %s | Unchanged %s | Failed %s",
        adjusted, unchanged, failed)
      if (!quiet) {
        if (failed) cli::cli_alert_danger(message) else success(message)
      }
    },
    finish = function(seconds) {
      if (!quiet) {
        cli::cli_text("")
        cli::cli_alert_info("Completed in {format(round(seconds, 1), nsmall = 1)} seconds")
      }
    }
  )
}
