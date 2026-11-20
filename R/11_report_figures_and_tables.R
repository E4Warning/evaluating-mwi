# Files written for the article: figures copied from the report, tables as CSV
#
#   * copy_report_section_figures() and its helpers find the figures that the
#     rendered analysis report shows in given sections and copy them to the
#     figures folder;
#   * write_tables() writes every summary table to a CSV file.


# Count how many times a pattern occurs in one string.
count_regex_matches <- function(text, pattern) {
  # Positions of all matches in the string (-1 if there is none).
  matches <- gregexpr(pattern, text, perl = TRUE)[[1]]
  # No match: the count is zero.
  if (identical(matches[1], -1L)) {
    return(0L)
  }
  # Otherwise the count is the number of positions found.
  length(matches)
}


# Find the HTML file among the files the report target says it produced.
resolve_report_html_path <- function(report_output) {
  # The report target returns a list of file paths (the HTML, its folder of
  # figures, the source file); flatten it to a plain character vector without
  # repeats.
  candidate_paths <- unique(as.character(unlist(report_output, recursive = TRUE, use.names = FALSE)))
  # Keep the paths ending in ".html".
  html_paths <- candidate_paths[grepl("\\.html$", candidate_paths)]

  # Stop if there is none.
  if (!length(html_paths)) {
    stop("Could not find a rendered HTML report in the supplied tar_quarto output.")
  }

  # Return the full path of the first one (stopping if the file does not exist).
  normalizePath(html_paths[1], winslash = "/", mustWork = TRUE)
}


# Cut one section out of the rendered report: the HTML from the line where
# `<section id="...">` opens to the line where that section closes.
extract_html_section <- function(html_lines, section_id) {
  # Line(s) where the section with this identifier opens.
  section_start <- grep(sprintf('<section id="%s"', section_id), html_lines, fixed = TRUE)

  # Stop if the section is not in the report.
  if (!length(section_start)) {
    stop(sprintf("Section '%s' was not found in the rendered report.", section_id))
  }

  # Start at the first such line, with no section open yet and nothing
  # collected.
  section_start <- section_start[1]
  section_depth <- 0L
  section_lines <- character()

  # Walk down the file from there, one line at a time.
  for (line_index in seq.int(section_start, length(html_lines))) {
    # Add the line to what has been collected.
    this_line <- html_lines[[line_index]]
    section_lines <- c(section_lines, this_line)
    # Keep count of how many sections are open: add one for each opening tag on
    # the line and take one off for each closing tag. Sections can sit inside
    # sections, so the count can go above one.
    section_depth <- section_depth +
      count_regex_matches(this_line, "<section\\b") -
      count_regex_matches(this_line, "</section>")

    # When the count returns to zero, the section we started in has closed:
    # stop.
    if (section_depth == 0L) {
      break
    }
  }

  # Return the collected lines as one string.
  paste(section_lines, collapse = "\n")
}


# Paths of the figure files that the given sections of the rendered report show.
extract_report_section_figure_paths <- function(report_output, section_ids = c("main-results", "data")) {
  # Path of the report's HTML file and the folder it is in.
  report_html_path <- resolve_report_html_path(report_output)
  report_dir <- dirname(report_html_path)
  # Read the HTML, one element per line.
  html_lines <- readLines(report_html_path, warn = FALSE, encoding = "UTF-8")

  # Cut out each requested section (one string per section).
  section_html <- vapply(section_ids, function(section_id) {
    extract_html_section(html_lines, section_id)
  }, character(1))

  # In each section, find every image reference that points into the report's
  # figure folder.
  figure_matches <- regmatches(
    section_html,
    gregexpr('src="analysis_report_files/figure-html/[^"]+"', section_html, perl = TRUE)
  )

  # Strip the `src="` and the closing quote from each, leaving the relative
  # path, and drop repeats.
  relative_paths <- unique(unlist(lapply(figure_matches, function(matches) {
    if (!length(matches)) {
      return(character())
    }
    sub('^src="', "", sub('"$', "", matches))
  })))

  # Stop if no figure was found.
  if (!length(relative_paths)) {
    stop("No figure files were found in the requested report sections.")
  }

  # Return the full paths (stopping if a file does not exist).
  normalizePath(file.path(report_dir, relative_paths), winslash = "/", mustWork = TRUE)
}


# Copy the figures shown in the given sections of the report into another
# folder, and return the paths of the copies.
copy_report_section_figures <- function(report_output, destination_dir, section_ids = c("main-results", "data")) {
  # Figure files to copy.
  source_paths <- extract_report_section_figure_paths(
    report_output = report_output,
    section_ids = section_ids
  )

  # Create the destination folder if needed, and stop if it cannot be written
  # to.
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(destination_dir) || file.access(destination_dir, mode = 2) != 0) {
    stop(sprintf("Unable to create or write to destination directory: %s", destination_dir))
  }

  # Each copy keeps its file name.
  destination_paths <- file.path(destination_dir, basename(source_paths))
  # Copy, replacing any older file of the same name.
  copied <- file.copy(source_paths, destination_paths, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)

  # Stop if any copy failed.
  if (!all(copied)) {
    stop("Failed to copy one or more report figures to the destination directory.")
  }

  # Return the full paths of the copies.
  normalizePath(destination_paths, winslash = "/", mustWork = TRUE)
}


#' Write summary tables to CSV
#'
#' List columns (pointwise values and the like) are dropped.
#'
#' @param tables Named list of data frames.
#' @param dir Output directory.
#' @return The paths written, for a `format = "file"` target.
write_tables <- function(tables, dir) {
  # Create the folder if it is not there.
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # Go through the tables by name.
  vapply(names(tables), function(nm) {
    # Make it a plain data frame.
    x <- as.data.frame(tables[[nm]])
    # Drop columns that hold lists (they cannot be written to a CSV file).
    x <- x[, !vapply(x, is.list, logical(1)), drop = FALSE]
    # The file is named after the table.
    path <- file.path(dir, paste0(nm, ".csv"))
    # Write it and return its path.
    readr::write_csv(x, path)
    path
  }, character(1), USE.NAMES = FALSE)
}
