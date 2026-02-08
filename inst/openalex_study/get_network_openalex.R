# Standalone get_network for OpenAlex citation data (for use in openalex_hawkes_study.R).
# Requires: httr, jsonlite, dplyr, purrr, tidyr, gender, genderdata, network.

library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)
library(network)

if (requireNamespace("gender", quietly = TRUE)) {
  library(gender)
  if (requireNamespace("genderdata", quietly = TRUE)) library(genderdata)
} else {
  message("Package 'gender' not installed; all predicted_gender will be 'unknown'.")
}

#' Fetch citation network from OpenAlex API.
#'
#' @param email Email for polite pool (optional).
#' @param pages Number of pages to fetch.
#' @param per_page Results per page.
#' @param string Search string for title_and_abstract.
#' @param min_date Minimum date (character YYYY-MM-DD).
#' @param max_date Maximum date (character YYYY-MM-DD).
#' @return List with \code{net} (network object with time_scaled, gender, etc.) and \code{edges}.
get_network <- function(email = "",
                        pages = 10,
                        per_page = 50,
                        string = "Hawkes",
                        min_date = "1971-04-01",
                        max_date = "2100-01-01") {
  base_url <- "https://api.openalex.org/works"
  params <- list(
    filter = paste0("title_and_abstract.search:", string),
    select = "id,display_name,referenced_works,publication_year,publication_date,authorships,cited_by_count,primary_topic,type,doi",
    per_page = per_page,
    mailto = email,
    cursor = "*"
  )
  all_results <- list()
  message(paste("Fetching", pages, "page(s) from OpenAlex..."))
  for (i in 1:pages) {
    res <- GET(base_url, query = params)
    data <- content(res)
    if (length(data$results) == 0) break
    all_results <- c(all_results, data$results)
    params$cursor <- data$meta$next_cursor
    message(paste("Page", i, "retrieved... Total works:", length(all_results)))
    if (pages > 1) Sys.sleep(0.1)
  }
  nodes <- map_df(all_results, function(x) {
    full_name <- NA
    first_name <- NA
    if (!is.null(x$authorships) && length(x$authorships) > 0) {
      first_author_obj <- x$authorships[[1]]
      if (!is.null(first_author_obj$author$display_name)) {
        full_name <- first_author_obj$author$display_name
        first_name <- sub(" .*", "", full_name)
      }
    }
    tibble(
      id = x$id,
      title = x$display_name,
      date_str = x$publication_date,
      year = x$publication_year,
      citations = x$cited_by_count,
      topic = ifelse(is.null(x$primary_topic$display_name), "Unknown", x$primary_topic$display_name),
      type = x$type,
      doi = ifelse(is.null(x$doi), NA, x$doi),
      first_author_name = full_name,
      first_name = first_name
    )
  }) %>%
    dplyr::mutate(date = as.Date(date_str)) %>%
    dplyr::filter(!is.na(date)) %>%
    dplyr::filter(date >= as.Date(min_date)) %>%
    dplyr::filter(date <= as.Date(max_date)) %>%
    dplyr::arrange(date)
  message("Predicting author genders...")
  unique_names <- unique(na.omit(nodes$first_name))
  unique_names <- unique_names[unique_names != ""]
  if (length(unique_names) > 0 && requireNamespace("gender", quietly = TRUE) && 
      requireNamespace("genderdata", quietly = TRUE)) {
    tryCatch({
      # Load gender function explicitly
      gender_func <- get("gender", envir = asNamespace("gender"))
      
      # Use publication years as proxy for birth years (assuming authors are ~30-50 years old)
      # Calculate approximate birth years from publication years
      pub_years <- unique(na.omit(nodes$year))
      if (length(pub_years) > 0) {
        # Estimate birth years: assume authors are 30-50 years old when publishing
        # Use a wide range to cover most cases
        birth_year_min <- max(1930, min(pub_years, na.rm = TRUE) - 50)
        birth_year_max <- max(pub_years, na.rm = TRUE) - 25
        years_range <- c(birth_year_min, birth_year_max)
      } else {
        # Default to wide range if no years available
        years_range <- c(1932, 2012)
      }
      
      # Call gender() with years parameter
      gender_preds <- gender_func(unique_names, years = years_range, method = "ssa") %>%
        dplyr::select(first_name = name, predicted_gender = gender) %>%
        # Handle case where multiple rows per name (shouldn't happen with unique names, but be safe)
        dplyr::distinct(first_name, .keep_all = TRUE)
      
      nodes <- nodes %>%
        dplyr::left_join(gender_preds, by = "first_name") %>%
        dplyr::mutate(predicted_gender = ifelse(is.na(predicted_gender), "unknown", predicted_gender))
      
      # Check if gender prediction actually worked
      gender_counts <- table(nodes$predicted_gender, useNA = "ifany")
      message("  Gender distribution: ", paste(names(gender_counts), "=", gender_counts, collapse = ", "))
      if (all(nodes$predicted_gender == "unknown", na.rm = TRUE)) {
        warning("⚠ All genders are 'unknown' - gender prediction may have failed. Check if 'gender' and 'genderdata' packages are installed.")
      }
    }, error = function(e) {
      warning("⚠ Gender prediction failed: ", e$message, "\n  Setting all genders to 'unknown'. Install 'gender' and 'genderdata' packages for gender prediction.")
      nodes$predicted_gender <<- "unknown"
    })
  } else {
    warning("⚠ 'gender' or 'genderdata' packages not available. All genders set to 'unknown'.")
    nodes$predicted_gender <- "unknown"
  }
  min_d <- min(nodes$date, na.rm = TRUE)
  max_d <- max(nodes$date, na.rm = TRUE)
  nodes <- nodes %>%
    dplyr::mutate(time_scaled = as.numeric(date - min_d) / as.numeric(max_d - min_d))
  message("Building edge list with temporal attributes...")
  edges <- map_df(all_results, function(x) {
    if (is.null(x$referenced_works) || length(x$referenced_works) == 0) return(NULL)
    tibble(citing_id = x$id, cited_id = unlist(x$referenced_works))
  }) %>%
    dplyr::filter(citing_id %in% nodes$id & cited_id %in% nodes$id) %>%
    dplyr::left_join(nodes %>% dplyr::select(id, date_str, time_scaled), by = c("citing_id" = "id")) %>%
    dplyr::rename(edge_time = date_str, edge_time_scaled = time_scaled)
  nodes$index <- 1:nrow(nodes)
  id_map <- nodes$index
  names(id_map) <- nodes$id
  net <- network::network.initialize(nrow(nodes), directed = FALSE)
  if (nrow(edges) > 0) {
    network::add.edges(net, tail = id_map[edges$citing_id], head = id_map[edges$cited_id])
    network::set.edge.attribute(net, "time", edges$edge_time)
    network::set.edge.attribute(net, "time_scaled", edges$edge_time_scaled)
  }
  network::set.vertex.attribute(net, "title", nodes$title)
  network::set.vertex.attribute(net, "entry_time", nodes$date_str)
  network::set.vertex.attribute(net, "time_scaled", nodes$time_scaled)
  network::set.vertex.attribute(net, "citations", nodes$citations)
  network::set.vertex.attribute(net, "topic", nodes$topic)
  network::set.vertex.attribute(net, "type", nodes$type)
  network::set.vertex.attribute(net, "author_name", nodes$first_author_name)
  network::set.vertex.attribute(net, "gender", nodes$predicted_gender)
  return(list(net = net, edges = edges))
}
