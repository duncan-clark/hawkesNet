# Standalone get_network for OpenAlex citation data (for use in openalex_hawkes_study.R).
# Requires: httr, jsonlite, dplyr, purrr, tidyr, gender, genderdata, network.

library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)
library(network)

# Load gender packages if available
gender_available <- FALSE
genderdata_available <- FALSE
if (requireNamespace("gender", quietly = TRUE)) {
  library(gender)
  gender_available <- TRUE
  if (requireNamespace("genderdata", quietly = TRUE)) {
    library(genderdata)
    genderdata_available <- TRUE
  } else {
    message("⚠ Package 'genderdata' not installed. Gender prediction will fail.")
    message("  Install with: install.packages('genderdata')")
    message("  NOTE: This may require internet access to download historical name data.")
  }
} else {
  message("⚠ Package 'gender' not installed; all predicted_gender will be 'unknown'.")
  message("  Install with: install.packages(c('gender', 'genderdata'))")
  message("  NOTE: 'genderdata' may require internet access to download data.")
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
                        max_date = "2100-01-01",
                        min_topic = NULL,
                        remove_topics = NULL,
                        topics_include = NULL
                        ) {
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
  
  if(!is.null(min_topic)){
    topic_counts <- table(nodes$topic)
    nodes <- nodes %>%
      dplyr::mutate(topic = ifelse(topic_counts[topic] < 5, "Other", topic)) %>%
      filter(!topic == "Other")
  }
  
  if(!is.null(remove_topics)){
    remove <- sapply(unique(nodes$topic), function(t) any(grepl(remove_topics, t, ignore.case = TRUE)))
    nodes <- nodes %>%
      filter(!topic %in% remove)
  }
  
  if(!is.null(topics_include)){
    nodes <- nodes %>%
      filter(topic %in% topics_include)
  }
  
  message("Predicting author genders...")
  unique_names <- unique(na.omit(nodes$first_name))
  unique_names <- unique_names[unique_names != ""]
  
  if (length(unique_names) == 0) {
    message("  No first names found; setting all genders to 'unknown'.")
    nodes$predicted_gender <- "unknown"
  } else if (!gender_available || !genderdata_available) {
    warning("⚠ Gender prediction skipped: packages not available.")
    if (!gender_available) {
      warning("  - 'gender' package not installed. Install with: install.packages('gender')")
    }
    if (!genderdata_available) {
      warning("  - 'genderdata' package not installed. Install with: install.packages('genderdata')")
      warning("  - NOTE: Installing 'genderdata' requires internet access to download historical name data.")
      warning("  - On clusters without internet: install 'genderdata' before submitting SLURM jobs.")
    }
    nodes$predicted_gender <- "unknown"
  } else {
    # Both packages are available, try to predict genders
    tryCatch({
      # Load gender function explicitly
      gender_func <- get("gender", envir = asNamespace("gender"))
      
      # Test if gender() function works with a simple test name
      test_result <- tryCatch({
        gender_func("Mary", years = c(1932, 2012), method = "ssa")
      }, error = function(e) {
        NULL
      })
      
      if (is.null(test_result) || nrow(test_result) == 0) {
        warning("⚠ Gender prediction test failed - 'genderdata' package may not have data loaded.")
        warning("  This often happens on clusters without internet access.")
        warning("  Solution: Install 'genderdata' package BEFORE submitting SLURM jobs (requires internet).")
        warning("  Or: Copy genderdata package from a machine with internet to cluster.")
        nodes$predicted_gender <- "unknown"
      } else {
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
        
        message("  Predicting genders for ", length(unique_names), " unique names...")
        message("  Using birth year range: ", years_range[1], "-", years_range[2])
        
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
        
        n_predicted <- sum(nodes$predicted_gender != "unknown", na.rm = TRUE)
        if (n_predicted == 0) {
          warning("⚠ All genders are 'unknown' - gender prediction returned no results.")
          warning("  Possible causes:")
          warning("  1. Names not in historical database (try wider year range)")
          warning("  2. 'genderdata' package data not properly loaded")
          warning("  3. Network/internet issue preventing data download")
        } else {
          message("  Successfully predicted ", n_predicted, " genders (", 
                  round(100 * n_predicted / nrow(nodes), 1), "%)")
        }
      }
    }, error = function(e) {
      warning("⚠ Gender prediction failed with error: ", e$message)
      warning("  Common causes:")
      warning("  1. 'genderdata' package not installed or data not available")
      warning("  2. Cluster has no internet access (install packages before submitting jobs)")
      warning("  3. Package version incompatibility")
      warning("  Setting all genders to 'unknown'.")
      nodes$predicted_gender <<- "unknown"
    })
  }
  min_d <- min(nodes$date, na.rm = TRUE)
  max_d <- max(nodes$date, na.rm = TRUE)
  # One row per work (same id can appear in multiple API pages). Required so the edge join
  # is many-to-one and the network has no multiedges; otherwise structural fit, nodeMatch, GOF and saved edges are wrong.
  nodes <- nodes %>%
    dplyr::mutate(time_scaled = as.numeric(date - min_d) / as.numeric(max_d - min_d)) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  message("Building edge list with temporal attributes...")
  edges <- map_df(all_results, function(x) {
    if (is.null(x$referenced_works) || length(x$referenced_works) == 0) return(NULL)
    tibble(citing_id = x$id, cited_id = unlist(x$referenced_works))
  }) %>%
    dplyr::filter(citing_id %in% nodes$id & cited_id %in% nodes$id) %>%
    dplyr::left_join(nodes %>% dplyr::select(id, date_str, time_scaled), by = c("citing_id" = "id"), relationship = "many-to-one") %>%
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
  network::set.vertex.attribute(net, "topic", as.vector(nodes$topic))
  network::set.vertex.attribute(net, "type", nodes$type)
  network::set.vertex.attribute(net, "author_name", nodes$first_author_name)
  network::set.vertex.attribute(net, "gender", nodes$predicted_gender)
  return(list(net = net, edges = edges,nodes = nodes))
}
