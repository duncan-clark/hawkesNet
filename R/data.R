# Documentation for data files in data/
# These files are not loaded as R objects; they are read by package scripts or the user.

#' CollegeMsg messaging network data
#'
#' @docType data
#' @name CollegeMsg
#' @description
#' College messaging network: directed edges (sender, receiver) with Unix timestamps.
#' Stored as a space-separated text file in the package \code{data/} directory.
#' Used by the message-network and related examples to fit or illustrate Hawkes network growth.
#' @format
#' A text file with one line per message: \code{sender_id receiver_id timestamp} (space-separated).
#' Timestamps are Unix time.
#' @source
#' See references in the package vignettes or \code{inst/message_network/} for data origin.
#' @seealso \code{\link{ht09_contact_list}}, \code{inst/message_network/message_network.R}
NULL

#' Hypertext 2009 contact list (temporal edges)
#'
#' @docType data
#' @name ht09_contact_list
#' @description
#' Contact list / proximity network with temporal edge list: time and node pair.
#' Stored as a tab-separated file in the package \code{data/} directory.
#' Used by the hypertext-conference example to fit or illustrate Hawkes network growth.
#' @format
#' A tab-separated text file: \code{time node_i node_j}.
#' Time is in seconds; node IDs identify individuals.
#' @source
#' Hypertext 2009 conference data; see \code{inst/hypertext_conference/} for usage.
#' @seealso \code{\link{CollegeMsg}}, \code{inst/hypertext_conference/hypertext_conference.R}
NULL
