
#' Standard Message
#'
#' @param text String to print out as message.
#'
#' @export

print_message <- function(text) {
  message_text <- str_c(
    "###", str_c("[", Sys.time(), "]:"),
    "...", text, sep = " "
  )

  message(message_text)
}
