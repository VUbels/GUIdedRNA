# global.R

# Global message queue
message("Loading global.R file...")

# Global message queue - using a hidden environment variable for storage
.message_env <- new.env()
.message_env$queue <- character(0)

# Global message handling functions
send_message <- function(msg) {
  # Add to queue
  .message_env$queue <- c(.message_env$queue, msg)
  # Also log to console
  message(paste("Message queued:", msg))
}

get_messages <- function() {
  msgs <- .message_env$queue
  .message_env$queue <- character(0)
  return(msgs)
}