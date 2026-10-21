# 2026-10-19
# Visualize CARP consensus annotations from Censor and WU-BLAST
# This script does not read RepeatMasker genome-coordinate output
library(Gviz)
library(grid)

############################################
# PART 1: Define shared annotation settings
############################################

# Keep annotation sources in the same order in every PNG
DATABASE_ORDER <- c(
   "RepBase",
   "our-known",
   "Swiss-Prot",
   "GB_TE",
   "Retrovirus"
)

# Define match classes using the Java classification criteria
MATCH_LABELS <- c(
   "Known match",
   "Partial match"
)

# Use one fixed color for each match class
MATCH_COLORS <- c(
   "Known match" = "red",
   "Partial match" = "orange"
)

# Return the shared column structure used by every annotation reader
empty_annotation_data <- function() {
   data.frame(
      Consensus = character(),
      query_from = numeric(),
      query_to = numeric(),
      reference = character(),
      reference_length = numeric(),
      reference_from = numeric(),
      reference_to = numeric(),
      direction = character(),
      identity = numeric(),
      score = numeric(),
      score_type = character(),
      evalue = numeric(),
      database = character(),
      stringsAsFactors = FALSE
   )
}


############################################
# PART 2: Read FASTA information
############################################

# file_path:
#   Path to a nucleotide FASTA file
#
# file_label:
#   Human-readable name used in error messages
#   Example: file_label = "RepBase"
read_fasta_lengths <- function(
   file_path,
   file_label) {
   if (!file.exists(file_path)) {
      stop(
         "[ERROR] ",
         file_label,
         " FASTA file does not exist: ",
         file_path,
         call. = FALSE
      )
   }

   fasta_lines <- readLines(
      file_path,
      warn = FALSE
   )

   header_lines <- fasta_lines[
      startsWith(fasta_lines, ">")
   ]

   if (length(header_lines) == 0L) {
      stop(
         "[ERROR] ",
         file_label,
         " FASTA contains no headers",
         call. = FALSE
      )
   }

   # Remove the leading greater-than symbol and header description
   # Example: ">BovB_BosInd RTE" becomes "BovB_BosInd"
   fasta_ids <- sub(
      "[[:space:]].*$",
      "",
      substring(header_lines, 2L)
   )

   if (any(!nzchar(fasta_ids))) {
      stop(
         "[ERROR] ",
         file_label,
         " FASTA contains an empty ID",
         call. = FALSE
      )
   }

   duplicate_ids <- unique(
      fasta_ids[duplicated(fasta_ids)]
   )

   if (length(duplicate_ids) > 0L) {
      stop(
         "[ERROR] Duplicate ID in ",
         file_label,
         " FASTA: ",
         paste(duplicate_ids, collapse = ", "),
         call. = FALSE
      )
   }

   # Read the length of every sequence from the FASTA file
   # seqtype = "B" counts all sequence characters without interpreting
   # them as DNA, RNA, or protein characters
   # This prevents non-standard RepBase symbols from being discarded
   sequence_lengths <- tryCatch(
      Biostrings::fasta.seqlengths(
         filepath = file_path,
         seqtype = "B",
         use.names = TRUE
      ),
      error = function(error) {
         stop(
            "[ERROR] Could not read ",
            file_label,
            " FASTA: ",
            conditionMessage(error),
            call. = FALSE
         )
      }
   )

   # The number of calculated lengths must match the number of FASTA IDs
   if (length(sequence_lengths) != length(fasta_ids)) {
      stop(
         "[ERROR] Could not match all sequence lengths in ",
         file_label,
         " FASTA",
         call. = FALSE
      )
   }

   if (any(sequence_lengths < 1L)) {
      stop(
         "[ERROR] ",
         file_label,
         " FASTA contains a zero-length sequence",
         call. = FALSE
      )
   }

   # Use FASTA IDs as names so a length can be retrieved by reference ID
   # Example: sequence_lengths[["BovB_BosInd"]] returns 3769
   names(sequence_lengths) <- fasta_ids

   return(sequence_lengths)
}


# file_path:
#   Path to the consensus FASTA associated with the annotation results
#   Sequence IDs are taken from the first word of each FASTA header
read_consensus_lengths <- function(file_path) {
   if (!file.exists(file_path)) {
      stop(
         "[ERROR] Consensus FASTA file does not exist: ",
         file_path,
         call. = FALSE
      )
   }

   sequence_lengths <- tryCatch(
      Biostrings::fasta.seqlengths(
         filepath = file_path,
         seqtype = "DNA",
         use.names = TRUE
      ),
      error = function(error) {
         stop(
            "[ERROR] Could not read consensus FASTA: ",
            conditionMessage(error),
            call. = FALSE
         )
      }
   )

   if (length(sequence_lengths) == 0L) {
      stop(
         "[ERROR] Consensus FASTA contains no sequences",
         call. = FALSE
      )
   }

   # Keep only the first word of each complete FASTA header
   consensus_ids <- sub(
      "[[:space:]].*$",
      "",
      names(sequence_lengths)
   )

   if (any(!nzchar(consensus_ids))) {
      stop(
         "[ERROR] Consensus FASTA contains an empty ID",
         call. = FALSE
      )
   }

   duplicate_ids <- unique(
      consensus_ids[duplicated(consensus_ids)]
   )

   if (length(duplicate_ids) > 0L) {
      stop(
         "[ERROR] Duplicate consensus FASTA ID: ",
         paste(duplicate_ids, collapse = ", "),
         call. = FALSE
      )
   }

   if (any(sequence_lengths < 1L)) {
      stop(
         "[ERROR] Consensus FASTA contains a zero-length sequence",
         call. = FALSE
      )
   }

   data.frame(
      Consensus = consensus_ids,
      length = unname(sequence_lengths),
      stringsAsFactors = FALSE
   )
}


############################################
# PART 3: Read the 12-column Censor map
############################################

# file_path:
#   Path to the Censor repeat-library map
#
# repbase_file:
#   RepBase FASTA used by Censor
#
# our_known_file:
#   our-known FASTA used by Censor
#
# The map reference ID is matched exactly against both FASTA files
# An ID present in both or neither library is rejected
read_map_file <- function(
    file_path,
    repbase_file,
    our_known_file) {
   if (!file.exists(file_path)) {
      stop(
         "[ERROR] Censor map file does not exist: ",
         file_path,
         call. = FALSE
      )
   }

   map_lines <- readLines(
      file_path,
      warn = FALSE
   )

   map_lines <- map_lines[
      nzchar(trimws(map_lines))
   ]

   # An empty map is valid when Censor finds no repeat-library hits
   if (length(map_lines) == 0L) {
      return(empty_annotation_data())
   }

   fields <- strsplit(
      trimws(map_lines),
      "[[:space:]]+"
   )

   bad_rows <- which(
      lengths(fields) != 12L
   )

   if (length(bad_rows) > 0L) {
      stop(
         "[ERROR] Expected 12 Censor map columns at non-empty row: ",
         bad_rows[1],
         call. = FALSE
      )
   }
   map_data <- read.table(
      text = paste(map_lines, collapse = "\n"),
      header = FALSE,
      colClasses = "character",
      comment.char = "",
      quote = "",
      stringsAsFactors = FALSE
   )

   # Columns 11 and 12 are retained as Censor-specific metrics
   # They are not interpreted as BLAST bit scores
   names(map_data) <- c(
      "Consensus",
      "query_from",
      "query_to",
      "reference",
      "reference_from",
      "reference_to",
      "direction",
      "identity",
      "metric_9",
      "censor_score",
      "metric_11",
      "metric_12"
   )

   numeric_columns <- c(
      "query_from",
      "query_to",
      "reference_from",
      "reference_to",
      "identity",
      "metric_9",
      "censor_score",
      "metric_11",
      "metric_12"
   )

   for (column in numeric_columns) {
      map_data[[column]] <- suppressWarnings(
         as.numeric(map_data[[column]])
      )

      if (any(!is.finite(map_data[[column]]))) {
         stop(
            "[ERROR] Invalid numeric value in Censor map column: ",
            column,
            call. = FALSE
         )
      }
   }

   coordinate_columns <- c(
      "query_from",
      "query_to",
      "reference_from",
      "reference_to"
   )

   for (column in coordinate_columns) {
      values <- map_data[[column]]

      if (any(values < 1 | values != floor(values))) {
         stop(
            "[ERROR] Censor map coordinates must be positive integers: ",
            column,
            call. = FALSE
         )
      }
   }

   if (any(map_data$query_from > map_data$query_to)) {
      stop(
         "[ERROR] A Censor map query interval starts after it ends",
         call. = FALSE
      )
   }

   if (any(map_data$reference_from > map_data$reference_to)) {
      stop(
         "[ERROR] A Censor map reference interval starts after it ends",
         call. = FALSE
      )
   }

   if (any(map_data$identity < 0 | map_data$identity > 1)) {
      stop(
         "[ERROR] Censor map identity must be in [0, 1]",
         call. = FALSE
      )
   }

   if (any(map_data$censor_score < 0)) {
      stop(
         "[ERROR] Censor score must be non-negative",
         call. = FALSE
      )
   }

   if (any(!map_data$direction %in% c("d", "c"))) {
      stop(
         "[ERROR] Censor map direction must be d or c",
         call. = FALSE
      )
   }

   # Named vectors contain both the reference IDs and their full lengths
   # Example: repbase_lengths[["BovB_BosInd"]] returns 3769
   repbase_lengths <- read_fasta_lengths(
      repbase_file,
      "RepBase"
   )

   our_known_lengths <- read_fasta_lengths(
      our_known_file,
      "our-known"
   )

   # names() extracts the sequence IDs from each named length vector
   repbase_ids <- names(repbase_lengths)
   our_known_ids <- names(our_known_lengths)

   in_repbase <- map_data$reference %in% repbase_ids
   in_our_known <- map_data$reference %in% our_known_ids

   ambiguous_rows <- which(
      in_repbase & in_our_known
   )

   if (length(ambiguous_rows) > 0L) {
      stop(
         "[ERROR] Map reference occurs in both RepBase and our-known: ",
         map_data$reference[ambiguous_rows[1]],
         call. = FALSE
      )
   }

   unknown_rows <- which(
      !in_repbase & !in_our_known
   )

   if (length(unknown_rows) > 0L) {
      stop(
         "[ERROR] Map reference is absent from RepBase and our-known: ",
         map_data$reference[unknown_rows[1]],
         call. = FALSE
      )
   }

   database <- ifelse(
      in_repbase,
      "RepBase",
      "our-known"
   )
   # Match each reference ID to its full database-sequence length
   # Example: BovB_BosInd is 3769 bp long in RepBase
   reference_length <- ifelse(
      in_repbase,
      repbase_lengths[map_data$reference],
      our_known_lengths[map_data$reference]
   )

   # Convert Censor direction symbols into plotting strand symbols
   direction <- ifelse(
      map_data$direction == "d",
      "+",
      "-"
   )

   # Display database coordinates opposite to the query arrow direction
   # A right-pointing arrow shows the database end on the left
   # A left-pointing arrow shows the database start on the left
   reference_left <- ifelse(
      map_data$direction == "d",
      map_data$reference_to,
      map_data$reference_from
   )

   reference_right <- ifelse(
      map_data$direction == "d",
      map_data$reference_from,
      map_data$reference_to
   )

   data.frame(
      Consensus = map_data$Consensus,
      query_from = map_data$query_from,
      query_to = map_data$query_to,
      reference = map_data$reference,

      # unname() removes FASTA IDs stored as vector names
      reference_length = unname(reference_length),

      reference_from = reference_left,
      reference_to = reference_right,
      direction = direction,
      identity = map_data$identity,
      score = map_data$censor_score,
      score_type = rep(
         "Censor score",
         nrow(map_data)
      ),
      evalue = rep(
         NA_real_,
         nrow(map_data)
      ),
      database = database,
      stringsAsFactors = FALSE
   )
}


############################################
# PART 4: Read one WU-BLAST text report
############################################

# Convert values such as e-120 into valid R scientific notation
# Example: e-120 becomes 1e-120
parse_wublast_evalue <- function(value) {
   value <- trimws(value)

   if (grepl("^[eE][+-]?[0-9]+$", value)) {
      value <- paste0("1", value)
   }

   numeric_value <- suppressWarnings(
      as.numeric(value)
   )

   if (length(numeric_value) != 1L || is.na(numeric_value)) {
      stop(
         "[ERROR] Invalid WU-BLAST E-value: ",
         value,
         call. = FALSE
      )
   }

   return(numeric_value)
}


# file_path:
#   Path to one WU-BLAST text report such as a .spwb file
#
# database:
#   Source name displayed beside every parsed hit
#   Expected value is Swiss-Prot, GB_TE, or Retrovirus
read_wublast_hits <- function(
    file_path,
    database) {
   if (!file.exists(file_path)) {
      stop(
         "[ERROR] WU-BLAST report does not exist: ",
         file_path,
         call. = FALSE
      )
   }

   if (!database %in% DATABASE_ORDER) {
      stop(
         "[ERROR] Unsupported database name: ",
         database,
         call. = FALSE
      )
   }

   report_lines <- readLines(
      file_path,
      warn = FALSE
   )

   parsed_hits <- list()
   current_query <- NULL
   current_reference <- NULL

   # Store the Length value belonging to the current database sequence
   # Example: "Length = 1281" becomes numeric value 1281
   current_reference_length <- NULL

   current_hit <- NULL

   # Store the current HSP before a new HSP, subject, or query begins
   finish_hit <- function() {
      if (is.null(current_hit)) {
         return(invisible(NULL))
      }

      required_values <- c(
         current_hit$query_start,
         current_hit$query_end,
         current_hit$subject_start,
         current_hit$subject_end,
         current_hit$subject_length,
         current_hit$identity,
         current_hit$bit_score,
         current_hit$evalue
      )

      if (any(is.na(required_values))) {
         stop(
            "[ERROR] Incomplete WU-BLAST HSP for query ",
            current_hit$query_id,
            " and reference ",
            current_hit$subject_id,
            call. = FALSE
         )
      }

      # Record whether the original query coordinates increase left to right
      query_is_forward <-
         current_hit$query_start <= current_hit$query_end

      # Query coordinates determine the arrow direction on the consensus
      query_direction <- ifelse(
         query_is_forward,
         "+",
         "-"
      )

      # Display subject coordinates opposite to the query arrow direction
      # A forward query shows the subject end on the left
      # A reverse query shows the subject start on the left
      reference_left <- ifelse(
         query_is_forward,
         current_hit$subject_end,
         current_hit$subject_start
      )

      reference_right <- ifelse(
         query_is_forward,
         current_hit$subject_start,
         current_hit$subject_end
      )

      parsed_hits[[length(parsed_hits) + 1L]] <<- data.frame(
         Consensus = current_hit$query_id,
         query_from = min(
            current_hit$query_start,
            current_hit$query_end
         ),
         query_to = max(
            current_hit$query_start,
            current_hit$query_end
         ),
         reference = current_hit$subject_id,

         # Store the complete subject length, not the alignment length
         reference_length = current_hit$subject_length,

         reference_from = reference_left,
         reference_to = reference_right,
         direction = query_direction,
         identity = current_hit$identity,
         score = current_hit$bit_score,
         score_type = "Bit score",
         evalue = current_hit$evalue,
         database = database,
         stringsAsFactors = FALSE
      )

      current_hit <<- NULL
      invisible(NULL)
   }

   for (line_number in seq_along(report_lines)) {
      line <- report_lines[[line_number]]

      # Begin a new query record
      if (grepl("^Query=", line)) {
         finish_hit()

         query_match <- regexec(
            "^Query=[[:space:]]*([^[:space:]]+)",
            line
         )

         query_parts <- regmatches(
            line,
            query_match
         )[[1]]

         if (length(query_parts) != 2L) {
            stop(
               "[ERROR] Could not parse WU-BLAST query ID at line ",
               line_number,
               call. = FALSE
            )
         }

         current_query <- query_parts[[2]]
         current_reference <- NULL
         current_reference_length <- NULL
         next
      }

      # Begin a new reference record
      if (startsWith(line, ">")) {
         finish_hit()

         subject_match <- regexec(
            "^>([^[:space:]]+)",
            line
         )

         subject_parts <- regmatches(
            line,
            subject_match
         )[[1]]

         if (length(subject_parts) != 2L) {
            stop(
               "[ERROR] Could not parse WU-BLAST reference ID at line ",
               line_number,
               call. = FALSE
            )
         }

         current_reference <- subject_parts[[2]]

         # The next Length line must belong to this new reference
         current_reference_length <- NULL
         next
      }

      # Read the full length of the current database sequence
      # blastx reports protein lengths, while tblastx reports nucleotide lengths
      if (!is.null(current_reference) &&
         grepl("^[[:space:]]*Length[[:space:]]*=", line)) {
         length_match <- regexec(
            "^[[:space:]]*Length[[:space:]]*=[[:space:]]*([0-9,]+)[[:space:]]*$",
            line
         )

         length_parts <- regmatches(
            line,
            length_match
         )[[1]]

         if (length(length_parts) != 2L) {
            stop(
               "[ERROR] Could not parse WU-BLAST reference length at line ",
               line_number,
               call. = FALSE
            )
         }

         # Remove thousands separators before numeric conversion
         # Example: "7,430" becomes 7430
         current_reference_length <- as.numeric(
            gsub(",", "", length_parts[[2]])
         )

         if (!is.finite(current_reference_length) ||
            current_reference_length < 1 ||
            current_reference_length != floor(current_reference_length)) {
            stop(
               "[ERROR] Invalid WU-BLAST reference length at line ",
               line_number,
               call. = FALSE
            )
         }

         next
      }
      # Begin a new high-scoring segment pair
      if (grepl("^[[:space:]]*Score[[:space:]]*=", line)) {
         finish_hit()

         # A valid HSP requires the query ID, reference ID, and reference length
         if (is.null(current_query) ||
            is.null(current_reference) ||
            is.null(current_reference_length)) {
            stop(
               "[ERROR] WU-BLAST score appears before a query, reference, or reference length at line ",
               line_number,
               call. = FALSE
            )
         }

         bit_score_match <- regexec(
            "\\(([0-9.eE+-]+)[[:space:]]+bits\\)",
            line,
            perl = TRUE
         )

         bit_score_parts <- regmatches(
            line,
            bit_score_match
         )[[1]]

         evalue_match <- regexec(
            "Expect(?:\\([0-9]+\\))?[[:space:]]*=[[:space:]]*([^,[:space:]]+)",
            line,
            perl = TRUE
         )

         evalue_parts <- regmatches(
            line,
            evalue_match
         )[[1]]

         if (length(bit_score_parts) != 2L ||
            length(evalue_parts) != 2L) {
            stop(
               "[ERROR] Could not parse WU-BLAST score at line ",
               line_number,
               call. = FALSE
            )
         }

         bit_score <- suppressWarnings(
            as.numeric(bit_score_parts[[2]])
         )

         if (!is.finite(bit_score) || bit_score < 0) {
            stop(
               "[ERROR] Invalid WU-BLAST bit score at line ",
               line_number,
               call. = FALSE
            )
         }

         current_hit <- list(
            query_id = current_query,
            subject_id = current_reference,
            query_start = NA_real_,
            query_end = NA_real_,
            subject_start = NA_real_,
            subject_end = NA_real_,

            # Copy the complete reference length into this HSP
            subject_length = current_reference_length,

            identity = NA_real_,
            bit_score = bit_score,
            evalue = parse_wublast_evalue(
               evalue_parts[[2]]
            )
         )

         next
      }

      # Read the exact identity fraction instead of the rounded percentage
      if (!is.null(current_hit) &&
         grepl("^[[:space:]]*Identities[[:space:]]*=", line)) {
         identity_match <- regexec(
            "Identities[[:space:]]*=[[:space:]]*([0-9]+)/([0-9]+)",
            line
         )

         identity_parts <- regmatches(
            line,
            identity_match
         )[[1]]

         if (length(identity_parts) != 3L) {
            stop(
               "[ERROR] Could not parse WU-BLAST identity at line ",
               line_number,
               call. = FALSE
            )
         }

         identical_positions <- as.numeric(
            identity_parts[[2]]
         )

         aligned_positions <- as.numeric(
            identity_parts[[3]]
         )

         if (aligned_positions < 1 ||
            identical_positions < 0 ||
            identical_positions > aligned_positions) {
            stop(
               "[ERROR] Invalid WU-BLAST identity at line ",
               line_number,
               call. = FALSE
            )
         }

         current_hit$identity <- identical_positions /
            aligned_positions

         next
      }

      # Record the first and last query coordinates across alignment blocks
      if (!is.null(current_hit) && startsWith(line, "Query:")) {
         query_coordinate_match <- regexec(
            "^Query:[[:space:]]*([0-9,]+)[[:space:]]+[^[:space:]]+[[:space:]]+([0-9,]+)[[:space:]]*$",
            line
         )

         query_coordinate_parts <- regmatches(
            line,
            query_coordinate_match
         )[[1]]

         if (length(query_coordinate_parts) != 3L) {
            stop(
               "[ERROR] Could not parse WU-BLAST query coordinates at line ",
               line_number,
               call. = FALSE
            )
         }

         block_start <- as.numeric(
            gsub(",", "", query_coordinate_parts[[2]])
         )

         block_end <- as.numeric(
            gsub(",", "", query_coordinate_parts[[3]])
         )

         if (is.na(current_hit$query_start)) {
            current_hit$query_start <- block_start
         }

         current_hit$query_end <- block_end
         next
      }

      # Record the first and last reference coordinates across alignment blocks
      if (!is.null(current_hit) && startsWith(line, "Sbjct:")) {
         subject_coordinate_match <- regexec(
            "^Sbjct:[[:space:]]*([0-9,]+)[[:space:]]+[^[:space:]]+[[:space:]]+([0-9,]+)[[:space:]]*$",
            line
         )

         subject_coordinate_parts <- regmatches(
            line,
            subject_coordinate_match
         )[[1]]

         if (length(subject_coordinate_parts) != 3L) {
            stop(
               "[ERROR] Could not parse WU-BLAST reference coordinates at line ",
               line_number,
               call. = FALSE
            )
         }

         block_start <- as.numeric(
            gsub(",", "", subject_coordinate_parts[[2]])
         )

         block_end <- as.numeric(
            gsub(",", "", subject_coordinate_parts[[3]])
         )

         if (is.na(current_hit$subject_start)) {
            current_hit$subject_start <- block_start
         }

         current_hit$subject_end <- block_end
      }
   }

   # Store the final HSP after the last report line
   finish_hit()

   if (length(parsed_hits) == 0L) {
      return(empty_annotation_data())
   }

   wublast_data <- do.call(
      rbind,
      parsed_hits
   )

   if (any(
      wublast_data$query_from < 1 |
      wublast_data$query_to < 1 |
      wublast_data$reference_from < 1 |
      wublast_data$reference_to < 1
   )) {
      stop(
         "[ERROR] WU-BLAST coordinates must be positive integers",
         call. = FALSE
      )
   }

   return(wublast_data)
}


############################################
# PART 5: Build one Gviz track for one hit
############################################

# hit:
#   One row from the combined annotation data frame
#
# hit_color:
#   Color selected from the shared identity interval
create_hit_track <- function(
   hit,
   hit_color) {
   database <- hit$database[[1]]
   reference <- hit$reference[[1]]

   # Full database-sequence length
   # Example: BovB_BosInd has reference_length = 3769
   reference_length <- hit$reference_length[[1]]


   identity <- hit$identity[[1]]
   score <- hit$score[[1]]
   score_type <- hit$score_type[[1]]
   evalue <- hit$evalue[[1]]

   # blastx databases contain proteins, so their lengths are amino acids
   # Censor and tblastx reference libraries contain nucleotide sequences
   reference_unit <- ifelse(
      database %in% c("Swiss-Prot", "GB_TE"),
      "aa",
      "bp"
   )

   if (is.na(evalue)) {
      statistics_text <- sprintf(
         "%s: %.1f | Identity: %.1f%%",
         score_type,
         score,
         100 * identity
      )
   } else {
      statistics_text <- sprintf(
         "%s: %.1f | Identity: %.1f%% | E-value: %.2g",
         score_type,
         score,
         100 * identity,
         evalue
      )
   }

   # Reference coordinates are drawn separately at both ends of every arrow
   track_name <- sprintf(
      "%s | %s (%d %s)\n%s",
      database,
      reference,
      reference_length,
      reference_unit,
      statistics_text
   )

   # Use dark text on the lighter partial-match color
   text_color <- if (
      hit_color == MATCH_COLORS[["Partial match"]]
   ) {
      "black"
   } else {
      "white"
   }

   AnnotationTrack(
      name = "",
      start = hit$query_from[[1]],
      end = hit$query_to[[1]],
      strand = hit$direction[[1]],
      chromosome = "chrNA",
      id = track_name,
      group = track_name,
      showFeatureId = FALSE,
      showId = TRUE,
      just.group = "below",
      fontsize.group = 12,
      fontcolor.group = "black",
      rotation.group = 0,
      showTitle = FALSE,
      stacking = "squish",
      collapse = FALSE,
      shape = "fixedArrow",
      arrowHeadWidth = 10,
      stackHeight = 0.20,
      fill = hit_color,
      col = hit_color,
      fontsize.item = 9,
      fontcolor.item = text_color,
      rotation = 0
   )
}


############################################
# PART 6: Generate one PNG per consensus
############################################

# data:
#   Combined annotation data returned by the map and WU-BLAST readers
#
# consensus_lengths:
#   Data frame with Consensus and length columns
#
# output_dir:
#   Directory for PNG files and plot_manifest.tsv
#
# width_in:
#   PNG width in inches with a default of 14
#
# dpi:
#   PNG resolution with a default of 300 pixels per inch
generate_plots <- function(
    data,
    consensus_lengths,
    output_dir,
    width_in = 14,
    dpi = 300) {
   required_columns <- names(
      empty_annotation_data()
   )

   missing_columns <- setdiff(
      required_columns,
      names(data)
   )

   if (length(missing_columns) > 0L) {
      stop(
         "[ERROR] Annotation data is missing columns: ",
         paste(missing_columns, collapse = ", "),
         call. = FALSE
      )
   }

   unknown_databases <- setdiff(
      unique(data$database),
      DATABASE_ORDER
   )

   if (length(unknown_databases) > 0L) {
      stop(
         "[ERROR] Unsupported databases: ",
         paste(unknown_databases, collapse = ", "),
         call. = FALSE
      )
   }

   unknown_ids <- setdiff(
      unique(data$Consensus),
      consensus_lengths$Consensus
   )

   if (length(unknown_ids) > 0L) {
      stop(
         "[ERROR] Annotation ID is absent from consensus FASTA: ",
         paste(unknown_ids, collapse = ", "),
         call. = FALSE
      )
   }

   # Match sequence lengths by consensus ID rather than row order
   data$length <- consensus_lengths$length[
      match(
         data$Consensus,
         consensus_lengths$Consensus
      )
   ]

   bad_rows <- which(
      data$query_from < 1 |
      data$query_from > data$query_to |
      data$query_to > data$length
   )

   if (length(bad_rows) > 0L) {
      bad_row <- bad_rows[1]

      stop(
         "[ERROR] Annotation interval is outside the consensus: ",
         data$Consensus[bad_row],
         "; interval=",
         data$query_from[bad_row],
         "-",
         data$query_to[bad_row],
         "; length=",
         data$length[bad_row],
         call. = FALSE
      )
   }
   # Recover the lower and upper database coordinates
   # Display order may be reversed, so pmin and pmax are required
   reference_start <- pmin(
      data$reference_from,
      data$reference_to
   )

   reference_end <- pmax(
      data$reference_from,
      data$reference_to
   )

   # Reproduce Java Math.round for positive percentage values
   # Java Math.round(x) is equivalent to floor(x + 0.5)
   data$bpc <- floor(
      100 *
         (data$query_to - data$query_from + 1) /
         data$length +
         0.5
   )

   data$lspc <- floor(
      100 *
         reference_start /
         data$reference_length +
         0.5
   )

   data$lepc <- floor(
      100 *
         reference_end /
         data$reference_length +
         0.5
   )

   # A match is known only when all three Java thresholds pass
   data$match_class <- ifelse(
      data$bpc >= 85 &
         data$lspc <= 5 &
         data$lepc >= 95,
      "Known match",
      "Partial match"
   )
   # Restrict IDs to characters that are safe in PNG filenames
   if (any(!grepl(
      "^[A-Za-z0-9][A-Za-z0-9_.-]*$",
      consensus_lengths$Consensus
   ))) {
      stop(
         "[ERROR] A consensus ID contains unsafe filename characters",
         call. = FALSE
      )
   }

   dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
   )

   if (!dir.exists(output_dir)) {
      stop(
         "[ERROR] Could not create output directory: ",
         output_dir,
         call. = FALSE
      )
   }

   all_plots <- vector(
      "list",
      nrow(consensus_lengths)
   )

   for (
      consensus_index in seq_len(nrow(consensus_lengths))
   ) {
      consensus <- consensus_lengths$Consensus[
         consensus_index
      ]

      consensus_length <- consensus_lengths$length[
         consensus_index
      ]

      subset_data <- data[
         data$Consensus == consensus,
         ,
         drop = FALSE
      ]

      # Group hits by database and then place higher-identity hits first
      subset_data <- subset_data[
         order(
            match(
               subset_data$database,
               DATABASE_ORDER
            ),
            -subset_data$identity,
            -subset_data$score,
            subset_data$query_from,
            subset_data$query_to,
            subset_data$reference
         ),
         ,
         drop = FALSE
      ]

      n_hits <- nrow(subset_data)

      # Select colors from the Java-compatible match classification
      hit_colors <- unname(
         MATCH_COLORS[
            subset_data$match_class
         ]
      )

      track_list <- vector(
         "list",
         n_hits
      )

      for (hit_index in seq_len(n_hits)) {
         hit <- subset_data[
            hit_index,
            ,
            drop = FALSE
         ]

         track_list[[hit_index]] <- create_hit_track(
            hit = hit,
            hit_color = hit_colors[hit_index]
         )
      }

      # Increase image height as the number of hits increases
      height_in <- max(
         3.5,
         2.8 + 0.5 * n_hits
      )

      title_text <- sprintf(
         "%s | %d bp | %d matches",
         consensus,
         consensus_length,
         n_hits
      )

      png_name <- paste0(
         consensus,
         ".png"
      )

      png_path <- file.path(
         output_dir,
         png_name
      )

      png(
         filename = png_path,
         width = width_in,
         height = height_in,
         units = "in",
         res = dpi,
         bg = "white"
      )

      # Always close the PNG device if drawing raises an error
      tryCatch({
         grid.newpage()

         # Row 1 contains the match-classification legend
         # Row 2 contains the consensus and annotation tracks
         pushViewport(viewport(
            layout = grid.layout(
               nrow = 2,
               ncol = 1,
               heights = unit(
                  c(
                     0.55,
                     height_in - 0.55
                  ),
                  "inches"
               )
            )
         ))

         ####################################
         # Draw the match-classification legend
         ####################################

         pushViewport(viewport(
            layout.pos.row = 1
         ))

         grid.text(
            "Match classification",
            x = 0.02,
            y = 0.82,
            just = "left",
            gp = gpar(
               fontsize = 10
            )
         )

         legend_x <- seq(
            0.15,
            0.55,
            length.out = length(MATCH_LABELS)
         )

         for (
            legend_index in seq_along(MATCH_LABELS)
         ) {
            grid.rect(
               x = legend_x[legend_index],
               y = 0.29,
               width = 0.035,
               height = 0.25,
               gp = gpar(
                  fill = MATCH_COLORS[
                     MATCH_LABELS[legend_index]
                  ],
                  col = NA
               )
            )

            grid.text(
               MATCH_LABELS[legend_index],
               x = legend_x[legend_index] + 0.025,
               y = 0.35,
               just = "left",
               gp = gpar(
                  fontsize = 10
               )
            )
         }

         popViewport()

         ####################################
         # Draw the consensus and hit tracks
         ####################################

         pushViewport(viewport(
            layout.pos.row = 2,
            clip = "off"
         ))

         pushViewport(viewport(
            width = 0.94,
            height = 0.90,
            clip = "off",
            layout = grid.layout(
               nrow = 2,
               ncol = 1,
               heights = unit.c(
                  unit(
                     0.75,
                     "inches"
                  ),
                  unit(
                     1,
                     "null"
                  )
               )
            )
         ))

         # Add coordinate padding around both consensus ends
         coordinate_padding <- max(
            10,
            ceiling(
               consensus_length * 0.03
            )
         )

         view_from <- 1 - coordinate_padding
         view_to <- consensus_length +
            coordinate_padding

         ####################################
         # Draw the consensus bar
         ####################################

         pushViewport(viewport(
            layout.pos.row = 1
         ))

         pushViewport(viewport(
            xscale = c(
               view_from,
               view_to
            ),
            yscale = c(
               0,
               1
            ),
            clip = "off"
         ))

         grid.rect(
            x = unit(
               (consensus_length + 1) / 2,
               "native"
            ),
            y = 0.72,
            width = unit(
               consensus_length,
               "native"
            ),
            height = 0.32,
            gp = gpar(
               fill = "#52C3C5",
               col = NA
            )
         )

         grid.text(
            title_text,
            x = unit(
               (consensus_length + 1) / 2,
               "native"
            ),
            y = 0.72,
            gp = gpar(
               fontsize = 12,
               fontface = "bold"
            )
         )

         # Select readable coordinate tick positions
         tick_positions <- pretty(
            c(
               1,
               consensus_length
            ),
            n = 5
         )

         # Remove interior ticks too close to either endpoint
         tick_positions <- tick_positions[
            tick_positions >
               1 + 0.04 * consensus_length &
            tick_positions <
               consensus_length -
               0.04 * consensus_length
         ]

         # Always include the first and final consensus positions
         tick_positions <- sort(unique(c(
            1,
            tick_positions,
            consensus_length
         )))

         grid.segments(
            x0 = unit(
               tick_positions,
               "native"
            ),
            x1 = unit(
               tick_positions,
               "native"
            ),
            y0 = 0.56,
            y1 = 0.36,
            gp = gpar(
               col = "black",
               lwd = 1.5
            )
         )

         grid.text(
            format(
               tick_positions,
               trim = TRUE,
               scientific = FALSE
            ),
            x = unit(
               tick_positions,
               "native"
            ),
            y = 0.29,
            just = c(
               "center",
               "top"
            ),
            gp = gpar(
               fontsize = 12
            )
         )

         popViewport(2)

         ####################################
         # Draw every annotation track
         ####################################

         pushViewport(viewport(
            layout.pos.row = 2,
            clip = "off"
         ))

         if (n_hits > 0L) {
            plotTracks(
               track_list,
               sizes = rep(
                  1,
                  length(track_list)
               ),
               chromosome = "chrNA",
               from = view_from,
               to = view_to,
               panel.only = TRUE,
               showTitle = FALSE,
               background.panel = "white"
            )

            # Use the same horizontal scale as the annotation tracks
            # The vertical scale gives every hit one equal-height row
            pushViewport(viewport(
               xscale = c(
                  view_from,
                  view_to
               ),
               yscale = c(
                  0,
                  n_hits
               ),
               clip = "off"
            ))

            # Add endpoint coordinates to hits from every database
            for (hit_index in seq_len(n_hits)) {
               # The first hit is plotted in the top row
               # The 0.72 offset aligns text with the arrow body
               track_y <- n_hits -
                  hit_index +
                  0.72

               # Show the reference coordinate corresponding to query_from
               grid.text(
                  label = format(
                     subset_data$reference_from[hit_index],
                     trim = TRUE,
                     scientific = FALSE
                  ),
                  x = unit(
                     subset_data$query_from[hit_index],
                     "native"
                  ) - unit(
                     1.5,
                     "mm"
                  ),
                  y = unit(
                     track_y,
                     "native"
                  ),
                  just = c(
                     "right",
                     "center"
                  ),
                  gp = gpar(
                     fontsize = 10,
                     col = "black"
                  )
               )

               # Show the reference coordinate corresponding to query_to
               grid.text(
                  label = format(
                     subset_data$reference_to[hit_index],
                     trim = TRUE,
                     scientific = FALSE
                  ),
                  x = unit(
                     subset_data$query_to[hit_index],
                     "native"
                  ) + unit(
                     1.5,
                     "mm"
                  ),
                  y = unit(
                     track_y,
                     "native"
                  ),
                  just = c(
                     "left",
                     "center"
                  ),
                  gp = gpar(
                     fontsize = 10,
                     col = "black"
                  )
               )
            }

            popViewport()
         } else {
            grid.text(
               "No matches",
               gp = gpar(
                  fontsize = 12
               )
            )
         }

         popViewport()
         popViewport(2)
         popViewport()
      }, finally = {
         dev.off()
      })

      # Count hits from each source for the output manifest
      database_counts <- table(
         factor(
            subset_data$database,
            levels = DATABASE_ORDER
         )
      )

      all_plots[[consensus_index]] <- data.frame(
         Consensus = consensus,
         length = consensus_length,
         n_matches = n_hits,
         repbase_matches = unname(
            database_counts[["RepBase"]]
         ),
         our_known_matches = unname(
            database_counts[["our-known"]]
         ),
         swissprot_matches = unname(
            database_counts[["Swiss-Prot"]]
         ),
         gb_te_matches = unname(
            database_counts[["GB_TE"]]
         ),
         retrovirus_matches = unname(
            database_counts[["Retrovirus"]]
         ),
         png_file = png_name,
         stringsAsFactors = FALSE
      )
   }

   plot_manifest <- do.call(
      rbind,
      all_plots
   )

   write.table(
      plot_manifest,
      file = file.path(
         output_dir,
         "plot_manifest.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
   )

   return(plot_manifest)
}


############################################
# PART 7: Read pipeline arguments and run
############################################

main <- function() {
   # Expected argument order:
   # 1 consensus FASTA
   # 2 12-column Censor map
   # 3 RepBase FASTA
   # 4 our-known FASTA
   # 5 Swiss-Prot WU-BLAST report
   # 6 GB_TE WU-BLAST report
   # 7 Retrovirus WU-BLAST report
   # 8 output directory
   # 9 comma-separated consensus IDs or all
   arguments <- commandArgs(
      trailingOnly = TRUE
   )

   if (length(arguments) != 9L) {
      stop(
         paste(
            "[ERROR] Usage:",
            "Rscript src/plot_ins_map.R",
            "<consensus_fasta>",
            "<censor_map>",
            "<repbase_fasta>",
            "<our_known_fasta>",
            "<swissprot_wublast>",
            "<gb_te_wublast>",
            "<retrovirus_wublast>",
            "<output_dir>",
            "<consensus_ids>"
         ),
         call. = FALSE
      )
   }

   consensus_fasta <- arguments[[1]]
   map_file <- arguments[[2]]
   repbase_file <- arguments[[3]]
   our_known_file <- arguments[[4]]
   swissprot_file <- arguments[[5]]
   gb_te_file <- arguments[[6]]
   retrovirus_file <- arguments[[7]]
   output_dir <- arguments[[8]]
   consensus_ids_text <- arguments[[9]]

   # Split the ninth argument at literal commas
   # Example: "family1,family2" becomes c("family1", "family2")
   selected_consensus_ids <- strsplit(
      consensus_ids_text,
      split = ",",
      fixed = TRUE
   )[[1]]

   selected_consensus_ids <- trimws(
      selected_consensus_ids
   )

   if (any(!nzchar(selected_consensus_ids))) {
      stop(
         "[ERROR] Selected consensus IDs must not be empty",
         call. = FALSE
      )
   }

   if (anyDuplicated(selected_consensus_ids)) {
      stop(
         "[ERROR] Selected consensus IDs must not contain duplicates",
         call. = FALSE
      )
   }

   consensus_lengths <- read_consensus_lengths(
      consensus_fasta
   )

   # Replace all with every consensus ID from the FASTA
   if (identical(selected_consensus_ids, "all")) {
      selected_consensus_ids <- consensus_lengths$Consensus
   }

   missing_consensus_ids <- setdiff(
      selected_consensus_ids,
      consensus_lengths$Consensus
   )

   if (length(missing_consensus_ids) > 0L) {
      stop(
         "[ERROR] Selected consensus ID is absent from FASTA: ",
         paste(missing_consensus_ids, collapse = ", "),
         call. = FALSE
      )
   }

   # Keep selected consensus records in the requested order
   consensus_lengths <- consensus_lengths[
      match(
         selected_consensus_ids,
         consensus_lengths$Consensus
      ),
      ,
      drop = FALSE
   ]

   map_hits <- read_map_file(
      file_path = map_file,
      repbase_file = repbase_file,
      our_known_file = our_known_file
   )

   swissprot_hits <- read_wublast_hits(
      file_path = swissprot_file,
      database = "Swiss-Prot"
   )

   gb_te_hits <- read_wublast_hits(
      file_path = gb_te_file,
      database = "GB_TE"
   )

   retrovirus_hits <- read_wublast_hits(
      file_path = retrovirus_file,
      database = "Retrovirus"
   )

   # Every reader returns the same columns in the same order
   annotation_data <- rbind(
      map_hits,
      swissprot_hits,
      gb_te_hits,
      retrovirus_hits
   )

   annotation_data <- annotation_data[
      annotation_data$Consensus %in% selected_consensus_ids,
      ,
      drop = FALSE
   ]

   plot_manifest <- generate_plots(
      data = annotation_data,
      consensus_lengths = consensus_lengths,
      output_dir = output_dir,
      width_in = 14,
      dpi = 300
   )

   print(plot_manifest)
}


main()