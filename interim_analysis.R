### prepare_interim_data #######################################################
prepare_interim_data <-
  function(
    data,
    study_time,
    id_variable,
    covariates_variables,
    enrollment_time_variable,
    treatment_variable,
    outcome_variables,
    outcome_time_variables,
    observe_missing_times,
    outcomes_sequential = TRUE,
    time_to_event = FALSE
  ){
    missing_vars <-
      setdiff(
        x = c(id_variable, covariates_variables, enrollment_time_variable,
              outcome_variables, outcome_time_variables),
        y = names(data)
      )

    if(length(missing_vars) > 0){
      stop("Variables not contained in dataset: ",
           paste(missing_vars, collapse = ", "))
    }

    latest_event <-
      unlist(data[, c(enrollment_time_variable, outcome_time_variables)]) |>
      max(na.rm = TRUE)

    if(latest_event > study_time){
      stop("Latest enrollment/outcome event > `study_time`: ",
           latest_event, " > ", study_time)
    }

    n_outcomes <- length(outcome_variables)
    n_outcome_times <- length(outcome_time_variables)

    all_ids <- unique(data[, id_variable])

    if(!identical(x = n_outcomes, y = n_outcome_times)){
      stop("Number of outcomes (", n_outcomes, ") must be equal to the number ",
           "of outcome measurement times (", n_outcome_times, ").")
    }

    if(!identical(x = n_outcomes, y = length(observe_missing_times))){
      stop("Number of outcomes (", n_outcomes, ") must be equal to the length ",
           "of the times at which missingness is observed (",
           length(observe_missing_times), ").")
    }

    # Check data for potential errors: Missing or Duplicated IDs
    if(any(is.na(data[, id_variable]))){
      stop("ID variable `", id_variable, "` cannot contain NA values. ",
           "NA values observed in rows: ",
           paste(which(is.na(data[, id_variable])), collapse = ", "))
    }

    if(any(duplicated(data[, id_variable]))){
      stop("ID variable `", id_variable, "` cannot contain duplicates. ",
           "Duplicate values observed in rows: ",
           paste(which(duplicated(data[, id_variable])), collapse = ", "))
    }

    # Check for missing randomization times
    if(any(is.na(data[, enrollment_time_variable]))){
      stop(
        "Participants missing entry time in rows: ",
        which(is.na(data[, enrollment_time_variable])) |>
          paste(collapse = ", ")
      )
    }

    if(outcomes_sequential){
      # Check for inconsistencies in Outcome Timing
      outcome_time_differences <-
        data[, tail(x = outcome_time_variables, -1)] -
        data[, head(x = outcome_time_variables, -1)]

      inconsistencies <-
        which(x = outcome_time_differences < 0, arr.ind = TRUE) |>
        as.data.frame()

      if(nrow(inconsistencies) > 0){
        visits_affected <- inconsistencies$col
        msg <- "Inconsistent timing in outcomes: "
        for(i in visits_affected){
          rows_affected <- subset(x = inconsistencies, col == i)$row
          msg <-
            paste0(msg, "`", names(outcome_time_differences)[i], "`: rows (",
                   paste(rows_affected, collapse = ", "), "); ")
        }
      }
    }

    # Check for Outcomes with Missing Outcome Times
    missing_outcome_times <-
      which(
        x = !is.na(data[, outcome_variables]) &
          is.na(data[, outcome_time_variables]),
        arr.ind = TRUE
      ) |>
      data.frame()

    if(nrow(missing_outcome_times) > 0){
      stop("Missing outcome times in rows: ",
           paste0(unique(missing_outcome_times$row), collapse = ", "))
    }

    wide_data <-
      list(
        ".id" = data[, id_variable],
        data[, c(covariates_variables, treatment_variable)]
      )

    enrollment_time <- wide_data[[".e"]] <- data[, enrollment_time_variable]

    enrollment_to_study_time <- study_time - enrollment_time

    for(i in 1:n_outcomes){
      observed_time_i <- data[, outcome_time_variables[i]]
      outcome_i <- data[, outcome_variables[i]]

      missing_outcome_times_i <-
        which(is.na(observed_time_i) &
                enrollment_to_study_time > observe_missing_times[i])

      observed_time_i[missing_outcome_times_i] <-
        enrollment_time[missing_outcome_times_i] + observe_missing_times[i]

      enrollment_to_outcome <- observed_time_i - enrollment_time

      outcome_indicator <- rep(NA, length(outcome_i))
      outcome_indicator[which(!is.na(outcome_i))] <- TRUE
      outcome_indicator[
        which(is.na(outcome_i) & !is.na(observed_time_i))
      ] <- FALSE
      outcome_indicator[
        which(enrollment_to_outcome > observe_missing_times[i] &
                is.na(outcome_i))
      ] <- FALSE

      observed_time_i[
        which(
          is.na(observed_time_i) &
            study_time - enrollment_time > observe_missing_times[i]
        )
      ]

      wide_data[[outcome_variables[i]]] <- data[, outcome_variables[i]]
      wide_data[[paste0(".t_", i)]] <- observed_time_i
      wide_data[[paste0(".r_", i)]] <- 1*outcome_indicator
    }

    return(
      list(
        data = do.call(what = cbind, args = wide_data),
        original_data = data,
        variables = list(
          id_variable = id_variable,
          covariates_variables = covariates_variables,
          enrollment_time_variable = enrollment_time_variable,
          treatment_variable = treatment_variable,
          outcome_variables = outcome_variables,
          outcome_time_variables = outcome_time_variables,
          observe_missing_times = observe_missing_times,
          outcomes_sequential = outcomes_sequential,
          time_to_event = time_to_event,
          renamed_outcome_times = paste0(".t_", 1:n_outcomes),
          outcome_indicators = paste0(".r_", 1:n_outcomes)
        ),
        study_time = study_time
      )
    )
  }




### data_at_time_t #############################################################
data_at_time_t <-
  function(
    prepared_data,
    study_time
  ){

    data <- prepared_data$data

    # 1. Subset to participants enrolled by study time
    data <- data[which(data$`.e` <= study_time),]

    outcome_times <- prepared_data$variables$renamed_outcome_times

    wide_data <-
      list(
        data[, c(".id", prepared_data$variables$covariates_variables,
                 prepared_data$variables$treatment_variable, ".e")]
      )

    enrollment_time <- data$`.e`

    outcome_variables <- prepared_data$variables$outcome_variables
    observe_missing_times <- prepared_data$variables$observe_missing_times

    for(i in 1:length(outcome_times)){
      observed_time <- data[, outcome_times[i]]
      enrollment_to_outcome <- observed_time - enrollment_time
      outcome_i <- outcome_i_at_time_t <- data[, outcome_variables[i]]

      not_yet_observed <-
        is.na(observed_time) | (observed_time > study_time)

      is_observed_complete <-
        !(is.na(outcome_i) | not_yet_observed)

      is_observed_missing <-
        is.na(outcome_i) & (!not_yet_observed) &
        (enrollment_to_outcome >= c(observe_missing_times)[i])

      is_not_yet_observed <-
        (
          not_yet_observed | # Observed after study_time threshold
            (is.na(outcome_i) & (!not_yet_observed) &
               (enrollment_to_outcome < c(observe_missing_times)[i]))
        )

      outcome_determined <-
        (is_observed_complete + is_observed_missing + is_not_yet_observed)

      outcome_undetermined <- !(outcome_determined %in% 1)

      if(any(outcome_undetermined)){
        stop("Undetermined outcomes in rows: ",
             paste0(which(outcome_undetermined), collapse = ", "))
      }

      outcome_indicator <- NA*is_observed_complete
      outcome_indicator[which(is_observed_complete)] <- 1
      outcome_indicator[which(is_observed_missing)] <- 0

      wide_data[[paste0(".r_", i)]] <- outcome_indicator

      observed_time[is_not_yet_observed] <- NA
      wide_data[[paste0(".t_", i)]] <- pmin(observed_time, study_time)

      outcome_i_at_time_t[which(is_not_yet_observed)] <- NA
      wide_data[[outcome_variables[i]]] <- outcome_i_at_time_t
    }

    original_data = prepared_data$original_data
    variables = prepared_data$variables

    return(
      list(
        data = do.call(what = cbind, args = wide_data),
        original_data = original_data,
        variables = variables,
        study_time = study_time
      )
    )
  }




### count_outcomes #############################################################
count_outcomes <-
  function(
    prepared_data,
    study_time = NULL
  ){

    if(!is.null(study_time)){
      prepared_data <-
        data_at_time_t(
          prepared_data = prepared_data,
          study_time = study_time
        )
    }

    data <- prepared_data$data

    data <-
      data[rank(x = data$`.e`, ties.method = "first"),]

    id_variable <- prepared_data$variables$id_variable
    outcome_variables <- prepared_data$variables$outcome_variables
    outcome_time_variables <- prepared_data$variables$renamed_outcome_times

    binary_outcome <-
      all(unlist(data[, outcome_variables]) %in% c(NA, 0, 1))

    long_times <-
      list(
        data.frame(
          id = data[, id_variable],
          event = "randomization",
          value = NA,
          time = data$`.e`,
          count_total =
            rank(
              x = data$`.e`,
              na.last = "keep",
              ties = "first"
            )
        )
      )

    long_times[[1]]$count_complete <- long_times[[1]]$count_total

    if(binary_outcome){
      long_times[[1]]$count_events <- NA
    }


    for(i in 1:length(outcome_variables)){
      complete_times <- data[, outcome_time_variables[i]]

      if(binary_outcome){
        event_times <- complete_times
        event_times[which(!(data[, outcome_variables[i]] %in% 1))] <- NA
      }

      complete_times[which(is.na(data[, outcome_variables[i]]))] <- NA

      event_times_i <-
        data.frame(
          id = data[, id_variable],
          event = outcome_variables[i],
          value = data[, outcome_variables[i]],
          time = data[, outcome_time_variables[i]],
          count_total =
            rank(
              x = data[, outcome_time_variables[i]],
              na.last = "keep",
              ties = "first"
            ),
          count_complete =
            rank(
              x = complete_times,
              na.last = "keep",
              ties = "first"
            )
        )

      if(binary_outcome){
        event_times_i$count_events <-
          rank(
            x = event_times,
            na.last = "keep",
            ties = "first"
          )
      }

      long_times[[length(long_times) + 1]] <-
        event_times_i[order(event_times_i$time),]
    }

    long_times <-
      do.call(
        what = rbind,
        args = long_times
      )

    long_times$event <-
      factor(
        x = long_times$event,
        levels = c("randomization", outcome_variables)
      )

    long_times[, c("id", "value")] <- NULL

    return(
      subset(
        long_times,
        !is.na(time)
      )
    )
  }




### plot_outcome_counts ########################################################
plot_outcome_counts <-
  function(
    prepared_data,
    study_time = NULL,
    type = "tc",
    count_increment = 10,
    time_increment = 30,
    color_palette = NULL
  ) {

    if(!is.null(study_time)){
      prepared_data <-
        data_at_time_t(
          prepared_data = prepared_data,
          study_time = study_time
        )
      study_time <- prepared_data$study_time
    }

    type <- tolower(type)
    plot_total <- length(grep(pattern = "t", x =  type)) > 0
    plot_complete <- length(grep(pattern = "c", x =  type)) > 0
    plot_events <- length(grep(pattern = "e", x =  type)) > 0

    plot_elements <- plot_total + plot_complete + plot_events
    if(plot_elements < 1){
      stop("`type` must contain at least one of the following characters: ",
           "'t' (total observations), 'c' (complete observations), or 'e' ",
           "(event count).")
    }

    outcome_counts <-
      count_outcomes(
        prepared_data = prepared_final_data,
        study_time = study_time
      )

    if(plot_events){
      if(!("count_events" %in% names(outcome_counts))){
        stop("`type` contains 'e' (event counts), which is for time-to-event and ",
             "binary outcomes, but `count_outcomes()` did not detect a ",
             "time-to-event or binary outcome.")
      }
    }

    events <- levels(outcome_counts$event)

    if(length(color_palette) < length(events)){
      if(!is.null(color_palette)){
        warning("Length of `color_palette (", length(color_palette), ") must be ",
                "at least ", length(events), ".")
      }
      color_palette <- 1:length(events)
    }

    y_max <-
      if(plot_total){
        max(outcome_counts$count_total, na.rm = TRUE)
      } else if(plot_complete){
        max(outcome_counts$count_complete, na.rm = TRUE)
      } else if(plot_complete){
        max(outcome_counts$count_events, na.rm = TRUE)
      }

    par(mar = c(5, 4, 4, 4) + 0.25)

    plot(
      x = NULL,
      y = NULL,
      xlim = range(outcome_counts$time),
      ylim = c(0, y_max),
      xlab = "Study Time",
      ylab = "Count"
    )

    abline(
      h = seq(from = 0, to = y_max, by = count_increment),
      col = gray(level = 0, alpha = 0.125)
    )

    abline(
      v = seq(from = 0, to = max(outcome_counts$time), by = time_increment),
      col = gray(level = 0, alpha = 0.125)
    )

    for(i in 1:length(events)){
      plot_i_data <-
        subset(
          x = outcome_counts,
          event == events[i]
        )

      if(plot_total){
        lines(
          count_total ~ time,
          data = plot_i_data[order(plot_i_data$count_total),],
          col = i,
          type = "S"
        )
      }

      if(plot_complete){
        lines(
          count_complete ~ time,
          data = plot_i_data[order(plot_i_data$count_complete),],
          col = i,
          lty = 2,
          type = "S"
        )
      }

      if(plot_events){
        lines(
          count_events ~ time,
          data = plot_i_data[order(plot_i_data$count_events),],
          col = i,
          lty = 3,
          type = "S"
        )
      }
    }

    axis(side = 4, las = 1)

    legend_text <- legend_lty <- NULL

    if(plot_complete) {
      legend_text <- c(legend_text, "Complete")
      legend_lty <- c(legend_lty, 1)
    }
    if(plot_total) {
      legend_text <- c(legend_text, "Non-missing")
      legend_lty <- c(legend_lty, 2)
    }
    if(plot_events) {
      legend_text <- c(legend_text, "Events")
      legend_lty <- c(legend_lty, 3)
    }

    legend(
      x ="bottomright",
      legend = legend_text,
      lty = legend_lty,
      inset = c(0, 1),
      xpd = TRUE,
      horiz = TRUE,
      bty = "n"
    )

    events[1] <- "RND"

    legend(
      x = "topleft",
      legend = events,
      col = color_palette,
      lty = 1,
      lwd = 1.75
    )

  }




### count_outcomes_at_time_t ###################################################
count_outcomes_at_time_t <-
  function(
    prepared_data,
    study_times
  ){

    data <- prepared_data$data
    outcome_time_variables <- prepared_data$variables$renamed_outcome_times
    outcome_variables <- prepared_data$variables$outcome_variables

    event_names <- c("randomization", prepared_data$variables$outcome_variables)
    count_total <-
      data.frame(
        times = study_times
      )
    count_total[, event_names] <- NA
    count_complete <- count_total

    for(i in 1:length(study_times)){
      for(j in 1:length(event_names)){
        if(j == 1){
          count_total[i, event_names[j]] <-
            count_complete[i, event_names[j]] <-
            sum(data$`.e` <= study_times[i], na.rm = TRUE)
        } else {
          count_total[i, event_names[j]] <-
            sum(data[, outcome_time_variables[j - 1]] <=
                  study_times[i], na.rm = TRUE)
          count_complete[i, event_names[j]] <-
            sum(data[, outcome_time_variables[j - 1]] <= study_times[i] &
                  !is.na(data[, outcome_variables[j - 1]]), na.rm = TRUE)
        }
      }
    }

    return(
      list(
        count_total = count_total,
        count_complete = count_complete
      )
    )
  }




### calculate_estimate #########################################################
calculate_estimate <-
  function(
    data,
    estimation_function,
    estimation_arguments
  ){
    return(
      get(
        x = "estimate",
        pos =
          do.call(
            what = estimation_function,
            args = c(list(data = data), estimation_arguments)
          )
      )
    )
  }




### calculate_covariance #######################################################
calculate_covariance <-
  function(
    data,
    bootstrap_ids = NULL,
    bootstrap_results = NULL,
    estimation_function,
    estimation_arguments,
    rng_seed,
    control = monitored_design_control()
  ){

    current_ids <- unique(data$.id)
    n_current_ids <- length(current_ids)

    ### Issue: Stratify bootstrap on observed outcome pattern
    if(is.null(bootstrap_ids)){
      n_k <- length(unique(data$.id))
      ids_to_sample <- n_k
    } else {
      n_k <-
        c(sapply(X = bootstrap_ids[[1]], FUN = length),
          length(unique(data$.id)))
      ids_to_sample <- tail(x = diff(x = n_k), n = 1)

      if(ids_to_sample == 0){
        stop("No new IDs from prior dataset: Please check `.id` column.")
      }
    }

    set.seed(seed = rng_seed, kind = "L'Ecuyer")
    random_seed <- .Random.seed
    k <- length(n_k)

    # Use separate RNG stream for each interim analysis
    if(k > 1){
      for(i in 2:k) parallel::nextRNGStream(seed = random_seed)
    }

    ### Issue: Stratify bootstrap on observed outcome pattern
    # Generate resampled IDs
    new_bootstrap_ids <-
      matrix(
        data = NA,

        ncol = control$n_bootstrap,
        nrow = ids_to_sample,
      ) |>
      apply(
        MARGIN = 2,
        FUN = function(x)
          current_ids[
            sample(
              x = 1:n_current_ids,
              size = ids_to_sample,
              replace = TRUE
            )
          ]
      )

    # Combine previous and current bootstrap IDs
    all_bootstrap_ids <-
      rbind(
        sapply(
          X = bootstrap_ids,
          FUN = function(x) do.call(what = c, args = x)
        ),
        new_bootstrap_ids
      )

    # Create cluster, export objects
    cluster <- parallel::makeCluster(control$n_cores)

    # Run estimation function on bootstrapped data
    apply_function <-
      function(
    x,
    original_data = data,
    args = estimation_arguments,
    fun = estimation_function,
    run_function = calculate_estimate
      ) {
        boot_ids <- x[[1]]

        run_function(
          data =
            original_data[match(x = boot_ids, table = original_data$.id),],
          estimation_function = fun,
          estimation_arguments = args
        )
      }

    parallel::clusterExport(
      cl = cluster,
      varlist = c("data", "calculate_estimate",
                  "estimation_function", "estimation_arguments",
                  "apply_function"),
      envir = environment()
    )

    # Choose appropriate call depending on whether load balancing is desired
    par_apply_function <-
      ifelse(
        test = control$use_load_balancing,
        yes = parallel::parSapplyLB,
        no = parallel::parSapply
      )

    # Compute estimates in parallel
    estimates <-
      do.call(
        what = par_apply_function,
        args =
          list(
            cl = cluster,
            X =
              apply(
                X = all_bootstrap_ids,
                MARGIN = 2,
                FUN = list
              ),
            FUN = apply_function
          )
      ) |>
      as.numeric()

    # Stop the cluster
    parallel::stopCluster(cl = cluster)

    if(is.null(bootstrap_ids)){
      bootstrap_ids <- list()
      length(bootstrap_ids) <- control$n_bootstrap
    }

    # Add bootstrap IDs for analysis k
    for(i in 1:length(bootstrap_ids)){
      bootstrap_ids[[i]][[k]] <- new_bootstrap_ids[, i]
    }

    estimates <- cbind(bootstrap_results, estimates)

    return(
      list(
        covariance = cov(estimates),
        bootstrap_estimates = estimates,
        bootstrap_ids = bootstrap_ids
      )
    )
  }




### orthogonalize_estimates ####################################################
orthogonalize_estimates <-
  function(
    estimates,
    covariance
  ){
    k <- length(estimates)

    projection <-
      matrix(data = covariance[k, k], nrow = k, ncol = k) -
      kronecker(
        X = matrix(data = as.vector(c(covariance[1:k-1, k], 0)), nrow = 1),
        Y = matrix(data = 1, nrow = k)
      ) -
      t(
        kronecker(
          Y = matrix(data = 1, nrow = k),
          X = matrix(data = as.vector(c(covariance[1:k-1, k], 0)), nrow = 1)
        )
      ) +
      rbind(cbind(covariance[1:(k-1), 1:(k-1)], 0), 0)

    # A is the Cholesky decomposition of covMatrixProjection
    A <- chol(projection)
    # Atilde is a kx(k-1) matrix with the first k-1 columns of A
    Atilde = A[, -k]
    # Ak is a k x 1 vector equal to the kth column of A
    Ak = A[, k]

    # W equals the vector # (theta_k-theta_1, ..., theta_k-theta_{k-1}, theta_k)
    W = c(estimates[k] - estimates[1:k-1], estimates[k])

    # scale are the orthogonalizing 'coefficients'
    scale <-
      solve(A)%*%(diag(k) - Atilde%*%solve(t(Atilde)%*%Atilde)%*%t(Atilde))%*%Ak

    # updatedEstimate is the updated/orthogonalized estimate at analysis k
    estimate_orthogonal <- t(scale)%*%W

    # updatedVariance is the variance of the
    # updated/orthogonalized estimate at analysis k
    variance_orthogonal <- t(scale)%*%t(A)%*%A%*%scale

    return(
      list(
        estimate_orthogonal = estimate_orthogonal,
        covariance_orthogonal = variance_orthogonal
      )
    )
  }




### interim_decision_asymptotic ################################################
interim_decision_asymptotic <-
  function(
    test_statistics,
    trial_design,
    information_fraction,
    information_target
  ){

    k <- length(test_statistics)
    k_max <- trial_design$kMax

    if(k < k_max){
      # Change information level from a priori value to observed value
      information_rates <- trial_design$informationRates
      information_rates[k] <- information_fraction[k]

      trial_design_params_interim <-
        c("kMax", "alpha", "beta", "sided",
          "typeOfDesign", "gammaA","typeBetaSpending", "gammaB")

      trial_design_updated <-
        do.call(
          what = rpact::getDesignGroupSequential,
          args =
            c(as.list(trial_design)[trial_design_params_interim],
              list(informationRates = information_rates))
        )

    } else if (k == k_max){
      information_rates <- information_fraction*information_target
      information_rates <- information_rates/tail(x = information_rates, 1)

      trial_design_params_final <-
        c("kMax", "alpha", "beta", "sided", "typeBetaSpending", "gammaB")

      trial_design_updated <-
        do.call(
          what = rpact::getDesignGroupSequential,
          args =
            c(as.list(trial_design)[trial_design_params_final],
              list(typeOfDesign = "asUser",
                   informationRates = information_rates,
                   userAlphaSpending = trial_design$alphaSpent))
        )
    }

    bounds_efficacy <- trial_design_updated$criticalValues
    bounds_futility <- trial_design_updated$futilityBounds

    futility_bound_k <-
      ifelse(
        test = length(bounds_futility[k]) > 0,
        yes = bounds_futility[k],
        no = -Inf
      )

    if(test_statistics[k] < futility_bound_k){
      continue <- FALSE
      decision <- "Futility"
      decision_detail <- paste0("Futility - stage ", k, " of ", k_max)
      stopping_stage <- k
    } else if(test_statistics[k] > bounds_efficacy[k]){
      continue <- FALSE
      decision <- "Efficacy"
      decision_detail <- paste0("Efficacy - stage ", k, " of ", k_max)
      stopping_stage <- k
    } else if(k < k_max){
      continue <- TRUE
      decision <- "Continue"
      decision_detail <- paste0("Continue - stage ", k, " of ", k_max)
      stopping_stage <- NA
    } else {
      continue <- FALSE
      decision <- "Fail to reject"
      decision_detail <- paste0("Fail to reject: stage ", k, " of ", k_max)
      stopping_stage <- k_max
    }


    if(length(test_statistics) < k_max){
      test_statistics <-
        c(test_statistics, rep(NA, k_max - length(test_statistics)))
    }

    return(
      list(
        continue = continue,
        decision = decision,
        decision_detail = decision_detail,
        stopping_stage = stopping_stage,
        trial_design_updated = trial_design_updated,
        decision_data =
          data.frame(
            test_statistic = test_statistics,
            efficacy = bounds_efficacy,
            futility = c(bounds_futility, NA)
          )
      )
    )
  }




### monitored_design_checks ####################################################
monitored_design_checks <-
  function(
    data,
    monitored_design = NULL,
    information_target,
    estimation_function,
    estimation_arguments,
    correction_function = NULL
  ) {
    if(!(".id" %in% names(data))){
      stop("Data must contain a column `.id` identifying all rows belonging ",
           "to each individual.")
    }

    prior_analysis <- tail(x = monitored_design, 1)[[1]]

    if(is.null(trial_design) & is.null(monitored_design)){
      stop("For an initial analysis, `trial_design` must be specified. For ",
           "subsequent analyses, `monitored_design` must be specified.")
    } else if(!(is.null(trial_design) | is.null(monitored_design))){
      previous_design <- prior_analysis$trial_design
      compare_params <-
        c("kMax", "alpha", "beta", "sided",
          "typeOfDesign", "gammaA","typeBetaSpending", "gammaB")
      if(
        !identical(
          x = as.list(trial_design)[compare_params],
          y = as.list(previous_design)[compare_params]
        )
      ){
        stop("`trial_design` must be identical to `trial_design` in the latest ",
             "analysis in `monitored_design` for the following parameters: ",
             paste(compare_params, collapse = ", "))
      }
    } else if(is.null(trial_design) & !is.null(monitored_design)){
      trial_design <- prior_analysis$trial_design
    }

    # Check previous analyses against current design
    if(!is.null(monitored_design)){

      previous_decisions <-
        sapply(
          X = monitored_design[-1],
          FUN = function(x) get(x = "decision", pos = x)
        )

      if(any(!previous_decisions %in% c("continue"))){
        if(any(previous_decisions %in% c("futility"))){
          warning("Futility stopping boundary reached at analyses: ",
                  paste(which(previous_decisions == "futility"), collapse = ", ")
          )
        }

        if(any(previous_decisions %in% c("efficacy"))){
          stop("Efficacy stopping boundary reached at analyses: ",
               paste(which(previous_decisions == "efficacy"), collapse = ", "))
        }

        if(any(previous_decisions %in% c("fail to reject"))){
          stop("All pre-specified analyses completed: Failed to reject null ",
               "hypothesis at final analysis.")
        }
      }

      previous_designs <-
        lapply(
          X = monitored_design,
          FUN = function(y) get(x = "trial_design", pos = y)
        )

      k <-
        setdiff(
          x = names(monitored_design),
          y = "original_design"
        ) |> length()

      # Get information rates: Arrange each into column of matrix
      informationRates <-
        sapply(
          X = previous_designs,
          FUN = function(y) get(x = "informationRates", pos = y)
        )

      if(is.list(informationRates)){
        stop("Information rates should have the same length across all ",
             "analyses. `informationRates` in `monitored_design` have ",
             "lengths: ",
             paste0(sapply(X = informationRates, FUN = length), collapse = ", ")
        )
      }

      # Check consistency of observed information
      if(k > 1){
        for(i in 1:(k - 1)){
          if(length(unique(tail(x = informationRates[i, ], -i))) > 1){
            stop("Inconsistent observed information fractions for interim ",
                 "analysis ", i, ".")
          }
        }
      }


      # Check consistency of design parameters
      check_params <-
        c("alpha", "beta", "sided", "typeOfDesign",
          "typeBetaSpending", "bindingFutility")

      previous_params <-
        sapply(
          X = previous_designs,
          FUN = function(x)
            unlist(mget(x = check_params, envir = as.environment(x)))
        )

      previous_param_conflicts <-
        previous_params |>
        apply(
          MARGIN = 1,
          FUN = unique
        ) |>
        sapply(
          FUN = length
        )

      if(any(previous_param_conflicts > 1)){
        stop("Conflicting values in analysis parameters for ",
             paste(paste0("`", names(which(previous_param_conflicts > 1)), "`"),
                   collapse = ", "), ".")
      }

      if(k > 1){
        # Check consistency of design parameters
        check_results <-
          c("estimates", "estimatesOrthogonal",
            "variance", "varianceOrthogonal",
            "testStatistics", "testStatisticsOrthogonal",
            "decision", "decisionOrthogonal")

        # Determine which results are contained in the object
        check_results <-
          lapply(X = tail(x = monitored_design, -1), FUN = names) |>
          unlist() |>
          intersect(y = check_results)

        previous_results <-
          lapply(
            X = tail(x = monitored_design, -1),
            FUN = function(x)
              mget(x = check_results, envir = as.environment(x))
          )

        for(i in check_results){
          extracted_result <-
            sapply(
              X = previous_results,
              FUN =
                function(x, get_result = i, new_length = k){
                  extracted_result <-
                    get(x = i, pos = x)
                  length(extracted_result) <- new_length
                  return(extracted_result)
                }
            )

          if(any(!is.na(extracted_result[lower.tri(extracted_result)]))){
            stop("Inconsistency in monitored_design element `", result, "`.")
          }

          unique_results <-
            apply(
              X = extracted_result,
              MARGIN = 1,
              FUN = function(x) length(unique(na.omit(object = x)))
            )

          if(any(unique_results > 1)){
            stop("Inconsistency in monitored_design element `", result, "`.")
          }
        }
      }
    }
  }




### estimate_information #######################################################
estimate_information <-
  function(
    data,
    monitored_design,
    estimation_function,
    estimation_arguments,
    orthogonalize = NULL,
    rng_seed,
    return_results = FALSE,
    control = monitored_design_control()
  ) {

    # If previous analyses specified, assign most recent result to trial_design
    if(!is.null(monitored_design)){
      k <- setdiff(x = names(monitored_design), y = "original_design") |>
        length() + 1

      prior_analysis <- tail(monitored_design, 1)[[1]]

      extract_params <-
        c("trial_design", "orthogonalize", "estimates", "covariance",
          "estimates_orthogonal", "covariance_orthogonal",
          "bootstrap_ids", "bootstrap_results",
          "information", "information_target", "rng_seed")

      for(i in 1:length(extract_params)){
        assign(
          x = extract_params[i],
          value = prior_analysis[[extract_params[i]]]
        )
      }

      rm(prior_analysis)

    } else {
      k <- 1
      dataset_list <- list(data)
      bootstrap_ids <- NULL
      bootstrap_results <- NULL
      estimates <- NULL
      covariance <- NULL
      estimates_orthogonal <- NULL
      covariance_orthogonal <- NULL
    }

    n_analyses <- trial_design$kMax
    if(k == n_analyses){
      analysis_label <- "final_analysis"
    } else{
      analysis_label <- paste0("interim_analysis_", k)
    }

    k <- length(estimates) + 1

    estimate_k <-
      calculate_estimate(
        data = data,
        estimation_function = estimation_function,
        estimation_arguments = estimation_arguments
      )

    covariance_k <-
      calculate_covariance(
        data = data,
        bootstrap_ids = bootstrap_ids,
        bootstrap_results = bootstrap_results,
        estimation_function = estimation_function,
        estimation_arguments = estimation_arguments,
        rng_seed = rng_seed,
        control = control
      )

    estimates <- as.numeric(c(estimates, estimate_k))
    covariance <- covariance_k$covariance
    information <- 1/diag(covariance)

    # Add to Previous Estimates
    if(orthogonalize == TRUE){
      if(k == 1) {
        estimates_orthogonal <- estimate_k
        covariance_orthogonal <- covariance_k$covariance
      } else {
        orthogonalized_k <-
          orthogonalize_estimates(
            estimates = estimates,
            covariance = covariance_k$covariance
          )

        estimates_orthogonal[k] <- orthogonalized_k$estimate_orthogonal

        covariance_orthogonal <-
          c(diag(covariance_orthogonal), orthogonalized_k$covariance_orthogonal)

        covariance_orthogonal <-
          sqrt(covariance_orthogonal) %*% t(sqrt(covariance_orthogonal))
      }

      information_orthogonal <- 1/diag(covariance_orthogonal)
    } else {
      estimates_orthogonal <- NA*estimate
      covariance_orthogonal <- NA*covariance$covariance
    }

    if(return_results){
      bootstrap_ids <- covariance_k$bootstrap_ids
      bootstrap_results <- covariance_k$bootstrap_estimates
    } else {
      estimates <- estimates_orthogonal <-
        bootstrap_ids <- bootstrap_results <- NULL
    }

    return(
      list(
        estimates = estimates,
        covariance = covariance,
        information = information,
        estimates_orthogonal = estimates_orthogonal,
        covariance_orthogonal = covariance_orthogonal,
        information_orthogonal = information_orthogonal,
        bootstrap_results = bootstrap_results,
        bootstrap_ids = bootstrap_ids,
        analysis_label = analysis_label
      )
    )
  }




### monitored_analysis #########################################################
monitored_analysis <-
  function(
    data,
    monitored_design = NULL,
    estimation_function,
    estimation_arguments,
    correction_function = NULL,
    control = monitored_design_control(),
    ...
  ){

    all_args <- c(as.list(environment()), list(...))

    # Check for valid input
    do.call(
      what = monitored_design_checks,
      args =
        all_args[intersect(x = names(all_args),
                           y = formalArgs(monitored_design_checks))]
    )

    if(is.null(monitored_design)) {
      stop("Initialize an information monitored design with ",
           "initialize_monitored_design()")
    } else {
      information_target <- monitored_design$original_design$information_target
      trial_design <- tail(x = monitored_design, n = 1)[[1]]$trial_design
      orthogonalize <- monitored_design$original_design$orthogonalize
    }

    estimated_information <-
      estimate_information(
        data = data,
        monitored_design = monitored_design,
        estimation_function = estimation_function,
        estimation_arguments = estimation_arguments,
        orthogonalize = orthogonalize,
        rng_seed = monitored_design$original_design$rng_seed,
        control = control,
        return_results = TRUE
      )

    for(i in 1:length(estimated_information)){
      assign(
        x = names(estimated_information)[i],
        value = estimated_information[[i]]
      )
    }

    information_fraction <- information/information_target
    information_fraction_orthogonal <- information_orthogonal/information_target

    z_statistics <- estimates/diag(covariance)
    z_statistics_orthogonal <- estimates_orthogonal/diag(covariance_orthogonal)

    if(orthogonalize){
      interim_decisions <-
        interim_decision_asymptotic(
          test_statistics = z_statistics_orthogonal,
          trial_design = trial_design,
          information_fraction = information_fraction_orthogonal,
          information_target = information_target
        )
    } else {
      interim_decisions <-
        interim_decision_asymptotic(
          test_statistics = z_statistics,
          trial_design = trial_design,
          information_fraction = information_fraction,
          information_target = information_target
        )
    }

    for(i in 1:length(interim_decisions)){
      assign(
        x = names(interim_decisions)[i],
        value = interim_decisions[[i]]
      )
    }

    trial_design_updated <- interim_decisions$trial_design

    analysis_results <-
      c(monitored_design,
        setNames(
          object =
            list(
              list(
                data = data,
                trial_design = trial_design_updated,
                information_target = information_target,
                estimates = estimates,
                covariance = covariance,
                z_statistics = z_statistics,
                information = information,
                information_fraction = information_fraction,
                estimates_orthogonal = estimates_orthogonal,
                covariance_orthogonal = covariance_orthogonal,
                z_statistics_orthogonal = z_statistics_orthogonal,
                information_orthogonal = information_orthogonal,
                information_fraction_orthogonal =
                  information_fraction_orthogonal,
                bootstrap_results = bootstrap_results,
                bootstrap_ids = bootstrap_ids,
                continue = continue,
                decision = decision,
                decision_detail = decision_detail,
                decision_data = decision_data,
                session_info = sessionInfo()
              )
            ),
          nm = estimated_information$analysis_label
        )
      )

    class(analysis_results) <- "monitored_design"

    return(analysis_results)
  }




### information_trajectory #####################################################
information_trajectory <-
  function(
    prepared_data,
    monitored_design = NULL,
    estimation_function,
    estimation_arguments,
    orthogonalize,
    n_min = 30,
    n_increment = 5,
    rng_seed,
    control = monitored_design_control()
  ){
    outcome_counts <-
      count_outcomes(
        prepared_data =
          data_at_time_t(
            prepared_data = prepared_data,
            study_time = prepared_data$study_time
          )
      )

    primary_outcome <-
      tail(x = prepared_data$variables$outcome_variables, n = 1)

    primary_outcome_counts <-
      subset(
        x = outcome_counts,
        event == primary_outcome
      )

    n_current <- max(primary_outcome_counts$count_complete, na.rm = TRUE)

    information_times <-
      subset(
        x = primary_outcome_counts,
        count_complete %in%
          seq(
            from = n_min,
            to = n_current,
            by = n_increment
          )
      )$time

    information <-
      data.frame(
        count_outcomes_at_time_t(
          prepared_data = prepared_data,
          study_times = information_times
        )$count_complete,
        information = NA
      )

    for(i in 1:nrow(information)){
      data_i <-
        data_at_time_t(
          prepared_data = prepared_data,
          study_time = information$times[i]
        )

      information_i <-
        estimate_information(
          data = data_i$data,
          monitored_design = monitored_design,
          estimation_function = estimation_function,
          estimation_arguments = estimation_arguments,
          orthogonalize = orthogonalize,
          rng_seed = rng_seed,
          control = control,
          return_results = FALSE
        )

      information$information[i] <-
        if(orthogonalize){
          information_i$information_orthogonal
        } else {
          information_i$information
        }
    }

    information$information_lag_1 <-
      with(information, c(NA, head(x = information, n = -1)))

    information$information_change =
      with(information, information - information_lag_1)

    information$information_pct_change =
      with(information, 100*information_change/information)

    return(information)
  }




### initialize_monitored_design ################################################
initialize_monitored_design <-
  function(
    trial_design,
    maximum_sample_size,
    information_target,
    orthogonalize,
    rng_seed_analysis
  ){
    monitored_design <-
      setNames(
        object = list(
          list(
            trial_design = trial_design,
            maximum_sample_size = maximum_sample_size,
            information_target = information_target,
            orthogonalize = orthogonalize,
            rng_seed_analysis = rng_seed_analysis,
            estimates = NULL,
            covariance = NULL,
            information = NULL,
            information_fraction = NULL,
            estimates_orthogonal = NULL,
            covariance_orthogonal = NULL,
            bootstrap_results = NULL,
            bootstrap_ids = NULL
          )
        ),
        nm = "original_design"
      )
    class(monitored_design) <- "monitored_design"
    return(monitored_design)
  }




### estimate_final_sample_size #################################################
estimate_final_sample_size <-
  function(
    data,
    monitored_design,
    estimation_function,
    estimation_arguments,
    control = monitored_design_control()
  ){

    if(is.null(monitored_design)) {
      stop("Initialize an information monitored design with ",
           "initialize_monitored_design()")
    } else {
      information_target <- monitored_design$original_design$information_target
      trial_design <- tail(x = monitored_design, n = 1)[[1]]$trial_design
    }

    # Get dataset of participants who have completed the study
    outcome_indicators <-
      names(data)[grep(pattern = "^\\.r_\\d{1,}$", x = names(data))]

    outcome_indicators <-
      outcome_indicators[
        gsub(pattern = "^\\.r_", replacement = "", x = outcome_indicators) |>
          as.numeric() |> order()
      ]

    final_indicator <- tail(x = outcome_indicators, n = 1)

    data_completed <-
     data[which(!is.na(data[final_indicator])),]

    completed_n <- nrow(data_completed)

    completed_information <-
      estimate_information(
        data = data_completed,
        monitored_design = monitored_design,
        estimation_function = estimation_function,
        estimation_arguments = estimation_arguments,
        return_results = FALSE,
        control = control
      )

    information <-
      ifelse(
        test = monitored_design$original_design$orthogonalize,
        yes = completed_information$information_orthogonal,
        no = completed_information$information
      )

    n_max_design <-
      monitored_design$original_design$maximum_sample_size
    n_max_information <-
      ceiling(x = completed_n*information_target/information)

    n_final <-
      min(n_max_design, n_max_information)

    return(
      list(
        n_final = n_final,
        observed_outcomes = completed_n,
        information_fraction = information/information_target,
        exceed_maximum = n_max_information > n_max_design
      )
    )
  }




### monitored_design_control ###################################################
monitored_design_control <-
  function(
    n_bootstrap = 1000,
    n_cores = 1,
    use_load_balancing = FALSE
  ){
    return(
      list(
        n_bootstrap = n_bootstrap,
        n_cores = n_cores,
        use_load_balancing = use_load_balancing
      )
    )
  }




### summary.monitored_design ###################################################
summary.monitored_design <-
  function(
    object
  ){
    analysis_names <- names(object)
    stage <- length(analysis_names) - 1
    total_stages <- object$original_design$trial_design$kMax
    information_targets <-
      with(
        data = object$original_design,
        expr = {
          round(
            x = information_target*(trial_design$informationRates),
            digits = 4
          )
        }
      )
    orthogonalize <- object$original_design$orthogonalize

    latest_analysis <- tail(x = object, n = 1)[[1]]

    current_information <-
      get(
        x = ifelse(
          test = orthogonalize,
          yes = "information_orthogonal",
          no = "information"
        ),
        pos = latest_analysis
      )

    current_information_fraction <-
      get(
        x = ifelse(
          test = orthogonalize,
          yes = "information_fraction_orthogonal",
          no = "information_fraction"
        ),
        pos = latest_analysis
      )

    cat("Information Monitoring Design: Stage", stage, "of", total_stages, "\n")
    cat("Decision at latest analysis:", latest_analysis$decision_detail,"\n")
    cat("Maximum Sample Size:", object$original_design$maximum_sample_size, "\n")
    cat("Information:\n")
    cat("\tTargets:", paste(information_targets, collapse = ", "), "\n")
    cat("\tLatest Information Level:", current_information, "\n")
    cat("\tLatest Information Fraction:", current_information_fraction, "\n")
  }




### regression_wrapper #########################################################
regression_wrapper <-
  function(
    data = data,
    regression_function,
    formula,
    ...
  ){
    args <- as.list(substitute(list(...)))[-1L]

    fitted_regression <-
      do.call(
        what = regression_function,
        args =
          c(list(), args,
            list(
              formula = formula,
              data = data
            )
          )
      )

    return(fitted_regression)
  }




### inverse_regression #########################################################
inverse_regression <-
  function(
    data = data,
    indices = NULL,
    y_target,
    formula = formula,
    smooth = FALSE,
    smooth_kernel = c("normal", "uniform")[1],
    ...
  ) {

    if(!is.null(indices)){
      if(smooth){
        n <- nrow(data)

        if(smooth_kernel == "normal"){
          smooth_sds <-
            diag(x = diag(cov(data))/sqrt(n))

          noise <-
            mvtnorm::rmvnorm(
              n = n,
              mean = c(0, 0),
              sigma = t(smooth_sds) %*% cov2cor(cov(data)) %*% smooth_sds
            )
        } else if(smooth_kernel == "uniform"){
          noise <-
            matrix(
              data =
                runif(
                  n = 2*n,
                  min = -1/sqrt(n),
                  max = 1/sqrt(n)
                ),
              nrow = n,
              ncol = 2
            )
        }

        data <- data[indices, ] + noise
      } else {
        data <- data[indices, ]
      }
    }

    fitted_regression <-
      regression_wrapper(data = data, formula = formula, ...)

    coefs <- as.numeric(coef(fitted_regression))

    return(
      c("estimate" =
          (y_target - coefs[1])/coefs[2],
        "(intercept)" = coefs[1],
        "slope" = coefs[2]
      )
    )
  }




### inverse_regression_bca #####################################################
inverse_regression_bca <-
  function(
    data,
    formula,
    regression_function = deming::theilsen,
    y_target,
    n_bootstrap = 10000,
    confidence = c(0.95),
    smooth = FALSE,
    ...
  ){

    data <- data[, all.vars(expr = formula)]

    inverse_regression_boot <-
      boot::boot(
        data = data,
        statistic = inverse_regression,
        R = n_bootstrap,
        regression_function = regression_function,
        formula = formula,
        y_target = y_target,
        smooth = smooth,
        ...
      )

    fitted_regression <-
      regression_wrapper(
        data = data,
        formula = formula,
        regression_function = regression_function,
        ...
      )

    intervals <-
      data.frame(
        expand.grid(
          i = 1:length(inverse_regression_boot$t0),
          parameter = NA, estimate = NA, lcl = NA, ucl = NA,
          confidence = confidence
        )
      )

    for(i in 1:nrow(intervals)){
      intervals$parameter[i] <-
        c("target", "intercept", "slope")[intervals$i[i]]
      intervals$estimate[i] <-
        inverse_regression_boot$t0[intervals$i[i]]
      intervals[i, c("lcl", "ucl")] <-
        tail(
          x =
            boot::boot.ci(
              boot.out = inverse_regression_boot,
              conf = intervals$confidence[i],
              type = "bca",
              index = intervals$i[i]
            )$bca[1,],
          n = 2
        )
    }

    intervals$i <- NULL

    return(
      list(
        fitted_regression = fitted_regression,
        ci = intervals,
        boot_object = inverse_regression_boot
      )
    )
  }
