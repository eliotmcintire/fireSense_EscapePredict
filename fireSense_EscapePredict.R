# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_EscapePredict",
  description = paste("Predicts fire escape probabilities using a model fitted with",
                      "the `fireSense_EscapeFit` module. These can be used to feed the",
                      "fire escape component of a landscape fire model (e.g., fireSense)."),
  keywords = c("escape probability", "fire frequency", "logistic", "fireSense"),
  authors = person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.0", fireSense_EscapePredict = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_EscapePredict.Rmd"),
  loadOrder = list(after = "fireSense_dataPrepPredict"),
  reqdPkgs = list("magrittr", "raster", "stats"),
  parameters = bindrows(
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = paste("optional. Interval between two runs of this module",
                                 "expressed in units of simulation time. By default, 1 year.")),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA,
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA,
                    desc = "optional. Interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated? This is generally intended",
                                 "for data-type modules where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "fireSense_EscapeFitted", objectClass = "fireSense_EscapeFit",
                 desc = "An object of class fireSense_EscapeFit created with the fireSense_EscapeFit module."),
    expectsInput(objectName = "fireSense_IgnitionAndEscapeCovariates", objectClass = "data.frame",
                 desc = "An object of class RasterStack or data.frame with prediction variables"),
    expectsInput(objectName = "flammableRTM", objectClass = "SpatRaster",
                 desc = "a raster with values of 1 for flammable pixels and 0 for nonflammable pixels")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "fireSense_EscapePredicted", objectClass = "SpatRaster",
                  desc = "A raster with values representing escape probability")
  )
))
doEvent.fireSense_EscapePredict = function(sim, eventTime, eventType, debug = FALSE) {

  switch(
    eventType,
    init = {
      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, "fireSense_EscapePredict", "run",
                           eventPriority = 5.12)

      if (!is.na(P(sim)$.saveInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "fireSense_EscapePredict", "save", .last())
    },
    run = {
      sim <- escapePredictRun(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, "fireSense_EscapePredict", "run",
                             eventPriority = 5.12)
    },
    save = {
      sim <- escapePredictSave(sim)

      if (!is.na(P(sim)$.saveInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval,"fireSense_EscapePredict", "save", .last())
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )

  invisible(sim)
}

escapePredictRun <- function(sim) {
  
  EscapeCovariates <- copy(sim$fireSense_IgnitionAndEscapeCovariates)
  moduleName <- currentModule(sim)
  
  ## Toolbox: set of functions used internally by escapePredictRun
  escapePredictRaster <- function(model, data, sim) {
    model %>%
      model.matrix(data) %>%
      `%*%` (coef(sim$fireSense_EscapeFitted)) %>%
      drop %>% sim$fireSense_EscapeFitted$family$linkinv(.)
  }
  
  # Load inputs in the data container
  # list2env(as.list(envir(sim)), envir = mod)
  
  mod_env <- new.env(parent = baseenv())
  list2env(EscapeCovariates, env = mod_env)
  
  terms <- delete.response(terms.formula(sim$fireSense_EscapeFitted$formula))
  
  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)
  
  if (all(unlist(lapply(allxy, function(x) is.vector(mod_env[[x]]))))) {
    mod_env[["predictEscapeFun"]] <- function() {
      formula %>%
        model.matrix(mod_env) %>%
        `%*%` (coef(sim$fireSense_EscapeFitted)) %>%
        drop %>% sim$fireSense_EscapeFitted$family$linkinv(.)
    }
  }
  else if (all(unlist(lapply(allxy, function(x) is(mod_env[[x]], "SpatRaster"))))) {
    mod_env[["predictEscapeFun"]] <- function() {
      mget(allxy, envir = mod_env, inherits = FALSE) %>%
        stack %>% raster::predict(model = formula, fun = escapePredictRaster, na.rm = TRUE, sim = sim)
    }
  }
  else {
    missing <- !allxy %in% ls(mod_env, all.names = TRUE)
    
    if (s <- sum(missing))
      stop(moduleName, "> '", allxy[missing][1L], "'",
           if (s > 1) paste0(" (and ", s-1L, " other", if (s>2) "s", ")"),
           " not found in data objects.")
    
    badClass <- unlist(lapply(allxy, function(x) is.vector(mod_env[[x]]) || is(mod_env[[x]], "SpatRaster")))
    
    if (any(badClass)) {
      stop(moduleName, "> Data objects of class 'data.frame' cannot be mixed with objects of other classes.")
    } else    {
      stop(moduleName, "> '", paste(allxy[which(!badClass)], collapse = "', '"),
           "' does not match a data.frame's column, a SpatRaster or a layer from a RasterStack or RasterBrick.")
    }
  }
  
  sim$fireSense_EscapePredicted <- mod_env[["predictEscapeFun"]]()
  
  #force EscapePredicted a raster
  if (!inherits(sim$fireSense_EscapePredicted, "SpatRaster")){
    EscapeRas <- rast(sim$flammableRTM)
    EscapeRas[EscapeCovariates$pixelID] <- sim$fireSense_EscapePredicted
    sim$fireSense_EscapePredicted <- EscapeRas
  }
  
  return(invisible(sim))
}

.inputObjects <- function(sim){
  
  if (!suppliedElsewhere("fireSense_IgnitionAndEscapeCovariates", sim) & 
      !suppliedElsewhere("fireSense_EscapeFitted", sim) &
      !suppliedElsewhere("flammableRTM", sim)) {
    
    ignitions <- c(490, 14, 2)
    escapes <- c(498, 7, 1)
    class2 <- runif(n = 506, min = 0, max = 1) * 0.55
    class3 <- runif(n = 506, min = 0, max = 1) * 0.45
    nonForest_highFlammable <- 1 - c(class2 + class3)
    pseudoData <- data.table(ignitions = rep(0:2, times = ignitions),
                             escapes = rep(0:2, times = escapes), 
                             class2 = class2, class3 = class3, 
                             NF_HF = nonForest_highFlammable)
    pseudoData[, MDC := runif(506, 120, 200)]
    sim$fireSense_IgnitionAndEscapeCovariates <- pseudoData
    sim$fireSense_IgnitionAndEscapeCovariates$pixelID <- sample(1:600, size = 506, replace = FALSE)
    fireSense_EscapeFormula <- "cbind(escapes, ignitions - escapes) ~ MDC + NF_HF + class2 + class3 - 1"
    sim$fireSense_EscapeFitted <- glm(formula = as.formula(fireSense_EscapeFormula), 
                                      data = sim$fireSense_IgnitionAndEscapeCovariates, family = "binomial")
    class(sim$fireSense_EscapeFitted) <- c("fireSense_EscapeFit", "glm", "lm")
    flammableRTM <- rast(nrow = 30, ncol = 20)
    flammableRTM[sim$fireSense_IgnitionAndEscapeCovariates$pixelID] <- 1
    sim$flammableRTM <- flammableRTM
  }

  #TODO: ideally this would retrieve necessary inputObjects  
  return(invisible(sim))
}
