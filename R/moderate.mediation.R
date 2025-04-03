#' NEWWS Riverside data
#' 
#' The National Evaluation of Welfare-to-Work Strategies (NEWWS) was conducted before the nationwide welfare reform in the mid-1990s. 
#' Participants of the study consist of applicants to the Aid to Families with Dependent Children (AFDC) program and current AFDC recipients who were not working for 30 or more hours per week. 
#' They were randomly assigned to the labor force attachment (LFA) program, which aimed at moving low-income parents from welfare to work by providing employment-focused incentives and services, and the control group, which received aid from AFDC without requirement for working. 
#' The data is from the Riverside sample. It includes 694 participants with preschool-age children, 208 randomly assigned to the LFA program and 486 randomly assigned to the control group. 
#' The National Evaluation of Welfare-to-Work Strategies (NEWWS) Riverside study. 
#' @docType data
#' @name newws
#' @return A list containing
#' \item{treat}{Treatment is a binary variable, taking the value of 1 if a participant was assigned to LFA and 0 otherwise}
#' \item{emp}{Mediator is a binary variable, taking the value of 1 if one was ever employed during the two-year period after randomization and 0 if not.}
#' \item{depression}{Outcome, maternal depression at the end of the second year after randomization, is a summary score of 12 items measuring depressive symptoms during the past week on a 0-3 scale.}
#' \item{nevmar}{1 if never married and 0 otherwise.}
#' \item{emp_prior}{1 if employed and 0 otherwise.}
#' \item{hispanic}{1 if Hispanic and 0 otherwise.}
#' \item{CHCNT}{1 if had 1 child, 2 if had 2 children, and 3 if had 3 or more children before randomization.}
#' \item{ADCPC}{welfare amount in the year before randomization}
#' \item{attitude}{a composite score of two attitude items - so many family problems that I cannot work at a full time or part time job; so much to do during the day that I cannot go to a school or job training program - measured on the scale of 1-4.}
#' \item{depress_prior}{a composite score of three depressive symptom items - sad, depressed, blues, and lonely - in the week before randomization measured on the scale of 1-4.}
#' \item{workpref}{one's level of preference for taking care of family full time than working on the scale of 1-4.}
#' \item{nohsdip}{1 if had never obtained a high school diploma or a General Educational Development certificate and 0 otherwise.}
NULL

#' Causal Moderated Mediation Analysis
#' 
#' 'modmed' is used to fit mediator and outcome models and estimate and test causal effects for causal moderated mediation analysis. It is applicable to a treatment of any scale, a binary or continuous mediator and outcome, one or more moderators of any scale, and a wide range of scenarios of moderated mediation. 
#' @param data A data frame containing the variables in the analysis. Users need to impute any missing values in the data before running the function.
#' @param treatment A character string indicating the name of the treatment variable. If a treatment is binary, please code it as 0 and 1.
#' @param mediator A character string indicating the name of the mediator variable.
#' @param outcome A character string indicating the name of the outcome variable.
#' @param moderators.disc A name vector of all the discrete moderators under study (string). The default is NULL.  
#' @param moderators.cont A name vector of all the continuous moderators under study (string). The default is NULL.
#' @param covariates.disc A name vector of the discrete pretreatment covariates that are not moderators (string). The default is NULL. Users do not need to reformat discrete variables. If treatment is randomized, it should contain confounders of the mediator-outcome relationship. If treatment is not randomized, it should also contain confounders of the treatment-mediator and treatment-outcome relationships. Please do not include moderators here. The main effects of moderators are accounted for by including the moderators in the model of the main model intercept. 
#' @param covariates.cont A name vector of the continuous pretreatment covariates that are not moderators (string). The default is NULL. If treatment is randomized, it should contain confounders of the mediator-outcome relationship. If treatment is not randomized, it should also contain confounders of the treatment-mediator and treatment-outcome relationships. Please do not include moderators here. The main effects of moderators are accounted for by including the moderators in the model of the main model intercept. 
#' @param m.model A list. The name of each element in the list is a predictor in the main mediator model.  The names must include intercept, treatment, covariates.disc, and covariates.cont. Each element of the list is a vector of the names of the moderators (string) of the coefficient of the main model predictor as represented by the name of the element. Each moderator specified in moderators.disc and moderators.cont must moderate at least one slope in either the mediator model or the outcome model. The moderators of the intercept must include all the moderators in the mediator model. If a main model coefficient is not moderated, then the corresponding vector should be specified as NULL. The set of moderators in m.model and that in y.model are not necessarily the same. In other words, for a moderator that moderates coefficient(s) in one model, it is possible that it does not occur in the other model. It is also possible that it only has its main effect in the other model, in which case it will be only included in the vector of the intercept in the list of the other model. The union of the two sets of moderators should be contained in the union of moderators.disc, and moderators.cont. If no moderators are specified, only the population average mediation effects are estimated and tested. Check examples for details.
#' @param y.model A list. The name of each element in the list is a predictor in the main outcome model. The names must include intercept, treatment, mediator, covariates.disc, and covariates.cont. If the treatment is assumed to interact with the mediator when affecting the outcome, an additional element should be added to y.model, named as "tm". Each element of the list is a vector of the names of the moderators (string) of the coefficient of the main model predictor as represented by the name of the element. Each moderator specified in moderators.disc and moderators.cont must moderate at least one slope in either the main mediator model or the main outcome model. The moderators of the intercept must include all the moderators in the mediator model. If a main model coefficient is not moderated, then the corresponding vector should be specified as NULL. The set of moderators in m.model and that in y.model are not necessarily the same. In other words, for a moderator that moderates coefficient(s) in one model, it is possible that it does not occur in the other model. It is also possible that it only has its main effect in the other model, in which case it will be only included in the vector of the intercept in the list of the other model. The union of the two sets of moderators should be contained in the union of moderators.disc, and moderators.cont. If no moderators are specified, only the population average mediation effects are estimated and tested. Check examples for details.
#' @param comp.treatment.value If the total treatment effect for each individual is defined as Y(t) - Y(t'), comp.treatment.value refers to t. The default is 1.
#' @param ref.treatment.value If the total treatment effect for each individual is defined as Y(t) - Y(t'), ref.treatment.value refers to t'. The default is 0.
#' @param comp.mod.disc.values A vector of compare values of the discrete moderators given which the conditional causal effects are estimated. The default is NULL.  The length and order of comp.mod.disc.values should be the same as moderators.disc. If one does not want to condition some moderators on specific values, one may specify their values to be NA. If the discrete moderators take the same compare value, comp.moderator.disc.values can be specified as a single value. If not NULL, results contrasting the causal effects given compare moderator values with those given reference moderator values will be reported.
#' @param ref.mod.disc.values A vector of reference values of the discrete moderators given which the conditional causal effects are estimated. The default is NULL. The length and order of ref.mod.disc.values should be the same as moderators.disc. If one does not want to condition some moderators on specific values, one may specify their values to be NA. If the discrete moderators take the same reference value, ref.moderator.disc.values can be specified as a single value. The moderators whose reference values are specified must be included in either the mediator model or the outcome model or both.
#' @param comp.mod.cont.values A vector of compare values of the continuous moderators given which the conditional causal effects are estimated. The default is NULL. The length and order of comp.mod.cont.values should be the same as moderators.cont. If one does not want to condition some moderators on specific values, one may specify their values to be NA. If the continuous moderators take the same compare value, comp.moderator.cont.values can be specified as a single value. If not NULL, results contrasting the causal effects given compare moderator values with those given reference moderator values will be reported. 
#' @param ref.mod.cont.values A vector of reference values of the continuous moderators given which the conditional causal effects are estimated. The default is NULL. The length and order of ref.mod.cont.values should be the same as moderators.cont. If one does not want to condition some moderators on specific values, one may specify their values to be NA. If the continuous moderators take the same reference value, ref.moderator.cont.values can be specified as a single value. The moderators whose reference values are specified must be included in either the mediator model or the outcome model or both.
#' @param m.scale A character string indicating the scale of the mediator. "continuous" if the mediator is continuous. "binary" if the mediator is binary. Probit link is used in this case. The default is "continuous".
#' @param y.scale A character string indicating the scale of the outcome. "continuous" if the outcome is continuous. "binary" if the outcome is binary. Probit link is used in this case. The default is "continuous".
#' @param method Estimation and inference method. if 'mc', the Monte Carlo method is used; if 'boot', bootstrap is used. Default is 'mc'. When sample size is relatively small and the mediator or binary is binary, 'boot' is preferred. Otherwise, 'mc' is more recommended, mainly because its running speed is much faster.
#' @param nmc Number of simulations involved in the Monte Carlo algorithm. Used if method = 'mc'. The default value is 1000.
#' @param nboot Number of bootstrapped samples involved in the bootstrapping algorithm. Used if method = 'boot'. The default value is 1000.
#' @param conf.level Level of the returned two-sided confidence intervals. The default is 0.95, which returns the 2.5 and 97.5 percentiles. 
#' @param seed The seed for the random number generator. The default value is NULL.
#' @param object Output from the original analysis only if is.U = TRUE. NULL in the original analysis. The default is NULL.
#' @param is.U A logical value. 'FALSE' in the original analysis. 'TRUE' in the sensitivity analysis which adjusts for a simulated unmeasured pretreatment confounder U. The default is FALSE.
#' @param sens.effect The default is NULL. It is used only when modmed.sens is run.
#' @param b.m Sensitivity parameter that represents the conditional association between the unobserved pretreatment confounder and the mediator. It is NULL in the original analysis. The default is NULL. It is used only when modmed.sens is run.
#' @param b.y Sensitivity parameter that represents the conditional association between the unobserved pretreatment confounder and the outcome. It is NULL in the original analysis. The default is NULL. It is used only when modmed.sens is run.
#' @param plot.effect The effect to be plotted. The default is NULL. It is used only when modmed.plot is run.
#' 
#' @return The output of the function modmed contains a list containing the following elements. The function 'summary' or 'summary_modmed' can be used to obtain a table of the results.
#' \item{effects}{Estimation results of the causal effects. "TE" indicates the total treatment effect. "TIE" indicates the total indirect effect, "PDE" indicates the pure direct effect, and "INT" indicates the natural treatment-by-mediator interaction effect. "TE.ref", "TIE.ref", "PDE.ref", and "INT.ref" indicate the corresponding effects when the moderators take the reference values. "TE.dif", TIE.dif", "PDE.dif", and "INT.dif" each indicates the difference in the corresponding effect between the compare levels and the reference levels of the moderators.}
#' \item{m.model}{Estimation results of the mediator model}
#' \item{y.model}{Estimation results of the outcome model}
#' \item{results}{1000 draws from the sampling distribution of the causal effects.}
#' \item{args}{A list of the arguments specified in the function of modmed except for the default ones}
#' \item{args.full}{The full list of all the arugments specified in the function of modmed}
#' \item{l.m}{The mediator model}
#' \item{l.y}{The outcome model}
#' \item{formula.m}{The formula for the mediator model}
#' \item{formula.y}{The formula for the outcome model}
#' \item{m.predictors}{All the original predictors in the main mediator model}
#' \item{y.predictors}{All the original predictors in the main outcome model}
#' \item{m.moderators}{All the original moderators in the mediator model}
#' \item{y.moderators}{All the original moderators in the outcome model}
#' \item{m.predictors.new}{All the predictors in the main mediator model, in which the discrete confounders are transformed into multiple indicator codes.}
#' \item{y.predictors.new}{All the predictors in the main outcome model, in which the discrete confounders are transformed into multiple indicator codes.}
#' \item{m.moderators.new}{All the moderators in the mediator model, in which the discrete confounders are transformed into multiple indicator codes.}
#' \item{y.moderators.new}{All the moderators in the outcome model, in which the discrete confounders are transformed into multiple indicator codes.}
#' \item{predict.m.data.ref}{Data for prediction of the conditional potential mediator at the reference level of the moderators.}
#' \item{predict.m.data.comp}{Data for prediction of the conditional potential mediator at the compare level of the moderators.}
#' \item{predict.m.data}{Data for prediction of the marginal potential mediator}
#' \item{predict.y.data.ref}{Data for prediction of the conditional potential outcome at the reference level of the moderators.}
#' \item{predict.y.data.comp}{Data for prediction of the conditional potential outcome at the compare level of the moderators.}
#' \item{predict.y.data}{Data for prediction of the marginal potential outcome}
#' \item{data}{The data used in the analysis}
#' @author Xu Qin and Lijuan Wang
#' @references Qin, X., & Wang, L. (2023). Causal moderated mediation analysis: Methods and software
#' @export
#' @importFrom stats vcov quantile as.formula coef fitted glm lm predict model.matrix 
#' @importFrom mvtnorm rmvnorm
#' @importFrom earth earth
#' @examples
#' data(newws)
#' modmed.results = modmed(data = newws, 
#'                         treatment = "treat", 
#'                         mediator = "emp", 
#'                         outcome = "depression", 
#'                         covariates.disc = NULL,
#'                         covariates.cont = NULL, 
#'                         moderators.disc = "CHCNT", 
#'                         moderators.cont = NULL, 
#'                         m.model = list(intercept = "CHCNT", treatment = "CHCNT"), 
#'                         y.model = list(intercept = "CHCNT", treatment = "CHCNT", mediator = "CHCNT", tm = "CHCNT"),
#'                         comp.mod.disc.values = 3, 
#'                         ref.mod.disc.values = 2, 
#'                         comp.mod.cont.values = NULL, 
#'                         ref.mod.cont.values = NULL, 
#'                         m.scale = "binary", 
#'                         y.scale = "continuous", 
#'                         seed = 1) 
#' \donttest{
#' modmed.results = modmed(data = newws, 
#'                         treatment = "treat", 
#'                         mediator = "emp", 
#'                         outcome = "depression", 
#'                         covariates.disc = c("emp_prior", "nevmar", "hispanic", "nohsdip"),
#'                         covariates.cont = c("workpref", "attitude", "depress_prior"), 
#'                         moderators.disc = "CHCNT", 
#'                         moderators.cont = "ADCPC", 
#'                         m.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), 
#'                         y.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), mediator = c("ADCPC", "CHCNT"), tm = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL),
#'                         comp.mod.disc.values = 3, 
#'                         ref.mod.disc.values = 2, 
#'                         comp.mod.cont.values = 5050, 
#'                         ref.mod.cont.values = 5050, 
#'                         m.scale = "binary", 
#'                         y.scale = "continuous", 
#'                         seed = 1) 
#' summary(modmed.results)
#'}

modmed = function(
    data,
    treatment,
    mediator,
    outcome,
    covariates.disc = NULL,
    covariates.cont = NULL,
    moderators.disc = NULL,
    moderators.cont = NULL,
    m.model, 
    y.model,
    comp.treatment.value = 1,
    ref.treatment.value = 0,
    comp.mod.disc.values = NULL,
    ref.mod.disc.values = NULL,
    comp.mod.cont.values = NULL,
    ref.mod.cont.values = NULL,
    m.scale = "continuous",
    y.scale = "continuous",
    method = "mc",
    nmc = 1000, 
    nboot = 1000,
    conf.level = 0.95,
    seed = NULL, 
    object = NULL,
    is.U = FALSE,
    sens.effect = NULL,
    b.m = NULL,
    b.y = NULL,
    plot.effect = NULL
){
  if(sum(is.na(data)))
    stop("Users need to impute any missing values in the data before running the function.")
  if(is.character(unique(data[, treatment])))
    stop("If a treatment is binary, please code it as 0 and 1.")
  data[, treatment] = as.numeric(as.character(data[, treatment])) # Otherwise, if the treatment is categorical, its name may change after model.matrix(), which will cause problems. Without as.character, the values would change.
  Unmeasure = NULL # Otherwise no visible binding for global variable
  
  model.covariates = c(names(m.model)[-which(names(m.model)%in%c("intercept", "treatment"))], names(y.model)[-which(names(y.model)%in%c("intercept", "treatment", "mediator", "tm"))])
  obs.covariates = c(covariates.disc, covariates.cont)
  if(any(!model.covariates%in%obs.covariates)|any(!obs.covariates%in%model.covariates))
  stop("The union of the covariates in m.model and y.model should be equal to the union of discrete and continuous covariates.")
  
  for(i in 1:length(covariates.disc)){
    if(length(unique(data[, covariates.disc[i]])) > 10)
      stop(paste("Please move", covariates.disc[i], "from covariates.disc to covariates.cont."))
  }
    
  confounders = c(covariates.disc, covariates.cont)
  moderators = c(moderators.disc, moderators.cont)
  if(any(moderators %in% confounders))
    stop("Please do not include moderators in covariates.disc or covariates.cont")
  ref.mod.values = c(ref.mod.disc.values, ref.mod.cont.values)
  comp.mod.values = c(comp.mod.disc.values, comp.mod.cont.values)
  
  if(!is.U){
    m.predictors = names(m.model)
    if(!"intercept" %in% m.predictors)
      stop("The element names in m.model must include 'Intercept'.")  
    m.predictors = m.predictors[-which(m.predictors == "intercept")]
    if(treatment %in% m.predictors)
      stop("Please specify the element name that corresponds to the treatment in m.model as 'treatment' rather than its actual variable name.")
    if(!"treatment" %in% m.predictors)
      stop("The element names in m.model must include 'treatment'. In other words, predictors in the main model of the mediator must include the treatment.")  
    m.predictors[which(m.predictors == "treatment")] = treatment
    
    y.predictors = names(y.model)
    if(!"intercept" %in% y.predictors)
      stop("The element names in y.model must include 'Intercept'.") 
    y.predictors = y.predictors[-which(y.predictors == "intercept")]
    if(treatment %in% y.predictors)
      stop("Please specify the element name that corresponds to the treatment in y.model as 'treatment' rather than its actual variable name.")  
    if(!"treatment" %in% y.predictors)
      stop("The element names in y.model must include 'treatment'. In other words, predictors in the main model of the outcome must include the treatment.")
    if(mediator %in% y.predictors)
      stop("Please specify the element name that corresponds to the mediator in y.model as 'mediator' rather than its actual variable name.") 
    if(!"mediator" %in% y.predictors)
      stop("The element names in y.model must include 'mediator'. In other words, predictors in the main model of the outcome must include the mediator.")
    y.predictors[which(y.predictors == "treatment")] = treatment
    y.predictors[which(y.predictors == "mediator")] = mediator
    
    if(any(moderators %in% m.predictors)|any(moderators %in% y.predictors)){
      stop("The element names in m.model or y.model cannot include moderators. In other words, please do not include moderators in the main model.")
    } else {
      if(any(!m.predictors %in% c(treatment, confounders)))
        stop("The pretreatment covariates that are not moderators in the mediator model (i.e., the element names that are not intercept or treatment in m.model) should be contained in the union of covariates.disc and covariates.cont.")
      if(any(!y.predictors %in% c(treatment, mediator, confounders))){
        if(any(!y.predictors %in% c(treatment, mediator, "tm", confounders)))
          stop("The pretreatment covariates that are not moderators in the outcome model (i.e., the element names that are not intercept, treatment, or mediator in y.model) should be contained in the union of covariates.disc and covariates.cont. If there is treatment-by-mediator interaction, it must be named as 'tm'.")
      }
    }
    
    if("tm" %in% y.predictors){
      y.predictors[which(y.predictors == "tm")] = paste0(treatment, ":", mediator)
    }   
    
    m.model = append(m.model[which(names(m.model) == "intercept")], m.model)
    m.model = m.model[-which(names(m.model) == "intercept")[2]]
    m.moderators = m.model
    y.model = append(y.model[which(names(y.model) == "intercept")], y.model)
    y.model = y.model[-which(names(y.model) == "intercept")[2]]
    y.moderators = y.model
    m.moderators.ori = unique(unlist(m.moderators))
    y.moderators.ori = unique(unlist(y.moderators))
    
    if(any(!unique(c(m.moderators.ori, y.moderators.ori)) %in% moderators)|any(!moderators %in% unique(c(m.moderators.ori, y.moderators.ori))))
      stop('The union of the two sets of moderators in m.model and y.model should be the same as the union of moderators.disc and moderators.cont.')
    
    if(!is.null(c(moderators.disc, moderators.cont))){
      for(i in 1:length(moderators)){
        if(!moderators[i] %in% unique(unlist(m.moderators)) & !moderators[i] %in% unique(unlist(y.moderators)))
          stop("A moderator specified in moderators.disc and moderators.cont must moderate at least one slope in either the mediator model or the outcome model.")
      }
    } else {
      if(!is.null(c(unique(unlist(m.moderators)), unique(unlist(y.moderators)))))
        stop("A moderator in the model needs to be specified in moderators.disc/moderators.cont.")
    }

    
    if(!all(is.null(unlist(m.moderators)))){
      if(!all(unique(unlist(m.moderators[-which(names(m.moderators) == "intercept")])) %in% m.moderators[[which(names(m.moderators) == "intercept")]]))
        stop("The moderators of the intercept of m.model must include all the moderators in the mediator model.")
    }
    if(!all(is.null(unlist(y.moderators)))){
      if(!all(unique(unlist(y.moderators[-which(names(y.moderators) == "intercept")])) %in% y.moderators[[which(names(y.moderators) == "intercept")]]))
        stop("The moderators of the intercept of y.model must include all the moderators in the outcome model.")
    }
    
    # if(is.null(ref.mod.disc.values) & !is.null(moderators.disc)){
    #   stop("If moderators.disc is not NULL, ref.mod.disc.values should not be NULL. The length and order of ref.mod.disc.values should be the same as moderators.disc. If one does not want to condition some moderators on specific values, one may specify their values to be NA.")
    # }
    if(!is.null(ref.mod.disc.values) & is.null(moderators.disc)){
      stop("If moderators.disc is NULL, ref.mod.disc.values should also be NULL. The length and order of ref.mod.disc.values should be the same as moderators.disc.")
    }
    # if(is.null(ref.mod.cont.values) & !is.null(moderators.cont)){
    #   stop("If moderators.cont is not NULL, ref.mod.cont.values should not be NULL. The length and order of ref.mod.cont.values should be the same as moderators.cont. If one does not want to condition some moderators on specific values, one may specify their values to be NA.")
    # }
    if(!is.null(ref.mod.cont.values) & is.null(moderators.cont)){
      stop("If moderators.cont is NULL, ref.mod.cont.values should also be NULL. The length and order of ref.mod.cont.values should be the same as moderators.cont.")
    }
    
    if(!is.null(ref.mod.disc.values)){
      if(length(ref.mod.disc.values) != length(moderators.disc))
        stop("Please make sure that the length and order of ref.mod.disc.values should be the same as moderators.disc") 
    }
    if(!is.null(ref.mod.cont.values)){
      if(length(ref.mod.cont.values) != length(moderators.cont))
        stop("Please make sure that the length and order of ref.mod.cont.values should be the same as moderators.cont")
    }
    
    if(is.null(ref.mod.disc.values) & !is.null(comp.mod.disc.values))
      stop("Please make sure to specify ref.mod.disc.values before specifying comp.mod.disc.values")
    if(is.null(ref.mod.cont.values) & !is.null(comp.mod.cont.values))
      stop("Please make sure to specify ref.mod.cont.values before specifying comp.mod.cont.values")
    
    if(!is.null(comp.mod.disc.values) & is.null(moderators.disc)){
      stop("If moderators.disc is NULL, comp.mod.disc.values should also be NULL. The length and order of comp.mod.disc.values should be the same as moderators.disc.")
    }
    if(!is.null(comp.mod.cont.values) & is.null(moderators.cont)){
      stop("If moderators.cont is NULL, comp.mod.cont.values should also be NULL. The length and order of comp.mod.cont.values should be the same as moderators.cont.")
    }
    
    if(!is.null(comp.mod.disc.values)){
      if(length(comp.mod.disc.values) != length(moderators.disc))
        stop("Please make sure that the length and order of comp.mod.disc.values should be the same as moderators.disc")
    }
    if(!is.null(comp.mod.cont.values)){
      if(length(comp.mod.cont.values) != length(moderators.cont))
        stop("Please make sure that the length and order of comp.mod.cont.values should be the same as moderators.cont")
    }
    
    if(!is.null(ref.mod.disc.values) & !any(is.na(ref.mod.disc.values))){
      for(i in 1:length(ref.mod.disc.values)){
        if(!ref.mod.disc.values[i] %in% unique(data[, moderators.disc[i]]))
          stop("Please specify each of ref.mod.disc.values to be one of the values of the corresponding moderator in the data")
        if(!is.null(comp.mod.disc.values)){
          if(!comp.mod.disc.values[i] %in% unique(data[, moderators.disc[i]]))
            stop("Please specify each of comp.mod.disc.values to be one of the values of the corresponding moderator in the data")
        } 
      }
    }
    
    if(!is.null(covariates.disc)){
      for(i in 1:length(covariates.disc)){
        if(!is.factor(data[, covariates.disc[i]]))
          data[, covariates.disc[i]] = as.factor(data[, covariates.disc[i]])
      }
    }  
    if(!is.null(moderators.disc)){
      for(i in 1:length(moderators.disc)){
        if(!is.factor(data[, moderators.disc[i]]))
          data[, moderators.disc[i]] = as.factor(data[, moderators.disc[i]])
      }
    }
    
    m.moderators.new = m.moderators
    y.moderators.new = y.moderators
    if(!is.null(m.moderators.new)){
      for(i in 1:length(m.moderators.new)){
        if(is.null(m.moderators.new[[i]]))
          m.moderators.new[[i]] = 1
      }
    }
    if(!is.null(y.moderators.new)){
      for(i in 1:length(y.moderators.new)){
        if(is.null(y.moderators.new[[i]]))
          y.moderators.new[[i]] = 1
      }
    }
    
    if(!is.null(m.moderators.new)){
      formula.m = paste(c(m.predictors, m.moderators.new[[1]]), collapse = " + ")
      for(i in 2:length(m.moderators.new))
        formula.m = paste(formula.m, "+", paste(m.moderators.new[[i]], ":", m.predictors[i - 1], collapse = " + "))
    } else {
      formula.m = paste(m.predictors, collapse = " + ")
    }
    
    if(!is.null(y.moderators.new)){
      formula.y = paste(c(y.predictors, y.moderators.new[[1]]), collapse = " + ")
      for(i in 2:length(y.moderators.new))
        formula.y = paste(formula.y, "+", paste(y.moderators.new[[i]], ":", y.predictors[i - 1], collapse = " + "))
    } else {
      formula.y = paste(y.predictors, collapse = " + ")
    }
    
    new.moderators.predictors = function(predictors, moderators.ori.list){
      if(!is.null(unique(unlist(moderators.ori.list)))){
        v1 = colnames(model.matrix(as.formula(paste("~", paste(moderators.ori.list[[1]], collapse = " + "))), data = data))[-1]
        k = 1
        for(i in 1:length(moderators.ori.list)){
          if(!is.null(moderators.ori.list[[i]]))
            moderators.ori.list[[i]] = colnames(model.matrix(as.formula(paste("~", paste(moderators.ori.list[[i]], collapse = " + "))), data = data))[-1]
        }
      }
      predictors.new = NULL
      for(i in 1:length(predictors)){
        predictors.new = c(predictors.new, colnames(model.matrix(as.formula(paste("~", paste(predictors[i], collapse = " + "))), data = data))[-1])
        if(!is.null(unique(unlist(moderators.ori.list)))){
          if(!predictors[i] %in% covariates.disc){
            assign(paste0("v", k + 1), moderators.ori.list[[i + 1]])
            k = k + 1
          }
          if(predictors[i] %in% covariates.disc){
            predictors.cat = colnames(model.matrix(as.formula(paste("~", paste(predictors[[i]], collapse = " + "))), data = data))[-1]
            for(j in 1:length(predictors.cat)){
              assign(paste0("v", k + j), moderators.ori.list[[i + 1]])
            }
            k = k + length(predictors.cat)
          }
        }
      }
      if(!is.null(unique(unlist(moderators.ori.list)))){
        moderators.new = rep(list(NULL), k)
        for(i in 1:k){
          if(!is.null(get(paste0("v", i))))
            moderators.new[[i]] = get(paste0("v", i))
        }
      }
      if(!is.null(unique(unlist(moderators.ori.list)))){
        return(list(moderators.new = moderators.new, predictors.new = predictors.new))
      } else {
        return(list(predictors.new = predictors.new))
      }
    }
    m.new.moderators.predictors = new.moderators.predictors(m.predictors, m.moderators)
    m.moderators.new = m.new.moderators.predictors$moderators.new
    m.predictors.new = m.new.moderators.predictors$predictors.new
    y.new.moderators.predictors = new.moderators.predictors(y.predictors, y.moderators)
    y.moderators.new = y.new.moderators.predictors$moderators.new
    y.predictors.new = y.new.moderators.predictors$predictors.new
    
    ref.m.mod.values = NULL
    comp.m.mod.values = NULL
    if(!is.null(m.moderators.ori)){
      for(i in 1:length(m.moderators.ori)){
        if(!is.null(ref.mod.values))
          ref.m.mod.values = c(ref.m.mod.values, ref.mod.values[which(moderators %in% m.moderators.ori[i])])
        if(!is.null(comp.mod.values))
          comp.m.mod.values = c(comp.m.mod.values, comp.mod.values[which(moderators %in% m.moderators.ori[i])])
      }
    }
    ref.y.mod.values = NULL
    comp.y.mod.values = NULL
    if(!is.null(y.moderators.ori)){
      for(i in 1:length(y.moderators.ori)){
        if(!is.null(ref.mod.values))
          ref.y.mod.values = c(ref.y.mod.values, ref.mod.values[which(moderators %in% y.moderators.ori[i])])
        if(!is.null(comp.mod.values))
          comp.y.mod.values = c(comp.y.mod.values, comp.mod.values[which(moderators %in% y.moderators.ori[i])])
      }
    }
    
    # Function for importing moderators.values to the data for prediction
    predict.data.mod = function(predict.data, predict.moderators, predict.moderators.values){
      if(!is.null(predict.moderators.values) & sum(is.na(predict.moderators.values)) > 0){
        predict.moderators = predict.moderators[-which(is.na(predict.moderators.values))]
        predict.moderators.values = predict.moderators.values[-which(is.na(predict.moderators.values))]
      }
      if(!is.null(predict.moderators.values)){
        for(i in 1:length(predict.moderators)){
          moderator = predict.moderators[i]
          if(moderator %in% colnames(predict.data)){
            predict.data[, moderator] = as.numeric(predict.moderators.values[i])
          } else {
            moderator.value = paste0(predict.moderators[i], predict.moderators.values[i])
            name.moderator.value = colnames(predict.data)[grep(moderator, colnames(predict.data))]
            if(moderator.value %in% colnames(predict.data)){
              predict.data[, moderator.value] = 1
              if(length(name.moderator.value) > 1){
                name.moderator.value.0 = name.moderator.value[which(!name.moderator.value %in% moderator.value)]
                predict.data[, name.moderator.value.0] = 0
              }
            } else {
              name.moderator.value.0 = name.moderator.value[which(!name.moderator.value %in% moderator.value)]
              predict.data[, name.moderator.value.0] = 0
            }
          }
        }
      }
      
      colnames(predict.data) = gsub("[^[:alnum:]\\:\\s]", "", colnames(predict.data))
      predict.data = as.data.frame(predict.data)
      
      return(predict.data)
    }
    
    predict.m.data = model.matrix(as.formula(paste("~", formula.m)), data = data)
    predict.y.data = model.matrix(as.formula(paste("~", formula.y)), data = data)
    
    predict.m.data.ref = predict.data.mod(predict.m.data, m.moderators.ori, ref.m.mod.values)
    predict.m.data.comp = predict.data.mod(predict.m.data, m.moderators.ori, comp.m.mod.values)
    predict.m.data = predict.data.mod(predict.m.data, m.moderators.ori, NULL)
    predict.y.data.ref = predict.data.mod(predict.y.data, y.moderators.ori, ref.y.mod.values)
    predict.y.data.comp = predict.data.mod(predict.y.data, y.moderators.ori, comp.y.mod.values)
    predict.y.data = predict.data.mod(predict.y.data, y.moderators.ori, NULL)
  } else {
    m.moderators.new = object$m.moderators.new
    y.moderators.new = object$y.moderators.new
    m.predictors.new = object$m.predictors.new
    y.predictors.new = object$y.predictors.new
    formula.m = object$formula.m
    formula.y = object$formula.y
    predict.m.data.ref = object$predict.m.data.ref
    predict.m.data.comp = object$predict.m.data.comp
    predict.m.data = object$predict.m.data
    predict.y.data.ref = object$predict.y.data.ref
    predict.y.data.comp = object$predict.y.data.comp
    predict.y.data = object$predict.y.data
  }
  
  # Fit mediator and outcome model based on the original data set first
  if(m.scale == "binary" & y.scale == "continuous"){
    # Step 1. Fit mediator and outcome models
    if(!is.U){
      l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data, family = binomial(link = "probit"))
      l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data)
    }
    if(is.U){
      l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data, family = binomial(link = "probit"))
      l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data)
    }
    
    if(any(is.na(coef(l.m)))){
      print(coef(l.m))
      stop("NA in coefficients of the mediator model. Please double check.")
    }
    if(any(is.na(coef(l.y)))){
      print(coef(l.y))
      stop("NA in coefficients of the outcome model. Please double check.")
    }
  }
  
  if(m.scale == "continuous" & y.scale == "continuous"){
    # Step 1. Fit mediator and outcome models
    if(!is.U){
      l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data)
      l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data)
    }
    if(is.U){
      l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data)
      l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data)
    }
    
    if(any(is.na(coef(l.m)))){
      print(coef(l.m))
      stop("NA in coefficients of the mediator model. Please double check.")
    }
    if(any(is.na(coef(l.y)))){
      print(coef(l.y))
      stop("NA in coefficients of the outcome model. Please double check.")
    }
  }
  
  if(m.scale == "binary" & y.scale == "binary"){
    # Step 1. Fit mediator and outcome models
    if(!is.U){
      l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data, family = binomial(link = "probit"))
      l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data, family = binomial(link = "probit"))
    }
    if(is.U){
      l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data, family = binomial(link = "probit"))
      l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data, family = binomial(link = "probit"))
    }
    
    if(any(is.na(coef(l.m)))){
      print(coef(l.m))
      stop("NA in coefficients of the mediator model. Please double check.")
    }
    if(any(is.na(coef(l.y)))){
      print(coef(l.y))
      stop("NA in coefficients of the outcome model. Please double check.")
    }
  }
  
  if(m.scale == "continuous" & y.scale == "binary"){
    # Step 1. Fit mediator and outcome models
    if(!is.U){
      l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data)
      sd.m = sigma(l.m) 
      l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data, family = binomial(link = "probit"))
    }
    if(is.U){
      l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data)
      sd.m = sigma(l.m) 
      l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data, family = binomial(link = "probit"))
    }
    
    if(any(is.na(coef(l.m)))){
      print(coef(l.m))
      stop("NA in coefficients of the mediator model. Please double check.")
    }
    if(any(is.na(coef(l.y)))){
      print(coef(l.y))
      stop("NA in coefficients of the outcome model. Please double check.")
    }
  }
  
  set.seed(seed)
  
  # Monte Carlo method
  if(method == 'mc'){
    if(m.scale == "binary" & y.scale == "continuous"){
      # Step 2. Simulate model parameters from their sampling distribution.
      coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
      coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
      
      # Step 3. Predict potential outcomes
      predict.m = function(predict.data, t){
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
        if(!is.U)
          return(pnorm(tcrossprod(predict.data, coef.m.sim)))
        if(is.U)
          return(pnorm((tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)))
      }
      
      predict.y = function(predict.data, t, m) {
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
        predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
        if(!is.U)
          return(tcrossprod(predict.data, coef.y.sim))
        if(is.U)
          return(tcrossprod(predict.data, coef.y.sim) + b.y * data$Unmeasure)
      }
      
      y1m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y1m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
      }
      y0m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y0m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
      }
    }
    
    if(m.scale == "continuous" & y.scale == "continuous"){
      # Step 2. Simulate model parameters from their sampling distribution.
      coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
      coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
      
      # Step 3. Predict potential outcomes
      predict.m = function(predict.data, t){
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
        if(!is.U)
          return(tcrossprod(predict.data, coef.m.sim))
        if(is.U)
          return(tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)
      }
      
      predict.y = function(predict.data, t, m) {
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        j.mediator = which(colnames(predict.data) == gsub("[^[:alnum:]\\:\\s]", "", mediator))
        predict.results = matrix(NA, nrow(predict.data), nmc)
        for(i in 1:nmc){
          predict.data[, j.mediator] = m[, i]
          if(!is.U)
            predict.results[, i] = tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ]))
          if(is.U)
            predict.results[, i] = tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ])) + b.y * data$Unmeasure
        }
        return(predict.results)
      }
      
      y1m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y1m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
      }
      y0m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y0m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
      }
    }
    
    if(m.scale == "binary" & y.scale == "binary"){
      # Step 2. Simulate model parameters from their sampling distribution.
      coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
      coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
      
      # Step 3. Predict potential outcomes
      predict.m = function(predict.data, t){
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
        if(!is.U)
          return(pnorm(tcrossprod(predict.data, coef.m.sim)))
        if(is.U)
          return(pnorm((tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)))
      }
      
      predict.y = function(predict.data, t, m) {
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
        predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
        if(!is.U){
          return(pnorm(tcrossprod(predict.data, coef.y.sim)))
        }
        if(is.U){
          return(pnorm((tcrossprod(predict.data, coef.y.sim) + b.y * data$Unmeasure)))
        }
      }
      
      y1m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y1m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
      }
      y0m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y0m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
      }
    }
    
    if(m.scale == "continuous" & y.scale == "binary"){
      # Step 2. Simulate model parameters from their sampling distribution.
      coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
      coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
      
      # Step 3. Predict potential outcomes
      predict.m = function(predict.data, t){
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
        if(!is.U)
          return(tcrossprod(predict.data, coef.m.sim))
        if(is.U)
          return(tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)
      }
      
      predict.y = function(predict.data, t, m) {
        predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
        j.mediator = which(colnames(predict.data) == gsub("[^[:alnum:]\\:\\s]", "", mediator))
        predict.results = matrix(NA, nrow(predict.data), nmc)
        for(i in 1:nmc){
          denominator = sqrt(sd.m^2 * (sum(coef.y.sim[i, colnames(coef.y.sim)[grepl(mediator, colnames(coef.y.sim))]]))^2 + 1)
          predict.data[, j.mediator] = m[, i]
          if(!is.U)
            predict.results[, i] = pnorm(tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ]))/denominator)
          if(is.U)
            predict.results[, i] = pnorm((tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ])) + b.y * data$Unmeasure)/denominator)
        }
        return(predict.results)
      }
      
      y1m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y1m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
      }
      y0m1 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
      }
      y0m0 = function(predict.m.data, predict.y.data){
        return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
      }
    }
    
    # Step 4. Calculate the final effects for each sample
    est.y1m1 = y1m1(predict.m.data, predict.y.data)
    est.y1m0 = y1m0(predict.m.data, predict.y.data)
    est.y0m0 = y0m0(predict.m.data, predict.y.data)
    est.y0m1 = y0m1(predict.m.data, predict.y.data)
    
    if(is.null(plot.effect)){
      if(is.U){
        sens.results = NULL
        if("TIE" %in% sens.effect){
          TIE = apply(est.y1m1 - est.y1m0, 2, mean)
          sens.results = cbind(sens.results, TIE)
        }
        if("PDE" %in% sens.effect){
          PDE = apply(est.y1m0 - est.y0m0, 2, mean)
          sens.results = cbind(sens.results, PDE)
        }
        if("PIE" %in% sens.effect){
          PIE = apply(est.y0m1 - est.y0m0, 2, mean)
          sens.results = cbind(sens.results, PIE)
        }
        if("TDE" %in% sens.effect){
          TDE = apply(est.y1m1 - est.y0m1, 2, mean)
          sens.results = cbind(sens.results, TDE)
        }
        if("INT" %in% sens.effect){
          INT = apply(est.y1m1 - est.y1m0, 2, mean) - apply(est.y0m1 - est.y0m0, 2, mean)
          sens.results = cbind(sens.results, INT)
        }
      } else {
        TIE = apply(est.y1m1 - est.y1m0, 2, mean)
        PDE = apply(est.y1m0 - est.y0m0, 2, mean)
        PIE = apply(est.y0m1 - est.y0m0, 2, mean)
        INT = TIE - PIE
        TDE = PDE + INT
        results = cbind(TIE, PIE, PDE, TDE, INT)
      }
    }
    
    if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
      if(is.null(confounders)|
         m.scale == "continuous" & y.scale == "continuous"|
         m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
        if(!is.null(ref.mod.values)){
          if(is.null(plot.effect)){
            if(is.U){
              if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              }
              if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              }
              if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              }
              if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              }
            } else {
              y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              TIE.ref = y1m1.ref - y1m0.ref
              PDE.ref = y1m0.ref - y0m0.ref
              PIE.ref = y0m1.ref - y0m0.ref
              INT.ref = TIE.ref - PIE.ref
              TDE.ref = PDE.ref + INT.ref
            }
          } else {
            if(plot.effect == "TE"){
              y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              results = y1m1.ref - y0m0.ref
            }
            if(plot.effect == "TIE"){
              y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              results = y1m1.ref - y1m0.ref
            }
            if(plot.effect == "PIE"){
              y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              results = y0m1.ref - y0m0.ref
            }
            if(plot.effect == "TDE"){
              y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              results = y1m1.ref - y0m1.ref
            }
            if(plot.effect == "PDE"){
              y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              results = y1m0.ref - y0m0.ref
            }
            if(plot.effect == "INT"){
              y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
              y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
              results = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
            }
          }
        }  
        if(!is.null(comp.mod.values)){
          if(is.null(plot.effect)){
            if(is.U){
              if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m1.dif = apply(y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref, 2, mean)
              }
              if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m0.dif = apply(y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref, 2, mean)
              }
              if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m0.dif = apply(y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref, 2, mean)
              }
              if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m1.dif = apply(y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref, 2, mean)
              }
            } else {
              y1m1.dif = apply(y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref, 2, mean)
              y1m0.dif = apply(y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref, 2, mean)
              y0m0.dif = apply(y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref, 2, mean)
              y0m1.dif = apply(y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref, 2, mean)
              TIE.dif = y1m1.dif - y1m0.dif
              PDE.dif = y1m0.dif - y0m0.dif
              PIE.dif = y0m1.dif - y0m0.dif
              INT.dif = TIE.dif - PIE.dif
              TDE.dif = PDE.dif + INT.dif
            }
          }
        }
      }
    }
  }
  
  # Bootstrap method
  if(method == 'boot'){
    est.y1m1 = est.y1m0 = est.y0m0 = est.y0m1 = y1m1.ref = y1m0.ref = y0m0.ref = y0m1.ref = y1m1.dif = y1m0.dif = y0m0.dif = y0m1.dif = NULL
    for(b in 1:nboot){
      data.boot = data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
      if(m.scale == "binary" & y.scale == "continuous"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data.boot, family = binomial(link = "probit"))
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data.boot)
        }
        if(is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot, family = binomial(link = "probit"))
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot)
        }
        
        # Step 2. model coefficient estimates
        coef.m = coef(l.m)
        coef.y = coef(l.y)
        
        # Step 3. Predict potential outcomes
        predict.m = function(predict.data, t){
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U)
            return(pnorm(tcrossprod(predict.data, t(coef.m))))
          if(is.U)
            return(pnorm((tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)))
        }
        
        predict.y = function(predict.data, t, m) {
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U)
            return(tcrossprod(predict.data, t(coef.y)))
          if(is.U)
            return(tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)
        }
        
        y1m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y1m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
        }
        y0m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y0m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
        }
      }
      
      if(m.scale == "continuous" & y.scale == "continuous"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data.boot)
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data.boot)
        }
        if(is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot)
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot)
        }
        
        # Step 2. Simulate model parameters from their sampling distribution.
        coef.m = coef(l.m)
        coef.y = coef(l.y)
        
        # Step 3. Predict potential outcomes
        predict.m = function(predict.data, t){
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U)
            return(tcrossprod(predict.data, t(coef.m)))
          if(is.U)
            return(tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)
        }
        
        predict.y = function(predict.data, t, m) {
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U)
            return(tcrossprod(predict.data, t(coef.y)))
          if(is.U)
            return(tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)
        }
        
        y1m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y1m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
        }
        y0m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y0m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
        }
      }
      
      if(m.scale == "binary" & y.scale == "binary"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data.boot, family = binomial(link = "probit"))
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data.boot, family = binomial(link = "probit"))
        }
        if(is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot, family = binomial(link = "probit"))
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot, family = binomial(link = "probit"))
        }
        
        # Step 2. Simulate model parameters from their sampling distribution.
        coef.m = coef(l.m)
        coef.y = coef(l.y)
        
        # Step 3. Predict potential outcomes
        predict.m = function(predict.data, t){
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U)
            return(pnorm(tcrossprod(predict.data, t(coef.m))))
          if(is.U)
            return(pnorm((tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)))
        }
        
        predict.y = function(predict.data, t, m) {
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U){
            return(pnorm(tcrossprod(predict.data, t(coef.y))))
          }
          if(is.U){
            return(pnorm((tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)))
          }
        }
        
        y1m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y1m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
        }
        y0m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y0m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
        }
      }
      
      if(m.scale == "continuous" & y.scale == "binary"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data.boot)
          sd.m = sigma(l.m) 
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data.boot, family = binomial(link = "probit"))
        }
        if(is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot)
          sd.m = sigma(l.m) 
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot, family = binomial(link = "probit"))
        }
        
        # Step 2. Simulate model parameters from their sampling distribution.
        coef.m = coef(l.m)
        coef.y = coef(l.y)
        
        # Step 3. Predict potential outcomes
        predict.m = function(predict.data, t){
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U)
            return(tcrossprod(predict.data, t(coef.m)))
          if(is.U)
            return(tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)
        }
        
        predict.y = function(predict.data, t, m) {
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
          predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
          predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
          if(!is.U){
            return(pnorm(tcrossprod(predict.data, t(coef.y))))
          }
          if(is.U){
            return(pnorm((tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)))
          }
        }
        
        y1m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y1m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
        }
        y0m1 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
        }
        y0m0 = function(predict.m.data, predict.y.data){
          return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
        }
      }
      
      # Step 4. Calculate the final effects for each sample
      est.y1m1 = cbind(est.y1m1, y1m1(predict.m.data, predict.y.data))
      est.y1m0 = cbind(est.y1m0, y1m0(predict.m.data, predict.y.data))
      est.y0m0 = cbind(est.y0m0, y0m0(predict.m.data, predict.y.data))
      est.y0m1 = cbind(est.y0m1, y0m1(predict.m.data, predict.y.data))
      
      if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
        if(is.null(confounders)|
           m.scale == "continuous" & y.scale == "continuous"|
           m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
          if(!is.null(ref.mod.values)){
            if(is.null(plot.effect)){
              if(is.U){
                if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                }
                if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                }
                if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                }
                if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
                }
              } else {
                y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
              }
            } else {
              if(plot.effect == "TE"){
                y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
              }
              if(plot.effect == "TIE"){
                y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
              }
              if(plot.effect == "PIE"){
                y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
                y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
              }
              if(plot.effect == "TDE"){
                y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
              }
              if(plot.effect == "PDE"){
                y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
              }
              if(plot.effect == "INT"){
                y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
              }
            }
            
          }  
          if(!is.null(comp.mod.values)){
            if(is.null(plot.effect)){
              if(is.U){
                if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y1m1.dif = cbind(y1m1.dif, y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref[, b])
                }
                if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y1m0.dif = cbind(y1m0.dif, y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref[, b])
                }
                if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y0m0.dif = cbind(y0m0.dif, y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref[, b])
                }
                if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                  y0m1.dif = cbind(y0m1.dif, y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref[, b])
                }
              } else {
                y1m1.dif = cbind(y1m1.dif, y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref[, b])
                y1m0.dif = cbind(y1m0.dif, y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref[, b])
                y0m0.dif = cbind(y0m0.dif, y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref[, b])
                y0m1.dif = cbind(y0m1.dif, y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref[, b])
              }
            }
          }
        }
      }
    }
    
    if(is.null(plot.effect)){
      if(is.U){
        sens.results = NULL
        if("TIE" %in% sens.effect){
          TIE = apply(est.y1m1 - est.y1m0, 2, mean)
          sens.results = cbind(sens.results, TIE)
        }
        if("PDE" %in% sens.effect){
          PDE = apply(est.y1m0 - est.y0m0, 2, mean)
          sens.results = cbind(sens.results, PDE)
        }
        if("PIE" %in% sens.effect){
          PIE = apply(est.y0m1 - est.y0m0, 2, mean)
          sens.results = cbind(sens.results, PIE)
        }
        if("TDE" %in% sens.effect){
          TDE = apply(est.y1m1 - est.y0m1, 2, mean)
          sens.results = cbind(sens.results, TDE)
        }
        if("INT" %in% sens.effect){
          INT = apply(est.y1m1 - est.y1m0, 2, mean) - apply(est.y0m1 - est.y0m0, 2, mean)
          sens.results = cbind(sens.results, INT)
        }
      } else {
        TIE = apply(est.y1m1 - est.y1m0, 2, mean)
        PDE = apply(est.y1m0 - est.y0m0, 2, mean)
        PIE = apply(est.y0m1 - est.y0m0, 2, mean)
        INT = TIE - PIE
        TDE = PDE + INT
        results = cbind(TIE, PIE, PDE, TDE, INT)
      }
    }
    
    if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
      if(is.null(confounders)|
         m.scale == "continuous" & y.scale == "continuous"|
         m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
        if(!is.null(ref.mod.values)){
          if(is.null(plot.effect)){
            if(is.U){
              if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m1.ref = apply(y1m1.ref, 2, mean)
              }
              if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m0.ref = apply(y1m0.ref, 2, mean)
              }
              if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m0.ref = apply(y0m0.ref, 2, mean)
              }
              if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m1.ref = apply(y0m1.ref, 2, mean)
              }
            } else {
              TIE.ref = apply(y1m1.ref - y1m0.ref, 2, mean)
              PDE.ref = apply(y1m0.ref - y0m0.ref, 2, mean)
              PIE.ref = apply(y0m1.ref - y0m0.ref, 2, mean)
              INT.ref = TIE.ref - PIE.ref
              TDE.ref = PDE.ref + INT.ref
            }
          } else {
            if(plot.effect == "TE"){
              results = apply(y1m1.ref - y0m0.ref, 2, mean)
            }
            if(plot.effect == "TIE"){
              results = apply(y1m1.ref - y1m0.ref, 2, mean)
            }
            if(plot.effect == "PIE"){
              results = apply(y0m1.ref - y0m0.ref, 2, mean)
            }
            if(plot.effect == "TDE"){
              results = apply(y1m1.ref - y0m1.ref, 2, mean)
            }
            if(plot.effect == "PDE"){
              results = apply(y1m0.ref - y0m0.ref, 2, mean)
            }
            if(plot.effect == "INT"){
              results = apply((y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref), 2, mean)
            }
          }
        }  
        if(!is.null(comp.mod.values)){
          if(is.null(plot.effect)){
            if(is.U){
              if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m1.dif = apply(y1m1.dif, 2, mean)
              }
              if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m0.dif = apply(y1m0.dif, 2, mean)
              }
              if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m0.dif = apply(y0m0.dif, 2, mean)
              }
              if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m1.dif = apply(y0m1.dif, 2, mean)
              }
            } else {
              TIE.dif = apply(y1m1.dif - y1m0.dif, 2, mean)
              PDE.dif = apply(y1m0.dif - y0m0.dif, 2, mean)
              PIE.dif = apply(y0m1.dif - y0m0.dif, 2, mean)
              INT.dif = TIE.dif - PIE.dif
              TDE.dif = PDE.dif + INT.dif
            }
          }
        }
      }
    }
  }
  
  eff.earth.fun = function(effect = NULL, plot.effect = NULL, sens.effect = NULL){
    if(!is.null(ref.mod.values)){
      if(any(is.na(ref.mod.values))){
        earth.moderators = colnames(model.matrix(as.formula(paste0("~", paste0(moderators[-which(is.na(ref.mod.values))], collapse = "+"))), data = data))[-1]
      } else {
        earth.moderators = unique(c(m.moderators.new[[1]], y.moderators.new[[1]]))
      }
      earth.moderators = gsub("[^[:alnum:]\\:\\s]", "", earth.moderators)
      earth.data = cbind(predict.m.data, predict.y.data)[, earth.moderators]
      earth.ref = cbind(predict.m.data.ref, predict.y.data.ref)[1, earth.moderators]
      ref.fun = function(x){
        smooth = earth(earth.data, x)
        ref = as.numeric(predict(smooth, earth.ref))
        return(ref)
      }
      if(is.null(plot.effect)){
        if(is.null(effect)){
          if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y1m1.ref = apply(est.y1m1, 2, ref.fun)
          }
          if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y1m0.ref = apply(est.y1m0, 2, ref.fun)
          }
          if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y0m0.ref = apply(est.y0m0, 2, ref.fun)
          }
          if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y0m1.ref = apply(est.y0m1, 2, ref.fun)
          }
          results = NULL
          if("TIE.ref" %in% sens.effect){
            TIE.ref = y1m1.ref - y1m0.ref
            results = cbind(results, TIE.ref)
          }
          if("PDE.ref" %in% sens.effect){
            PDE.ref = y1m0.ref - y0m0.ref
            results = cbind(results, PDE.ref)
          }
          if("PIE.ref" %in% sens.effect){
            PIE.ref = y0m1.ref - y0m0.ref
            results = cbind(results, PIE.ref)
          }
          if("TDE.ref" %in% sens.effect){
            TDE.ref = y1m1.ref - y0m1.ref
            results = cbind(results, TDE.ref)
          }
          if("INT.ref" %in% sens.effect){
            INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
            results = cbind(results, INT.ref)
          }
        } else {
          y1m1.ref = apply(est.y1m1, 2, ref.fun)
          y1m0.ref = apply(est.y1m0, 2, ref.fun)
          y0m0.ref = apply(est.y0m0, 2, ref.fun)
          y0m1.ref = apply(est.y0m1, 2, ref.fun)
          TIE.ref = y1m1.ref - y1m0.ref
          PDE.ref = y1m0.ref - y0m0.ref
          PIE.ref = y0m1.ref - y0m0.ref
          INT.ref = TIE.ref - PIE.ref
          TDE.ref = PDE.ref + INT.ref
          results = cbind(TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref)
        }
      } else {
        if(plot.effect == "TE"){
          y1m1.ref = apply(est.y1m1, 2, ref.fun)
          y0m0.ref = apply(est.y0m0, 2, ref.fun)
          results = y1m1.ref - y0m0.ref
        }
        if(plot.effect == "TIE"){
          y1m1.ref = apply(est.y1m1, 2, ref.fun)
          y1m0.ref = apply(est.y1m0, 2, ref.fun)
          results = y1m1.ref - y1m0.ref
        }
        if(plot.effect == "PIE"){
          y0m0.ref = apply(est.y0m0, 2, ref.fun)
          y0m1.ref = apply(est.y0m1, 2, ref.fun)
          results = y0m1.ref - y0m0.ref
        }
        if(plot.effect == "TDE"){
          y1m1.ref = apply(est.y1m1, 2, ref.fun)
          y0m1.ref = apply(est.y0m1, 2, ref.fun)
          results = y1m1.ref - y0m1.ref
        }
        if(plot.effect == "PDE"){
          y1m0.ref = apply(est.y1m0, 2, ref.fun)
          y0m0.ref = apply(est.y0m0, 2, ref.fun)
          results = y1m0.ref - y0m0.ref
        }
        if(plot.effect == "INT"){
          y1m1.ref = apply(est.y1m1, 2, ref.fun)
          y1m0.ref = apply(est.y1m0, 2, ref.fun)
          y0m0.ref = apply(est.y0m0, 2, ref.fun)
          y0m1.ref = apply(est.y0m1, 2, ref.fun)
          results = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
        }
        return(results)
      }
    }
    
    if(is.null(plot.effect)){
      if(!is.null(comp.mod.values)){
        if(any(is.na(comp.mod.values))){
          earth.moderators = colnames(model.matrix(as.formula(paste0("~", paste0(moderators[-which(is.na(comp.mod.values))], collapse = "+"))), data = data))[-1]
        } else {
          earth.moderators = unique(c(m.moderators.new[[1]], y.moderators.new[[1]]))
        }
        earth.moderators = gsub("[^[:alnum:]\\:\\s]", "", earth.moderators)
        earth.data = cbind(predict.m.data, predict.y.data)[, earth.moderators]
        earth.comp = cbind(predict.m.data.comp, predict.y.data.comp)[1, earth.moderators]
        comp.fun = function(x){
          smooth = earth(earth.data, x)
          comp = as.numeric(predict(smooth, earth.comp))
          return(comp)
        }
        if(is.null(effect)){
          if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y1m1.dif = apply(est.y1m1, 2, comp.fun) - y1m1.ref
          }
          if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y1m0.dif = apply(est.y1m0, 2, comp.fun) - y1m0.ref
          }
          if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y0m0.dif = apply(est.y0m0, 2, comp.fun) - y0m0.ref
          }
          if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
            y0m1.dif = apply(est.y0m1, 2, comp.fun) - y0m1.ref
          }
          if("TIE.dif" %in% sens.effect){
            TIE.dif = y1m1.dif - y1m0.dif
            results = cbind(results, TIE.dif)
          }
          if("PDE.dif" %in% sens.effect){
            PDE.dif = y1m0.dif - y0m0.dif
            results = cbind(results, PDE.dif)
          }
          if("PIE.dif" %in% sens.effect){
            PIE.dif = y0m1.dif - y0m0.dif
            results = cbind(results, PIE.dif)
          }
          if("TDE.dif" %in% sens.effect){
            TDE.dif = y1m1.dif - y0m1.dif
            results = cbind(results, TDE.dif)
          }
          if("INT.dif" %in% sens.effect){
            INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
            results = cbind(results, INT.dif)
          }
        } else {
          y1m1.dif = apply(est.y1m1, 2, comp.fun) - y1m1.ref
          y1m0.dif = apply(est.y1m0, 2, comp.fun) - y1m0.ref
          y0m0.dif = apply(est.y0m0, 2, comp.fun) - y0m0.ref
          y0m1.dif = apply(est.y0m1, 2, comp.fun) - y0m1.ref
          TIE.dif = y1m1.dif - y1m0.dif
          PDE.dif = y1m0.dif - y0m0.dif
          PIE.dif = y0m1.dif - y0m0.dif
          INT.dif = TIE.dif - PIE.dif
          TDE.dif = PDE.dif + INT.dif
          results = cbind(results, TIE.dif, PIE.dif, PDE.dif, TDE.dif, INT.dif)
        }
      }
    }
    if(!is.null(effect))
      return(results[, effect])
    if(!is.null(sens.effect))
      return(results)
  }
  
  if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
    if(is.null(confounders)){
      if(is.U){
        if("TIE.ref" %in% sens.effect){
          TIE.ref = y1m1.ref - y1m0.ref
          sens.results = cbind(sens.results, TIE.ref)
        }
        if("PDE.ref" %in% sens.effect){
          PDE.ref = y1m0.ref - y0m0.ref
          sens.results = cbind(sens.results, PDE.ref)
        }
        if("PIE.ref" %in% sens.effect){
          PIE.ref = y0m1.ref - y0m0.ref
          sens.results = cbind(sens.results, PIE.ref)
        }
        if("TDE.ref" %in% sens.effect){
          TDE.ref = y1m1.ref - y0m1.ref
          sens.results = cbind(sens.results, TDE.ref)
        }
        if("INT.ref" %in% sens.effect){
          INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
          sens.results = cbind(sens.results, INT.ref)
        }
        if("TIE.dif" %in% sens.effect){
          TIE.dif = y1m1.dif - y1m0.dif
          sens.results = cbind(sens.results, TIE.dif)
        }
        if("PDE.dif" %in% sens.effect){
          PDE.dif = y1m0.dif - y0m0.dif
          sens.results = cbind(sens.results, PDE.dif)
        }
        if("PIE.dif" %in% sens.effect){
          PIE.dif = y0m1.dif - y0m0.dif
          sens.results = cbind(sens.results, PIE.dif)
        }
        if("TDE.dif" %in% sens.effect){
          TDE.dif = y1m1.dif - y0m1.dif
          sens.results = cbind(sens.results, TDE.dif)
        }
        if("INT.dif" %in% sens.effect){
          INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
          sens.results = cbind(sens.results, INT.dif)
        }
      } else {
        if(!is.null(ref.mod.values)){
          if(!is.null(comp.mod.values)){
            results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref, TIE.dif, PIE.dif, PDE.dif, TDE.dif, INT.dif)
          } else {
            if(is.null(plot.effect)){
              results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref)
            } else {
              return(results)
            }
          }
        }
      }
    } else if(m.scale == "continuous" & y.scale == "continuous"){
      if(!paste0(treatment, ":", mediator) %in% y.predictors){
        if(is.U){
          if("TIE.ref" %in% sens.effect){
            TIE.ref = y1m1.ref - y1m0.ref
            sens.results = cbind(sens.results, TIE.ref)
          }
          if("PDE.ref" %in% sens.effect){
            PDE.ref = y1m0.ref - y0m0.ref
            sens.results = cbind(sens.results, PDE.ref)
          }
          if("PIE.ref" %in% sens.effect){
            PIE.ref = y0m1.ref - y0m0.ref
            sens.results = cbind(sens.results, PIE.ref)
          }
          if("TDE.ref" %in% sens.effect){
            TDE.ref = y1m1.ref - y0m1.ref
            sens.results = cbind(sens.results, TDE.ref)
          }
          if("INT.ref" %in% sens.effect){
            INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
            sens.results = cbind(sens.results, INT.ref)
          }
          if("TIE.dif" %in% sens.effect){
            TIE.dif = y1m1.dif - y1m0.dif
            sens.results = cbind(sens.results, TIE.dif)
          }
          if("PDE.dif" %in% sens.effect){
            PDE.dif = y1m0.dif - y0m0.dif
            sens.results = cbind(sens.results, PDE.dif)
          }
          if("PIE.dif" %in% sens.effect){
            PIE.dif = y0m1.dif - y0m0.dif
            sens.results = cbind(sens.results, PIE.dif)
          }
          if("TDE.dif" %in% sens.effect){
            TDE.dif = y1m1.dif - y0m1.dif
            sens.results = cbind(sens.results, TDE.dif)
          }
          if("INT.dif" %in% sens.effect){
            INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
            sens.results = cbind(sens.results, INT.dif)
          }
        } else {
          if(!is.null(ref.mod.values)){
            if(!is.null(comp.mod.values)){
              results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref, TIE.dif, PIE.dif, PDE.dif, TDE.dif, INT.dif)
            } else {
              if(is.null(plot.effect)){
                results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref)
              } else {
                return(results)
              }
            }
          }
        }
      } else {
        if(is.U){
          if("TIE.ref" %in% sens.effect){
            TIE.ref = y1m1.ref - y1m0.ref
            sens.results = cbind(sens.results, TIE.ref)
          }
          if("PIE.ref" %in% sens.effect){
            PIE.ref = y0m1.ref - y0m0.ref
            sens.results = cbind(sens.results, PIE.ref)
          }
          if("INT.ref" %in% sens.effect){
            INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
            sens.results = cbind(sens.results, INT.ref)
          }
          if("TIE.dif" %in% sens.effect){
            TIE.dif = y1m1.dif - y1m0.dif
            sens.results = cbind(sens.results, TIE.dif)
          }
          if("PIE.dif" %in% sens.effect){
            PIE.dif = y0m1.dif - y0m0.dif
            sens.results = cbind(sens.results, PIE.dif)
          }
          if("INT.dif" %in% sens.effect){
            INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
            sens.results = cbind(sens.results, INT.dif)
          }
          if(any(sens.effect %in% c("TDE.ref", "PDE.ref", "TDE.dif", "PDE.dif"))){
            if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TIE.ref", "PIE.ref", "INT.ref", "TIE.dif", "PIE.dif", "INT.dif"))){
              sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TIE.ref", "PIE.ref", "INT.ref", "TIE.dif", "PIE.dif", "INT.dif"))]))
            } else {
              sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
            }
          }
        } else {
          if(!is.null(ref.mod.values)){
            if(!is.null(comp.mod.values)){
              results = cbind(results, TIE.ref, PIE.ref, INT.ref, TIE.dif, PIE.dif, INT.dif, eff.earth.fun(effect = c("PDE.ref", "TDE.ref", "PDE.dif", "TDE.dif")))
            } else {
              if(is.null(plot.effect)){
                results = cbind(results, TIE.ref, PIE.ref, INT.ref, eff.earth.fun(effect = c("PDE.ref", "TDE.ref")))
              } else {
                if(plot.effect == "TIE"|plot.effect == "PIE"|plot.effect == "INT"){
                  return(results)
                } else {
                  return(eff.earth.fun(plot.effect = plot.effect))
                }
              }
            }
          }
        }
      }
    } else if(m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
      if(is.U){
        if("PDE.ref" %in% sens.effect){
          PDE.ref = y1m0.ref - y0m0.ref
          sens.results = cbind(sens.results, PDE.ref)
        }
        if("TDE.ref" %in% sens.effect){
          TDE.ref = y1m1.ref - y0m1.ref
          sens.results = cbind(sens.results, TDE.ref)
        }
        if("INT.ref" %in% sens.effect){
          INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
          sens.results = cbind(sens.results, INT.ref)
        }
        if("PDE.dif" %in% sens.effect){
          PDE.dif = y1m0.dif - y0m0.dif
          sens.results = cbind(sens.results, PDE.dif)
        }
        if("TDE.dif" %in% sens.effect){
          TDE.dif = y1m1.dif - y0m1.dif
          sens.results = cbind(sens.results, TDE.dif)
        }
        if("INT.dif" %in% sens.effect){
          INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
          sens.results = cbind(sens.results, INT.dif)
        }
        if(any(sens.effect %in% c("TIE.ref", "PIE.ref", "TIE.dif", "PIE.dif"))){
          if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TDE.ref", "PDE.ref", "INT.ref", "TDE.dif", "PDE.dif", "INT.dif"))){
            sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TDE.ref", "PDE.ref", "INT.ref", "TDE.dif", "PDE.dif", "INT.dif"))]))
          } else {
            sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
          }
        }
      } else {
        if(!is.null(ref.mod.values)){
          if(!is.null(comp.mod.values)){
            results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "TIE.dif", "PIE.dif")), PDE.ref, TDE.ref, INT.ref, PDE.dif, TDE.dif, INT.dif)
          } else {
            if(is.null(plot.effect)){
              results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref")), PDE.ref, TDE.ref, INT.ref)
            } else {
              if(plot.effect == "TDE"|plot.effect == "PDE"|plot.effect == "INT"){
                return(results)
              } else {
                return(eff.earth.fun(plot.effect = plot.effect))
              }
            }
          }
        }
      }
    } else {
      if(is.U){
        if(any(sens.effect %in% c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif"))){
          if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))){
            sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))]))
          } else {
            sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
          }
        }
      } else {
        if(!is.null(ref.mod.values)){
          if(!is.null(comp.mod.values)){
            results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif")))
          } else {
            if(is.null(plot.effect)){
              results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref")))
            } else {
              return(eff.earth.fun(plot.effect = plot.effect))
            }
          }
        }
      }
    }
  } else {
    if(is.U){
      if(any(sens.effect %in% c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif"))){
        if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))){
          sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))]))
        } else {
          sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
        }
      }
    } else {
      if(!is.null(ref.mod.values)){
        if(!is.null(comp.mod.values)){
          results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif")))
        } else {
          if(is.null(plot.effect)){
            results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref")))
          } else {
            return(eff.earth.fun(plot.effect = plot.effect))
          }
        }
      }
    }
  }
  
  if(is.null(plot.effect)){
    if(is.U){
      return(sens.results[, sens.effect])
    } else {
      # Step 5. Compute summary statistics such as point estimates and confidence intervals.
      results.ori = results
      TE = results.ori[, "TIE"] + results.ori[, "PDE"]
      results = cbind(TE, results.ori[, 1:5])
      if(!is.null(ref.mod.values)){
        TE.ref = results.ori[, "TIE.ref"] + results.ori[, "PDE.ref"]
        results = cbind(results, TE.ref, results.ori[, 6:10])
        if(!is.null(comp.mod.values)){
          TE.dif = results.ori[, "TIE.dif"] + results.ori[, "PDE.dif"]
          results = cbind(results, TE.dif, results.ori[, 11:15])
        }
      }
      
      est.results = apply(results, 2, mean)
      se.results = apply(results, 2, sd)
      ci.results = apply(results, 2, quantile, probs = c((1 - conf.level)/2, (1 + conf.level)/2))
      
      
      # Output
      # Mediator and outcome models
      summary.table = function(moderators, predictors, l){
        if(is.null(moderators)){
          return(summary(l))
        }
        
        if(!is.null(moderators)){
          model.coef = names(coef(l))
          
          if(all(predictors %in% y.predictors.new) & all(y.predictors.new %in% predictors)){
            if(paste0(treatment, ":", mediator) %in% model.coef)
              predictors[which(predictors == paste0(treatment, ":", mediator))] = paste0(treatment, ":", mediator)
            if(paste0(mediator, ":", treatment) %in% model.coef)
              predictors[which(predictors == paste0(treatment, ":", mediator))] = paste0(mediator, ":", treatment)
          }
          
          main = rep("Intercept", 1 + length(moderators[[1]]))
          moderation = c("Intercept", moderators[[1]])
          for(i in 1:length(predictors)){
            moderation = c(moderation, "Intercept", moderators[[i + 1]])
            main = c(main, rep(predictors[i], 1 + length(moderators[[i + 1]])))
          }
          
          coef = NULL
          for(i in 1:length(main)){
            if(main[i] == "Intercept"){
              coef = c(coef, paste0(moderation[i]))
            } else {
              if(moderation[i] == "Intercept"){
                coef = c(coef, paste0(main[i]))
              } else {
                main.i = unlist(strsplit(main[i], ":"))
                if(length(main.i) == 1){
                  if(paste0(main[i], ":", moderation[i]) %in% model.coef)
                    coef = c(coef, paste0(main[i], ":", moderation[i]))
                  if(paste0(moderation[i], ":", main[i]) %in% model.coef)
                    coef = c(coef, paste0(moderation[i], ":", main[i]))
                }
                if(length(main.i) == 2){
                  if(paste0(main.i[1], ":", main.i[2], ":", moderation[i]) %in% model.coef)
                    coef = c(coef, paste0(main.i[1], ":", main.i[2], ":", moderation[i]))
                  if(paste0(main.i[1], ":", moderation[i], ":", main.i[2]) %in% model.coef)
                    coef = c(coef, paste0(main.i[1], ":", moderation[i], ":", main.i[2]))
                  if(paste0(main.i[2], ":", main.i[1], ":", moderation[i]) %in% model.coef)
                    coef = c(coef, paste0(main.i[2], ":", main.i[1], ":", moderation[i]))
                  if(paste0(main.i[2], ":", moderation[i], ":", main.i[1]) %in% model.coef)
                    coef = c(coef, paste0(main.i[2], ":", moderation[i], ":", main.i[1]))
                  if(paste0(moderation[i], ":", main.i[1], ":", main.i[2]) %in% model.coef)
                    coef = c(coef, paste0(moderation[i], ":", main.i[1], ":", main.i[2]))
                  if(paste0(moderation[i], ":", main.i[2], ":", main.i[1]) %in% model.coef)
                    coef = c(coef, paste0(moderation[i], ":", main.i[2], ":", main.i[1]))
                }
              }
            }
          }
          
          main = c("Intercept", rep("", length(moderators[[1]])))
          for(i in 1:length(predictors)){
            main = c(main, predictors[i], rep("", length(moderators[[i + 1]])))
          }
          
          coef[which(coef == "Intercept")] = "(Intercept)"
          summary.table = as.data.frame(matrix(NA, length(main), 6))
          summary.table[, 1] = main
          summary.table[, 2] = moderation
          summary.table[, 3:6] = as.data.frame(summary(l)$coefficients[coef, ])
          mod_summary_sign = summary.table[, 6]  
          mod_summary_stars = NA                             # Named vector with significance stars
          mod_summary_stars[mod_summary_sign < 0.1] = "."
          mod_summary_stars[mod_summary_sign < 0.05] = "*"
          mod_summary_stars[mod_summary_sign < 0.01] = "**"
          mod_summary_stars[mod_summary_sign < 0.001] = "***"
          mod_summary_stars[is.na(mod_summary_stars)] = ""
          summary.table = cbind(summary.table, mod_summary_stars)
          colnames(summary.table) = c("main", "moderation", colnames(summary(l)$coefficients), "")
          rownames(summary.table) = coef
          return(summary.table)
        }
      }
      
      summary.m = summary.table(m.moderators.new, m.predictors.new, l.m)
      summary.y = summary.table(y.moderators.new, y.predictors.new, l.y)
      
      # Causal effects
      summary.effects = cbind(est.results, se.results, ci.results[1, ], ci.results[2, ])
        
      colnames(summary.effects) = c("Estimate", "Std. Error", paste(conf.level * 100, "% CI Lower ", (1-conf.level)/2 *100, "%", sep=""),
                                    paste(conf.level * 100, "% CI Upper ", (1-conf.level)/2 *100, "%", sep=""))
      
      args = as.list(match.call())[2:length(as.list(match.call()))]
      args.full = mget(names(formals()), sys.frame(sys.nframe()))
      
      if(y.scale == "continuous")
        res = list(effects = summary.effects, m.model = summary.m, y.model = summary.y, results = results, args = args, args.full = args.full, l.m = l.m, l.y = l.y, formula.m = formula.m, formula.y = formula.y, m.predictors = m.predictors, y.predictors = y.predictors, m.moderators = m.moderators, y.moderators = y.moderators, m.predictors.new = m.predictors.new, y.predictors.new = y.predictors.new, m.moderators.new = m.moderators.new, y.moderators.new = y.moderators.new, predict.m.data.ref = predict.m.data.ref, predict.m.data.comp = predict.m.data.comp, predict.m.data = predict.m.data, predict.y.data.ref = predict.y.data.ref, predict.y.data.comp = predict.y.data.comp, predict.y.data = predict.y.data, data = data)
      if(y.scale == "binary")
        res = list(effects = summary.effects, m.model = summary.m, y.model = summary.y, results = results, args = args, args.full = args.full, l.m = l.m, l.y = l.y, formula.m = formula.m, formula.y = formula.y, m.predictors = m.predictors, y.predictors = y.predictors, m.moderators = m.moderators, y.moderators = y.moderators, m.predictors.new = m.predictors.new, y.predictors.new = y.predictors.new, m.moderators.new = m.moderators.new, y.moderators.new = y.moderators.new, predict.m.data.ref = predict.m.data.ref, predict.m.data.comp = predict.m.data.comp, predict.m.data = predict.m.data, predict.y.data.ref = predict.y.data.ref, predict.y.data.comp = predict.y.data.comp, predict.y.data = predict.y.data, data = data)

      # giving the output a class to allow users to using summary() to get the results using the s3 method
      class(res) <- "addition"
      return(res)
    }
  }
}

#' Summarizing Output for Causal Moderated Mediation Analysis
#' 
#' 'summary_modmed' is used to report from causal moderated mediation analysis. Alternatively, 'summary' can also be used to generate the same output.
#' @param object 	output from \code{modmed} function
#' @param ... additional arguments affecting the summary produced.
#' @return \code{modmed} returns causal moderated mediation analysis results. The \code{summary_modmed} function provides summary tables of the results.
#' @author Xu Qin and Lijuan Wang
#' @references Qin, X., & Wang, L. (2023). Causal moderated mediation analysis: Methods and software
#' @importFrom stats df.residual pf anova AIC getCall
#' @examples
#' \donttest{
#' data(newws)
#' modmed.results = modmed(data = newws, treatment = "treat", mediator = "emp", outcome = "depression", covariates.disc = c("emp_prior", "nevmar", "hispanic", "nohsdip"), covariates.cont = c("workpref", "attitude", "depress_prior"), moderators.disc = "CHCNT", moderators.cont = "ADCPC", m.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), y.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), mediator = c("ADCPC", "CHCNT"), tm = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), comp.mod.disc.values = 3, ref.mod.disc.values = 2, comp.mod.cont.values = 5050, ref.mod.cont.values = 5050, m.scale = "binary", y.scale = "continuous", seed = 1) 
#' summary_modmed(modmed.results)
#' summary(modmed.results)
#' }

#' @export 
summary_modmed = function(object, ...){
  args = object$args
  args.full = object$args.full 
  m.scale = args$m.scale
  y.scale = args$y.scale
  method = args.full$method
  summary.effects = object$effects
  summary.m = object$m.model
  summary.y = object$y.model
  
  # Print causal effect estimation and inference results
  cat("Treatment:", args.full$treatment, "\n\n")
  cat("Mediator:", args.full$mediator, "\n\n")
  cat("Outcome:", args.full$outcome, "\n\n")
  cat("Pre-treatment confounders:", paste0(c(args.full$covariates.disc, args.full$covariates.cont), collapse = ", "), "\n\n")
  cat("Moderators:", paste0(c(args.full$moderators.disc, args.full$moderators.cont), collapse = ", "), "\n\n")
  cat("Compare values of the treatment:", args.full$comp.treatment.value, "\n\n")
  cat("Reference values of the treatment:", args.full$ref.treatment.value, "\n\n")
  cat("Compare values of the moderators:", paste0(c(args.full$comp.mod.disc.values, args.full$comp.mod.cont.values), collapse = ", "), "\n\n")
  cat("Reference values of the moderators:", paste0(c(args.full$ref.mod.disc.values, args.full$ref.mod.cont.values), collapse = ", "), "\n\n")
  if(method == 'mc')
    cat("Estimation method: Monte Carlo Method \n\n")
  if(method == 'boot')
    cat("Estimation method: Bootstrap Method \n\n")
  cat("Causal Effects:\n")
  print(summary.effects)
  cat("---")
  nmc = args.full$nmc
  cat("\nCI is confidence interval constructed based on simulation of mediator and outcome model parameters (number of simulations is ", nmc, ")", sep = "")
  
  l.m = object$l.m
  l.y = object$l.y
  cat("\n\nMediator Model:\n")
  print(getCall(l.m))
  print(summary.m)
  cat("---")
  cat("\nSignif.codes:")
  cat("\n0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  if(m.scale == "continuous"){
    cat("Residual standard error:", round(summary(l.m)[[6]], 2), "on", df.residual(l.m), "degrees of freedom")
    cat("\nMultiple R-squared:  ", round(summary(l.m)$r.squared, 4),"Adjusted R-squared: ", round(summary(l.m)$adj.r.squared, 5))
    cat("\nF-statistic:", round(summary(l.m)$fstatistic[1], 3), "on", summary(l.m)$fstatistic[2], "and", summary(l.m)$fstatistic[3], "DF,  p-value:", round(pf(summary(l.m)$fstatistic[1], summary(l.m)$fstatistic[2], summary(l.m)$fstatistic[3],lower.tail = FALSE), 3))
  }
  if(m.scale == "binary"){
    cat("Null deviance:", round(anova(l.m)$"Resid. Dev"[1], 2), "on", anova(l.m)$"Resid. Df"[1], "degrees of freedom")
    cat("\nResidual deviance:", round(anova(l.m)$"Resid. Dev"[length(anova(l.m)$"Resid. Dev")], 2), "on", anova(l.m)$"Resid. Df"[length(anova(l.m)$"Resid. Df")], "degrees of freedom")
    cat("\nAIC:", AIC(l.m))
  }
  cat("\n\nOutcome Model:\n")
  print(getCall(l.y))
  print(summary.y)
  cat("---")
  cat("\nSignif.codes:")
  cat("\n0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  if(y.scale == "continuous"){
    cat("Residual standard error:", round(summary(l.y)[[6]], 2), "on", df.residual(l.y), "degrees of freedom")
    cat("\nMultiple R-squared:  ", round(summary(l.y)$r.squared, 4),"Adjusted R-squared: ", round(summary(l.y)$adj.r.squared, 5))
    cat("\nF-statistic:", round(summary(l.y)$fstatistic[1], 3), "on", summary(l.y)$fstatistic[2], "and", summary(l.y)$fstatistic[3], "DF,  p-value:", round(pf(summary(l.y)$fstatistic[1], summary(l.y)$fstatistic[2], summary(l.y)$fstatistic[3],lower.tail = FALSE), 3))
  }
  if(y.scale == "binary"){
    cat("Null deviance:", round(anova(l.y)$"Resid. Dev"[1], 2), "on", anova(l.y)$"Resid. Df"[1], "degrees of freedom")
    cat("\nResidual deviance:", round(anova(l.y)$"Resid. Dev"[length(anova(l.y)$"Resid. Dev")], 2), "on", anova(l.y)$"Resid. Df"[length(anova(l.y)$"Resid. Df")], "degrees of freedom")
    cat("\nAIC:", AIC(l.y))
  }
  cat("\n\n") 
}

#' @export 
summary.addition = summary_modmed

#' Visual Representation of the Causal Moderated Mediation Analysis Results
#' 
#' 'modmed.plot' is used to visualize results from \code{modmed} function. This applies only if moderators.disc or moderators.cont is not NULL. The plot consists of two parts. The top represents the sampling distribution of the specified causal effect as a function of the specified moderator within the given levels of the other moderators. The bottom represents the distribution of the specified moderator on the x axis. 
#' @param object Output from the \code{modmed} function.
#' @param effect A character string indicating which causal effect to be plotted. effect can be specified as "TE", "TIE", "PIE", "PDE", "TDE", "INT", "TE.ref", "TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TE.dif", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", or "INT.dif".
#' @param moderator A character string indicating which moderator to be plotted. It must be one of the moderators specified in the function of \code{modmed}.
#' @param other.mod.disc.values A vector of values of the other discrete moderators given which the conditional effect at each value of the specified moderator is estimated. The order of other.mod.disc.values should be the same as moderators.disc specified in the function of \code{modmed}, with the specified moderator removed if the specified moderator is discrete. If one does not want to condition some moderators on specific values, one may specify their values to be NA. NULL if there are no other discrete moderators.
#' @param other.mod.cont.values A vector of values of the other continuous moderators given which the conditional effect at each value of the specified moderator is estimated. The order of other.mod.cont.values should be the same as moderators.cont specified in the function of \code{modmed}, with the specified moderator removedif the specified moderator is continuous. If one does not want to condition some moderators on specific values, one may specify their values to be NA. NULL if there are no other continuous moderators.
#' @param is.dist.moderator A logical value. "TRUE" if distribution of the specified moderator is to be plotted at the bottom, and "FALSE" if not. The default is "FALSE".
#' @param probs A vector of percentiles to be plotted on the distribution of the moderator if the moderator is continuous. NULL if is.dist.moderator = FALSE. The default is c(0.1, 0.25, 0.5, 0.75, 0.9).
#' @param ncore The number of cores for parallel computing. The default is the number of CPU cores on the current host minus 1. One core is saved for users to run other programs on the computer while running the R function. 
#' @return \code{modmed} returns causal moderated mediation analysis results. The \code{plot.modmed} function plots the results.
#' @author Xu Qin and Lijuan Wang
#' @references Qin, X., & Wang, L. (2023). Causal moderated mediation analysis: Methods and software
#' @export
#' @importFrom ggplot2 ggplot aes labs geom_boxplot geom_errorbar stat_summary geom_smooth geom_ribbon scale_fill_manual scale_x_continuous geom_vline geom_point geom_text geom_hline theme element_text geom_line scale_fill_brewer
#' @importFrom cowplot plot_grid
#' @importFrom scales pretty_breaks
#' @importFrom stats density loess
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel detectCores makeCluster stopCluster 
#' @importFrom foreach foreach %dopar%
#' @examples
#' \donttest{
#' data(newws)
#' modmed.results = modmed(data = newws, treatment = "treat", mediator = "emp", outcome = "depression", covariates.disc = c("emp_prior", "nevmar", "hispanic", "nohsdip"), covariates.cont = c("workpref", "attitude", "depress_prior"), moderators.disc = "CHCNT", moderators.cont = "ADCPC", m.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), y.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), mediator = c("ADCPC", "CHCNT"), tm = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), comp.mod.disc.values = 3, ref.mod.disc.values = 2, comp.mod.cont.values = 5050, ref.mod.cont.values = 5050, m.scale = "binary", y.scale = "continuous", seed = 1) 
#' modmed.plot(modmed.results, effect = "TIE", moderator = "ADCPC", other.mod.disc.values = 1, is.dist.moderator = TRUE, ncore = 1)
#' }
modmed.plot = function(object, effect, moderator, other.mod.disc.values = NULL, other.mod.cont.values = NULL, is.dist.moderator = FALSE, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), ncore = 2){
  dens.x = dens.y = qquant = NULL # Otherwise no visible binding for global variable
  args.full = object$args.full
  if(!moderator %in% args.full$moderators.disc & !moderator %in% args.full$moderators.cont)
    stop("moderator must be specified as one of the moderators in moderator.disc or moderator.cont specified in the modmed function")
  if(moderator %in% args.full$moderators.disc)
    moderator.scale = "discrete"
  if(moderator %in% args.full$moderators.cont)
    moderator.scale = "continuous"
  data = object$data
  moderators.disc = args.full$moderators.disc
  moderators.cont = args.full$moderators.cont
  conf.level = args.full$conf.level
  y.scale = args.full$y.scale
  
  if(!moderator %in% moderators.disc & !moderator %in% moderators.cont)
    stop("moderator must be one of the moderators specified in the function of modmed")
  
  if(moderator.scale == "discrete"){
    if(length(moderators.disc) != length(other.mod.disc.values) + 1)
      stop("If the specified moderator is discrete, other.mod.disc.values should be in the same order as the discrete moderators specified in the function of modmed, with the specified moderator removed")
    if(any(length(moderators.cont) != length(other.mod.cont.values)))
      stop("If the specified moderator is discrete, other.mod.cont.values should be in the same order as the continuous moderators specified in the function of modmed")
  } 
  if(moderator.scale == "continuous"){
    if(length(moderators.cont) != length(other.mod.cont.values) + 1)
      stop("If the specified moderator is continuous, other.mod.cont.values should be the values for the continuous moderators specified in the function of modmed, with the specified moderator removed")
    if(any(length(moderators.disc) != length(other.mod.disc.values)))
      stop("If the specified moderator is continuous, other.mod.disc.values should be in the same order as the discrete moderators specified in the function of modmed")
  }
  if(moderator.scale == "discrete" & is.dist.moderator == TRUE){
    stop("If the specified moderator is discrete, its sample distribution will not be plotted.")
  }
  
  if(moderator.scale == "discrete")
    unique.moderator.values = as.character(unique(data[, moderator]))
  if(moderator.scale == "continuous")
    unique.moderator.values = seq(min(data[, moderator]), max(data[, moderator]), length.out = 15)
  n.unique.moderator.values = length(unique.moderator.values)
  
  cl = makeCluster(ncore)
  on.exit(stopCluster(cl))
  registerDoSNOW(cl)
  
  if(interactive()){
    pb = txtProgressBar(min = 0,      # Minimum value of the progress bar
                        max = n.unique.moderator.values, # Maximum value of the progress bar
                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                        # width = 50,   # Progress bar width. Defaults to getOption("width")
                        char = "=")   # Character used to create the bar
  }
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  res = foreach(i = 1:n.unique.moderator.values, .options.snow = opts) %dopar% {
    args.full$comp.mod.disc.values = NULL
    args.full$comp.mod.cont.values = NULL
    
    upd.moderators.values = function(moderators.disc.values, moderators.disc, moderators.cont.values, moderators.cont, new.values){
      if(moderator.scale == "discrete"){
        moderators.disc.values[which(moderators.disc == moderator)] = new.values
        if(!is.null(other.mod.disc.values)){
          if(sum(is.na(other.mod.disc.values)) > 0){
            moderators.disc.values[which(moderators.disc != moderator)][-which(is.na(other.mod.disc.values))] = other.mod.disc.values[-which(is.na(other.mod.disc.values))]
          } else {
            moderators.disc.values[which(moderators.disc != moderator)] = other.mod.disc.values
          }
        }
        if(!is.null(other.mod.cont.values)){
          if(sum(is.na(other.mod.cont.values)) > 0){
            moderators.cont.values[-which(is.na(other.mod.cont.values))] = other.mod.cont.values[-which(is.na(other.mod.cont.values))]
          } else {
            moderators.cont.values = other.mod.cont.values
          }
        }
      }
      if(moderator.scale == "continuous"){
        moderators.cont.values[which(moderators.cont == moderator)] = new.values
        if(!is.null(other.mod.cont.values)){
          if(sum(is.na(other.mod.cont.values)) > 0){
            moderators.cont.values[which(moderators.cont != moderator)][-which(is.na(other.mod.cont.values))] = other.mod.cont.values[-which(is.na(other.mod.cont.values))]
          } else {
            moderators.cont.values[which(moderators.cont != moderator)] = other.mod.cont.values
          }
        }
        if(!is.null(other.mod.disc.values)){
          if(sum(is.na(other.mod.disc.values)) > 0){
            moderators.disc.values[-which(is.na(other.mod.disc.values))] = other.mod.disc.values[-which(is.na(other.mod.disc.values))]
          } else {
            moderators.disc.values = other.mod.disc.values
          }
        }
      }
      return(list(moderators.disc.values = moderators.disc.values, moderators.cont.values = moderators.cont.values))
    }
    
    results = do.call(upd.moderators.values, list(moderators.disc.values = args.full$ref.mod.disc.values, moderators.disc = args.full$moderators.disc, moderators.cont.values = args.full$ref.mod.cont.values, moderators.cont = args.full$moderators.cont, new.values = unique.moderator.values[i]))
    args.full$ref.mod.disc.values = results$moderators.disc.values
    args.full$ref.mod.cont.values = results$moderators.cont.values
    args.full$plot.effect = effect
    est.effect = do.call(modmed, args.full)
    moderator.values = rep(unique.moderator.values[i], length(est.effect))
    CIL = rep(quantile(est.effect, probs = (1 - conf.level)/2), length(est.effect))
    CIU = rep(quantile(est.effect, probs = (1 + conf.level)/2), length(est.effect))
    return(cbind.data.frame(est.effect = est.effect, moderator.values = moderator.values, CIL = CIL, CIU = CIU))
  }
  
  if(interactive()){
    close(pb)
  }
  
  est.effect = NULL
  moderator.values = NULL
  CIL = NULL
  CIU = NULL
  for(i in 1:n.unique.moderator.values){
    est.effect = c(est.effect, res[[i]]$est.effect)
    moderator.values = c(moderator.values, res[[i]]$moderator.values)
    CIL = c(CIL, res[[i]]$CIL)
    CIU = c(CIU, res[[i]]$CIU)
  }
  
  if(moderator.scale == "discrete")
    results = cbind.data.frame(factor(moderator.values, levels = levels(data[, moderator]), ordered = TRUE), est.effect, CIL, CIU)
  if(moderator.scale == "continuous")
    results = cbind.data.frame(as.factor(moderator.values), est.effect, CIL, CIU)
  colnames(results) = c(moderator, effect, "CIL", "CIU")
  
  ylab = paste("Conditional", effect)
  if(moderator.scale == "discrete"){
    pMain = ggplot(results, aes(x = results[, moderator], y = results[, effect])) +
      labs(x = moderator, y = ylab) +
      geom_boxplot() + 
      geom_errorbar(aes(ymin = results[, "CIL"], ymax = results[, "CIU"]), width=.2, colour = "blue") + 
      stat_summary(fun = mean, geom = "point", colour = "blue")
  }
  if(moderator.scale == "continuous"){
    data.loess.CIU = cbind.data.frame(mod = as.numeric(levels(results[, moderator]))[results[, moderator]], loess.y = results[, "CIU"])
    loess.CIU = loess(loess.y ~ mod, data = data.loess.CIU)
    predict.CIU = predict(loess.CIU, newdata = data.frame(mod = seq(min(data[, moderator]), max(data[, moderator]), length.out = 1000)))
    data.loess.CIL = cbind.data.frame(mod = as.numeric(levels(results[, moderator]))[results[, moderator]], loess.y = results[, "CIL"])
    loess.CIL = loess(loess.y ~ mod, data = data.loess.CIL)
    predict.CIL = predict(loess.CIL, newdata = data.frame(mod = seq(min(data[, moderator]), max(data[, moderator]), length.out = 1000)))
    if(any(is.na(predict.CIU)))
      predict.CIU = predict.CIU[-which(is.na(predict.CIU))]
    if(any(is.na(predict.CIL)))
      predict.CIL = predict.CIL[-which(is.na(predict.CIL))]
    CIU.0 = NULL
    for(i in 1:(length(predict.CIU) - 1)){
      if(predict.CIU[i] <= 0 & predict.CIU[i + 1] >= 0|predict.CIU[i] >= 0 & predict.CIU[i + 1] <= 0)
        CIU.0 = c(CIU.0, i)
    }
    CIL.0 = NULL
    for(i in 1:(length(predict.CIL) - 1)){
      if(predict.CIL[i] <= 0 & predict.CIL[i + 1] >= 0|predict.CIL[i] >= 0 & predict.CIL[i + 1] <= 0)
        CIL.0 = c(CIL.0, i)
    }
    xintercept = NULL
    if(!is.null(CIU.0)){
      xintercept = c(xintercept, seq(min(data[, moderator]), max(data[, moderator]), length.out = 1000)[CIU.0])
    }
    if(!is.null(CIL.0))
      xintercept = c(xintercept, seq(min(data[, moderator]), max(data[, moderator]), length.out = 1000)[CIL.0])
    
    if(is.null(xintercept)){
      pMain = ggplot(results, aes(x = results[, moderator], y = results[, effect])) +
        labs(x = moderator, y = ylab) +
        geom_smooth(aes(x = as.numeric(levels(results[, moderator]))[results[, moderator]], y = results[, effect]), formula = y ~ x, method = "loess", se = FALSE) +
        geom_ribbon(aes(ymin = predict(loess.CIL), ymax = predict(loess.CIU), x = as.numeric(levels(results[, moderator]))[results[, moderator]], fill = "band"), alpha = 0.5, show.legend = FALSE) +
        scale_fill_manual("",values="grey") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(min(as.numeric(levels(results[, moderator]))[results[, moderator]]), max(as.numeric(levels(results[, moderator]))[results[, moderator]]))) + 
        geom_vline(xintercept = xintercept, linetype = "dashed", col = "grey48")
    }
    if(length(xintercept) == 1){
      pMain = ggplot(results, aes(x = results[, moderator], y = results[, effect])) +
        labs(x = moderator, y = ylab) +
        geom_smooth(aes(x = as.numeric(levels(results[, moderator]))[results[, moderator]], y = results[, effect]), formula = y ~ x, method = "loess", se = FALSE) +
        geom_ribbon(aes(ymin = predict(loess.CIL), ymax = predict(loess.CIU), x = as.numeric(levels(results[, moderator]))[results[, moderator]], fill = "band"), alpha = 0.5, show.legend = FALSE) +
        scale_fill_manual("",values="grey") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(min(as.numeric(levels(results[, moderator]))[results[, moderator]]), max(as.numeric(levels(results[, moderator]))[results[, moderator]]))) + 
        geom_point(aes(x = xintercept, y = 0), col = "grey18") +
        geom_text(aes(x = xintercept, y = 0, label = paste0("(", round(xintercept, 2), ", 0)")), size = 3,  col = "grey18", hjust = 0.5, vjust = 1.5) +
        geom_vline(xintercept = xintercept, linetype = "dashed", col = "grey48")
    }
    if(length(xintercept) == 2){
      pMain = ggplot(results, aes(x = results[, moderator], y = results[, effect])) +
        labs(x = moderator, y = ylab) +
        geom_smooth(aes(x = as.numeric(levels(results[, moderator]))[results[, moderator]], y = results[, effect]), formula = y ~ x, method = "loess", se = FALSE) +
        geom_ribbon(aes(ymin = predict(loess.CIL), ymax = predict(loess.CIU), x = as.numeric(levels(results[, moderator]))[results[, moderator]], fill = "band"), alpha = 0.5, show.legend = FALSE) +
        scale_fill_manual("",values="grey") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(min(as.numeric(levels(results[, moderator]))[results[, moderator]]), max(as.numeric(levels(results[, moderator]))[results[, moderator]]))) + 
        geom_point(aes(x = xintercept[1], y = 0), col = "grey18") +
        geom_point(aes(x = xintercept[2], y = 0), col = "grey18") +
        geom_text(aes(x = xintercept[1], y = 0, label = paste0("(", round(xintercept[1], 2), ", 0)")), size = 3,  col = "grey18", hjust = 0.5, vjust = 1.5) +
        geom_text(aes(x = xintercept[2], y = 0, label = paste0("(", round(xintercept[2], 2), ", 0)")), size = 3,  col = "grey18", hjust = 0.5, vjust = -1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", col = "grey48") +
        geom_vline(xintercept = xintercept, linetype = "dashed", col = "grey48")
    }
  }
  if(is.dist.moderator){
    subsample.moderators = c(moderators.cont[which(moderators.cont != moderator)], moderators.disc)
    subsample.moderators.values = c(other.mod.cont.values, other.mod.disc.values)
    if(sum(is.na(subsample.moderators.values)) > 0){
      subsample.moderators = subsample.moderators[-which(is.na(subsample.moderators.values))]
      subsample.moderators.values = subsample.moderators.values[-which(is.na(subsample.moderators.values))]
    }
    if(length(moderators.cont) + length(moderators.disc) > 1){
      if(length(subsample.moderators) > 0){
        subsample = NULL
        for(i in 1:length(subsample.moderators)){
          subsample[[i]] = which(data[, subsample.moderators[i]] == subsample.moderators.values[i])
        }
        subsample = Reduce(intersect, subsample)
      } else {
        subsample = 1:nrow(data)
      }
      
      if(length(subsample) >= 30){
        dens = density(data[subsample, moderator])
        df = data.frame(dens.x = dens$x, dens.y = dens$y)
        df = df[-which(df$dens.x < min(as.numeric(levels(results[, moderator]))[results[, moderator]])|df$dens.x > max(as.numeric(levels(results[, moderator]))[results[, moderator]])), ]
        probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
        quantiles = quantile(data[subsample, moderator], prob = probs)
        df$qquant = factor(findInterval(df$dens.x, quantiles))
        if(length(subsample.moderators) == 0){
          caption = paste0("1. Top: Conditional ", effect, " as a function of ", moderator, ". \n2. Bottom: Sample distribution of ", moderator, ".")
        } else {
          caption = paste0("1. Top: Conditional ", effect, " as a function of ", moderator, " given ", paste(subsample.moderators, "=", subsample.moderators.values, collapse = ", "), ". \n2. Bottom: Sample distribution of ", moderator, " given ", paste(subsample.moderators, "=", subsample.moderators.values, collapse = ", "), ".")
        }
        pBottom = ggplot(df, aes(x = dens.x, y = dens.y)) +
          theme(plot.caption = element_text(hjust = 0, face = "italic")) + 
          labs(x = moderator, y = "density", caption = caption) +
          geom_line() + geom_ribbon(aes(ymin = 0, ymax = dens.y, fill = qquant)) + scale_x_continuous(breaks = quantiles, limits = c(min(as.numeric(levels(results[, moderator]))[results[, moderator]]), max(as.numeric(levels(results[, moderator]))[results[, moderator]]))) + scale_fill_brewer(guide = "none")
      } else {
        pBottom = NULL
        if(length(subsample.moderators) == 0){
          warning(paste0("The sample distribution of ", moderator, " is not generated at the bottom because the size of the sample is smaller than 30."))
        } else {
          warning(paste0("The sample distribution of ", moderator, " given ", paste(subsample.moderators, "=", subsample.moderators.values, collapse = ", "), " is not generated at the bottom because the size of the subsample is smaller than 30."))
        }
      }
    } else {
      dens = density(data[, moderator])
      df = data.frame(dens.x = dens$x, dens.y = dens$y)
      df = df[-which(df$dens.x < min(as.numeric(levels(results[, moderator]))[results[, moderator]])|df$dens.x > max(as.numeric(levels(results[, moderator]))[results[, moderator]])), ]
      probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
      quantiles = quantile(data[, moderator], prob = probs)
      df$qquant = factor(findInterval(df$dens.x, quantiles))
      caption = paste0("1. Top: ", effect, " as a function of ", moderator, ". \n2. Bottom: Sample distribution of ", moderator, ".")
      pBottom = ggplot(df, aes(x = dens.x, y = dens.y)) +
        theme(plot.caption = element_text(hjust = 0, face = "italic")) + 
        labs(x = moderator, y = "density", caption = caption) +
        geom_line() + geom_ribbon(aes(ymin = 0, ymax = dens.y, fill = qquant)) + scale_x_continuous(breaks = quantiles, limits = c(min(as.numeric(levels(results[, moderator]))[results[, moderator]]), max(as.numeric(levels(results[, moderator]))[results[, moderator]]))) + scale_fill_brewer(guide = "none")
    } 
    if(!is.null(pBottom)){
      plot_grid(pMain, pBottom, ncol = 1, rel_heights = c(3, 1), align = "v")
    } else {
      if(length(subsample.moderators) == 0){
        caption = paste0("Conditional ", effect, " as a function of ", moderator, ".")
      } else {
        caption = paste0("Conditional ", effect, " as a function of ", moderator, " given ", paste(subsample.moderators, "=", subsample.moderators.values, collapse = ", "), ".")
      }
      pMain + 
        labs(caption = caption) +
        theme(plot.caption = element_text(hjust = 0, face= "italic")) 
    }
  } else {
    if(moderator.scale == "discrete"){
      subsample.moderators = c(moderators.disc[which(moderators.disc != moderator)], moderators.cont)
      subsample.moderators.values = c(other.mod.disc.values, other.mod.cont.values)
      if(sum(is.na(subsample.moderators.values)) > 0){
        subsample.moderators = subsample.moderators[-which(is.na(subsample.moderators.values))]
        subsample.moderators.values = subsample.moderators.values[-which(is.na(subsample.moderators.values))]
      }
      if(length(moderators.cont) + length(moderators.disc) > 1){
        if(length(subsample.moderators) == 0){
          caption = paste0("Conditional ", effect, " as a function of ", moderator, ".")
        } else {
          caption = paste0("Conditional ", effect, " as a function of ", moderator, " given ", paste(subsample.moderators, "=", subsample.moderators.values, collapse = ", "), ".")
        }
      } else {
        caption = paste0(effect, " as a function of ", moderator, ".")
      }
    }
    if(moderator.scale == "continuous"){
      subsample.moderators = c(moderators.cont[which(moderators.cont != moderator)], moderators.disc)
      subsample.moderators.values = c(other.mod.cont.values, other.mod.disc.values)
      if(sum(is.na(subsample.moderators.values)) > 0){
        subsample.moderators = subsample.moderators[-which(is.na(subsample.moderators.values))]
        subsample.moderators.values = subsample.moderators.values[-which(is.na(subsample.moderators.values))]
      }
      if(length(moderators.cont) + length(moderators.disc) > 1){
        if(length(subsample.moderators) == 0){
          caption = paste0("Conditional ", effect, " as a function of ", moderator, ".")
        } else {
          caption = paste0("Conditional ", effect, " as a function of ", moderator, " given ", paste(subsample.moderators, "=", subsample.moderators.values, collapse = ", "), ".")
        }
      } else {
        caption = paste0(effect, " as a function of ", moderator, ".")
      }
    }
    pMain + 
      labs(caption = caption) +
      theme(plot.caption = element_text(hjust = 0, face= "italic")) 
  }
}

#' Simulation-Based Sensitivity Analysis Table for Causal Moderated Mediation Analysis
#' 
#' modmed.sens' is used to evaluate the sensitivity of the estimated causal effects obtained from \code{modmed} function to potential violations of the ignorability assumptions from the frequentist perspective. It estimates the causal effects after adjusting for an unmeasured pretreatment confounder, U, with a specified degree of confounding. In a randomized experiment, the degree of confounding is evaluated via two sensitivity parameters, the coefficient of U in the mediator model and that in the outcome model, given the specified prior distribution of U. When the treatment is not randomized, an additional sensitivity parameter is introduced -- the coefficient of U in the treatment model. The treatment, mediator, outcome, and unmeasured pretreatment confounder could be either binary or continuous. 
#' @param object Output from the \code{modmed} function.
#' @param sens.effect A vector of effects whose sensitivity will be assessed (string).  It is a vector of all the population average, conditional, and moderated mediation effects in the output of modmed function by default (i.e., if sens.effect is not specified). It can also be specified as a subvector of the default. For example, if a researcher is mainly interested in the sensitivity of moderated indirect effects, sens.effect can be specified as c("TIE.dif", "PIE.dif"). 
#' @param range.b.m The range of the sensitivity parameter that represents the slope of the unmeasured pretreatment confounder in the standardized mediator model, in which both the independent and dependent variables are standardized. If the dependent variable is binary, its latent index is standardized instead. E.g., it can be specified as c(-2, 2). If NULL, the upper bound of the range is determined by 2 times the maximum magnitude of the standardized conditional effect of the existing pretreatment confounders with the mediator. The lower bound is the negative of the upper bound. The default is NULL. 
#' @param range.b.y The range of the sensitivity parameter that represents the slope of the unmeasured pretreatment confounder in the standardized outcome model, in which both the independent and dependent variables are standardized. If the dependent variable is binary, its latent index is standardized instead. E.g., it can be specified as c(-2, 2). If NULL, the upper bound of the range is determined by 2 times the maximum magnitude of the standardized conditional effect of the existing pretreatment confounders with the outcome. The lower bound is the negative of the upper bound. The default is NULL.
#' @param grid.b.m The horizontal dimension of the grid. The default is 10. Increase the number for more smooth curves. 
#' @param grid.b.y The vertical dimension of the grid. The default is 10. Increase the number for more smooth curves.
#' @param U.scale The scale of the unobserved pretreatment confounder (string). Can be "continuous" or "binary". The default is "binary".
#' @param p.u This needs to be specified only if U.scale = "binary". The prior probability of the unobserved pretreatment confounder if it is binary. The default is 0.5.
#' @param sigma.u This needs to be specified only if U.scale = "continuous". The standard deviation of the prior distribution of the unobserved pretreatment confounder if it is continuous. The default is 1.
#' @param t.rand TRUE if the treatment is randomized, and FALSE if not. The default is TRUE. 
#' @param t.model A list. The default is NULL. This needs to be specified only if t.rand = FALSE because a treatment model is required for the derivation of the conditional distribution of the unmeasured pretreatment confounder. The names of the elements in t.model must include intercept and the pretreatment covariates that predict the treatment and are contained in the union of covariates.disc, and covariates.cont. Each element of the list is a vector of the names of the moderators (string) of the coefficient of the main model predictor as represented by the name of the element. The moderators of the intercept must include all the moderators in the treatment model. If a main model coefficient is not moderated, then the corresponding vector should be specified as NULL. The moderators should be contained in the union of moderators.disc, and moderators.cont. A moderator in the mediator and outcome models does not necessarily moderate any coefficient in the treatment model. If it only has its main effect in the treatment model, it will be included in the vector of the intercept in t.model.
#' @param t.scale The scale of the treatment (string). Can be "continuous" or "binary". The default is "binary".
#' @param b.t This needs to be specified only if t.rand = FALSE. The value of the sensitivity parameter that represents the slope of the unmeasured pretreatment confounder in the standardized treatment model, in which both the independent and dependent variables are standardized. Please make sure that 1 - b.t^2 is not negative. If the dependent variable is binary, its latent index is standardized instead. When t.rand = TRUE, b.t is NULL, because when the treatment is randomized, there is no need to specify this sensitivity parameter. Hence, the default of b.t is NULL.
#' @param iter The number of iterations in the stochastic EM algorithm for generating the unobserved pretreatment confounder. The default is 10. 
#' @param nsim The number of random draws of the unobserved pretreatment confounder in each cell. The default is 5. Increase the number for more smooth curves.
#' @param ncore The number of cores for parallel computing. The default is the number of CPU cores on the current host minus 1. One core is saved for users to run other programs on the computer while running the R function. 
#' @return A list containing
#' \item{X.coef.plot}{the sopes of the observed pretreatment confounders in the standardized mediator and outcome models (which are used to calibrate the strength of the sensitivity parameters in the sensitivity plots). }
#' \item{b.t}{same as specified in the \code{modmed.sens} function if b.t is not NULL.}
#' \item{range.b.m, range.b.y}{the value range of each sensitivity parameter.}
#' \item{b.m.all, b.y.all}{all the specific values of the sensitivity parameters for constructing the grids.}
#' \item{results.new}{tables for the causal effect estimates, standard errors, and p values in all the grids, after the adjustment of the unobserved pretreatment confounder.}
#' @author Xu Qin and Fan Yang
#' @references Qin, X., & Yang, F. (2020). Simulation-Based Sensitivity Analysis for Causal Mediation Studies.
#' @export
#' @importFrom stats as.formula binomial coef fitted glm lm pnorm dnorm rnorm rbinom predict model.matrix integrate qnorm sd sigma var
#' @importFrom distr AbscontDistribution r
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel detectCores makeCluster stopCluster 
#' @importFrom foreach foreach %dopar%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' \donttest{
#' data(newws)
#' modmed.results = modmed(data = newws, treatment = "treat", mediator = "emp", outcome = "depression", covariates.disc = c("emp_prior", "nevmar", "hispanic", "nohsdip"), covariates.cont = c("workpref", "attitude", "depress_prior"), moderators.disc = "CHCNT", moderators.cont = "ADCPC", m.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), y.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), mediator = c("ADCPC", "CHCNT"), tm = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), comp.mod.disc.values = 3, ref.mod.disc.values = 2, comp.mod.cont.values = 5050, ref.mod.cont.values = 5050, m.scale = "binary", y.scale = "continuous") 
#' sens.results = modmed.sens(modmed.results, sens.effect = "TIE.ref", U.scale = "binary", grid.b.m = 2, grid.b.y = 2, iter = 2, nsim = 2, ncore = 1)
#' sens.results = modmed.sens(modmed.results, sens.effect = "TIE.ref", U.scale = "binary", t.model = list(intercept = c("ADCPC", "CHCNT"), emp_prior = "ADCPC", nevmar = "CHCNT", hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), grid.b.m = 2, grid.b.y = 2, iter = 2, nsim = 2, ncore = 1)
#' }
#' 
modmed.sens = function(object, sens.effect = rownames(object$effects)[-which(rownames(object$effects) %in% c("TE", "TE.ref", "TE.dif"))], range.b.m = NULL, range.b.y = NULL, grid.b.m = 10, grid.b.y = 10, U.scale = "binary", p.u = 0.5, sigma.u = 1, t.rand = TRUE, t.model = NULL, t.scale = "binary", b.t = NULL, iter = 10, nsim = 5, ncore = 2){
  if(!is.null(b.t)){
    if(1 - b.t^2 < 0)
      stop("b.t is the slope of the unmeasured pretreatment confounder in the standardized treatment model, in which both the independent and dependent variables are standardized. Please make sure that 1 - b.t^2 is not negative.")
    b.t.ori = b.t
  }
  if("TE" %in% sens.effect|"TE.ref" %in% sens.effect|"TE.dif" %in% sens.effect)
    stop("Sensitivity analysis here is for mediation effects, not including the total treatment effect.")
  if(any(!sens.effect %in% rownames(object$effects)))
    stop("sens.effect must be a subvector of the effects in the output of modmed function.")
  args.full = object$args.full
  data = args.full$data
  outcome = args.full$outcome
  t = treatment = args.full$treatment
  m = mediator = args.full$mediator
  y.scale = args.full$y.scale
  m.scale = args.full$m.scale
  args.full$is.U = TRUE
  args.full$sens.effect = sens.effect
  args.full$object = object
  m.moderators = object$m.moderators
  y.moderators = object$y.moderators
  m.predictors = object$m.predictors
  y.predictors = object$y.predictors
  conf.level = args.full$conf.level
  formula.m = object$formula.m
  formula.y = object$formula.y
  m.predictors.new = object$m.predictors.new
  y.predictors.new = object$y.predictors.new
  
  if(!t.rand){
    if(is.null(t.model))
      stop("Please specify t.model if t.rand = FALSE")
    if(is.null(b.t))
      stop("Please specify b.t if t.rand = FALSE")
    
    t.predictors = names(t.model)
    if(!"intercept" %in% t.predictors)
      stop("The element names in t.model must include 'Intercept'.")  
    t.predictors = t.predictors[-which(t.predictors == "intercept")]
    
    t.model = append(t.model[which(names(t.model) == "intercept")], t.model)
    t.model = t.model[-which(names(t.model) == "intercept")[2]]
    t.moderators = t.model
    t.moderators.ori = unique(unlist(t.moderators))
    
    confounders = c(object$args.full$covariates.disc, object$args.full$covariates.cont)
    moderators = c(object$args.full$moderators.disc, object$args.full$moderators.cont)
    if(any(moderators %in% t.predictors)){
      stop("The element names in t.model cannot include moderators. In other words, please do not include moderators in the main model.")
    } else {
      if(any(!t.predictors %in% confounders))
        stop("The pretreatment covariates that are not moderators in the treatment model (i.e., the element names that are not intercept in t.model) should be contained in the union of covariates.disc and covariates.cont specified in modmed().")
    }
    if(any(!t.moderators.ori %in% moderators))
      stop("The moderators in the treatment model should be contained in the union of moderators.disc and moderators.cont.")
    if(!all(is.null(unlist(t.moderators)))){
      if(!all(unique(unlist(t.moderators[-which(names(t.moderators) == "intercept")])) %in% t.moderators[[which(names(t.moderators) == "intercept")]]))
        stop("The moderators of the intercept of t.model must include all the moderators in the treatment model.")
    }
    
    if(!is.null(t.moderators)){
      for(i in 1:length(t.moderators)){
        if(is.null(t.moderators[[i]]))
          t.moderators[[i]] = 1
      }
      formula.t = paste(c(t.predictors, t.moderators[[1]]), collapse = " + ")
      for(i in 2:length(t.moderators))
        formula.t = paste(formula.t, "+", paste(t.moderators[[i]], ":", t.predictors[i - 1], collapse = " + "))
    } else {
      formula.t = paste(t.predictors, collapse = " + ")
    }
  }
  
  # Generation of U based on the stochastic EM algorithm
  genU = function(b.y, b.m, b.t = NULL) {
    # Generate U from its prior distribution first
    if(U.scale == "binary")
      data$Unmeasure = rbinom(nrow(data), 1, p.u) 
    if(U.scale == "continuous")
      data$Unmeasure = rnorm(nrow(data), 0, sigma.u)
    
    coef.y.updated = NULL
    coef.m.updated = NULL
    if(t.rand == FALSE)
      coef.t.updated = NULL
    for(i in 1:iter) {
      if(t.rand == FALSE){
        if(t.scale == "continuous"){
          l.t = lm(as.formula(paste(t, "~", formula.t)), offset = b.t * Unmeasure, data = data)
          coef.t = c(l.t$coef, Unmeasure = b.t)
          sd.t = sigma(l.t) # This is the sd of t conditional on X and Unmeasure.
        }
        if(t.scale == "binary"){
          l.t = glm(as.formula(paste(t, "~", formula.t)), offset = b.t * Unmeasure, family = binomial(link = "probit"), data = data)
          coef.t = c(l.t$coef, Unmeasure = b.t)
        }
      }
      if(y.scale == "continuous"){
        l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data)
        coef.y = c(l.y$coef, Unmeasure = b.y)
        sd.y = sigma(l.y) # This is the sd of Y conditional on t, m, X, and Unmeasure.
      }
      if(y.scale == "binary"){
        l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, family = binomial(link = "probit"), data = data)
        coef.y = c(l.y$coef, Unmeasure = b.y)
      }
      if(m.scale == "continuous"){
        l.m = lm(as.formula(paste(m, "~", formula.m)), offset = b.m * Unmeasure, data = data)
        coef.m = c(l.m$coef, Unmeasure = b.m)
        sd.m = sigma(l.m) # This is the sd of M conditional on t, X, and Unmeasure.
      }
      if(m.scale == "binary"){
        l.m = glm(as.formula(paste(m, "~", formula.m)), offset = b.m * Unmeasure, family = binomial(link = "probit"), data = data)
        coef.m = c(l.m$coef, Unmeasure = b.m)
      }
      
      # For a binary Unmeasure
      if(U.scale == "binary"){
        if(t.rand == FALSE){
          # Conditional probability of T
          # If T is continuous
          if(t.scale == "continuous"){
            mean.tu1 = cbind(model.matrix(l.t), Unmeasure = 1) %*% coef.t
            sd.tu1 = sd.t
            mean.tu0 = cbind(model.matrix(l.t), Unmeasure = 0) %*% coef.t
            sd.tu0 = sd.t
            ptu1 = dnorm(data[, t], mean = mean.tu1, sd = sd.tu1)
            ptu0 = dnorm(data[, t], mean = mean.tu0, sd = sd.tu0)
          }
          # If T is binary
          if(t.scale == "binary"){
            pt1u1 = pnorm(cbind(model.matrix(l.t), Unmeasure = 1) %*% coef.t)
            pt1u0 = pnorm(cbind(model.matrix(l.t), Unmeasure = 0) %*% coef.t)
            ptu1 = pt1u1^data[, t] * (1 - pt1u1)^(1 - data[, t])
            ptu0 = pt1u0^data[, t] * (1 - pt1u0)^(1 - data[, t])
          }
        }
        
        # Conditional probability of Y
        # If Y is continuous
        if(y.scale == "continuous"){
          mean.yu1 = cbind(model.matrix(l.y), Unmeasure = 1) %*% coef.y
          sd.yu1 = sd.y
          mean.yu0 = cbind(model.matrix(l.y), Unmeasure = 0) %*% coef.y
          sd.yu0 = sd.y
          pyu1 = dnorm(data[, outcome], mean = mean.yu1, sd = sd.yu1)
          pyu0 = dnorm(data[, outcome], mean = mean.yu0, sd = sd.yu0)
        }
        # If Y is binary
        if(y.scale == "binary"){
          py1u1 = pnorm(cbind(model.matrix(l.y), Unmeasure = 1) %*% coef.y)
          py1u0 = pnorm(cbind(model.matrix(l.y), Unmeasure = 0) %*% coef.y)
          pyu1 = py1u1^data[, outcome] * (1 - py1u1)^(1 - data[, outcome])
          pyu0 = py1u0^data[, outcome] * (1 - py1u0)^(1 - data[, outcome])
        }
        
        # Conditional probability of M
        # If M is continuous
        if(m.scale == "continuous"){
          mean.mu1 = cbind(model.matrix(l.m), Unmeasure = 1) %*% coef.m
          sd.mu1 = sd.m
          mean.mu0 = cbind(model.matrix(l.m), Unmeasure = 0) %*% coef.m
          sd.mu0 = sd.m
          pmu1 = dnorm(data[, m], mean = mean.mu1, sd = sd.mu1)
          pmu0 = dnorm(data[, m], mean = mean.mu0, sd = sd.mu0)
        }
        # If M is binary
        if(m.scale == "binary"){
          pm1u1 = pnorm(cbind(model.matrix(l.m), Unmeasure = 1) %*% coef.m)
          pm1u0 = pnorm(cbind(model.matrix(l.m), Unmeasure = 0) %*% coef.m)
          pmu1 = pm1u1^data[, m] * (1 - pm1u1)^(1 - data[, m])
          pmu0 = pm1u0^data[, m] * (1 - pm1u0)^(1 - data[, m])
        }
        
        # Conditional probability of Unmeasure  
        if(t.rand == TRUE)
          p = pyu1 * pmu1 * p.u/(pyu1 * pmu1 * p.u + pyu0 * pmu0 * (1 - p.u))
        if(t.rand == FALSE)
          p = pyu1 * pmu1 * p.u * ptu1/(pyu1 * pmu1 * p.u * ptu1 + pyu0 * pmu0 * (1 - p.u) * ptu0)
        Unmeasure = NULL
        for (k in 1:nrow(data)) Unmeasure = c(Unmeasure, rbinom(1, 1, p[k]))
        data$Unmeasure = Unmeasure
      }
      
      # For a continuous Unmeasure
      if(U.scale == "continuous"){
        # If Y is continuous and M is binary
        if(y.scale == "continuous" & m.scale == "binary"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(Unmeasure){
                dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                  (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                integrand = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
              if(t.scale == "binary"){
                integrand = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          
          # Generate random values of Unmeasure
          Unmeasure = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of Unmeasure
            if(t.rand == TRUE){
              conditional.u = function(Unmeasure){
                dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                  (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                conditional.u = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
              if(t.scale == "binary"){
                conditional.u = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            Unmeasure = c(Unmeasure, rdist(1))
          } 
          data$Unmeasure = Unmeasure
        }
        
        if(y.scale == "continuous" & m.scale == "continuous"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(Unmeasure){
                dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                integrand = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
              if(t.scale == "binary"){
                integrand = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          # Generate random values of Unmeasure
          Unmeasure = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of Unmeasure
            if(t.rand == TRUE){
              conditional.u = function(Unmeasure){
                dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                conditional.u = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
              if(t.scale == "binary"){
                conditional.u = function(Unmeasure){
                  dnorm(data[i, outcome], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y, sd = sd.y) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            Unmeasure = c(Unmeasure, rdist(1))
          } 
          data$Unmeasure = Unmeasure       
        }
        
        if(y.scale == "binary" & m.scale == "binary"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(Unmeasure){
                (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                  (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                integrand = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
              if(t.scale == "binary"){
                integrand = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          # Generate random values of Unmeasure
          Unmeasure = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of Unmeasure
            if(t.rand == TRUE){
              conditional.u = function(Unmeasure){
                (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                  (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                conditional.u = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
              if(t.scale == "binary"){
                conditional.u = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    (pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^data[i, m] * (1 - pnorm(c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m))^(1 - data[i, m]) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            Unmeasure = c(Unmeasure, rdist(1))
          } 
          data$Unmeasure = Unmeasure
        }
        
        if(y.scale == "binary" & m.scale == "continuous"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(Unmeasure){
                (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                integrand = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
              if(t.scale == "binary"){
                integrand = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)
                }
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          # Generate random values of Unmeasure
          Unmeasure = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of Unmeasure
            if(t.rand == TRUE){
              conditional.u = function(Unmeasure){
                (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                  dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              if(t.scale == "continuous"){
                conditional.u = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    dnorm(data[i, t], mean = c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t, sd = sd.t) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
              if(t.scale == "binary"){
                conditional.u = function(Unmeasure){
                  (pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^data[i, outcome] * (1 - pnorm(c(model.matrix(l.y)[i, ] %*% l.y$coef) + Unmeasure * b.y))^(1 - data[i, outcome]) *
                    dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + Unmeasure * b.m, sd = sd.m) *
                    (pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^data[i, t] * (1 - pnorm(c(model.matrix(l.t)[i, ] %*% l.t$coef) + Unmeasure * b.t))^(1 - data[i, t]) *
                    dnorm(Unmeasure, mean = 0, sd = sigma.u)/integral[i]
                }
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            Unmeasure = c(Unmeasure, rdist(1))
          } 
          data$Unmeasure = Unmeasure
        }
      }
      
      # Check the convergence of the coefficients
      coef.y.updated = rbind(coef.y.updated, coef.y)
      coef.m.updated = rbind(coef.m.updated, coef.m)
    }
    
    return(list(Unmeasure = Unmeasure, coef.y.updated = coef.y.updated, coef.m.updated = coef.m.updated))
  }
  
  m.moderators.new = m.moderators
  y.moderators.new = y.moderators
  if(!is.null(m.moderators.new)){
    for(i in 1:length(m.moderators.new)){
      if(is.null(m.moderators.new[[i]]))
        m.moderators.new[[i]] = 1
    }
  }
  if(!is.null(y.moderators.new)){
    for(i in 1:length(y.moderators.new)){
      if(is.null(y.moderators.new[[i]]))
        y.moderators.new[[i]] = 1   
    }
  }
  if(is.null(m.moderators.new)){
    formula.m = paste(m.predictors, collapse = " + ")
    if(m.scale == "continuous"){
      l.m = lm(as.formula(paste(m, "~", formula.m)), data = data)
      sd.m = sd(data[, m]) 
    }
    if(m.scale == "binary"){
      l.m = glm(as.formula(paste(m, "~", formula.m)), family = binomial(link = "probit"), data = data)
      sd.m = sqrt(var(predict(l.m, type = "link")) + 1) # The variance of the latent index is the sum of the predictor variance and the assumed error variance (McKelvey and Zavoina, 1975). 
    }
    if(length(m.predictors.new[2:length(m.predictors.new)]) > 1){
      X.coef.plot.m = l.m$coef[m.predictors.new[2:length(m.predictors.new)]]/sd.m * apply(model.matrix(l.m)[, m.predictors.new[2:length(m.predictors.new)]], 2, sd) # X.coef.plot are standardized coefficients
    } else {
      X.coef.plot.m = l.m$coef[m.predictors.new[2:length(m.predictors.new)]]/sd.m * sd(model.matrix(l.m)[, m.predictors.new[2:length(m.predictors.new)]]) # X.coef.plot are standardized coefficients
    }
  } else {
    X.coef.plot.m = NULL
    for(k in 2:length(m.predictors)){
      m.moderators.new[[k + 1]] = 1
      formula.m = paste(c(m.predictors, m.moderators.new[[1]]), collapse = " + ")
      for(i in 2:length(m.moderators.new))
        formula.m = paste(formula.m, "+", paste(m.moderators.new[[i]], ":", m.predictors[i - 1], collapse = " + "))
      if(m.scale == "continuous"){
        l.m = lm(as.formula(paste(m, "~", formula.m)), data = data)
        sd.m = sd(data[, m])
      }
      if(m.scale == "binary"){
        l.m = glm(as.formula(paste(m, "~", formula.m)), family = binomial(link = "probit"), data = data)
        sd.m = sqrt(var(predict(l.m, type = "link")) + 1) # The variance of the latent index is the sum of the predictor variance and the assumed error variance (McKelvey and Zavoina, 1975). 
      }
      if(length(colnames(model.matrix(as.formula(paste("~", m.predictors[k])), data = data))[-1]) > 1){
        X.coef.plot.m = c(X.coef.plot.m, l.m$coef[colnames(model.matrix(as.formula(paste("~", m.predictors[k])), data = data))[-1]]/sd.m * apply(model.matrix(l.m)[, colnames(model.matrix(as.formula(paste("~", m.predictors[k])), data = data))[-1]], 2, sd)) # X.coef.plot are standardized coefficients
      } else {
        X.coef.plot.m = c(X.coef.plot.m, l.m$coef[colnames(model.matrix(as.formula(paste("~", m.predictors[k])), data = data))[-1]]/sd.m * sd(model.matrix(l.m)[, colnames(model.matrix(as.formula(paste("~", m.predictors[k])), data = data))[-1]])) # X.coef.plot are standardized coefficients
      }
    }
  }
  
  if(is.null(y.moderators.new)){
    formula.y = paste(y.predictors, collapse = " + ")
    if(y.scale == "continuous"){
      l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data)
      sd.y = sd(data[, outcome])
    }
    if(y.scale == "binary"){
      l.y = glm(as.formula(paste(outcome, "~", formula.y)), family = binomial(link = "probit"), data = data)
      sd.y = sqrt(var(predict(l.y, type = "link")) + 1) # The variance of the latent index is the sum of the predictor variance and the assumed error variance (McKelvey and Zavoina, 1975). 
    }
    if(!paste0(treatment, ":", mediator) %in% y.predictors){
      if(length(y.predictors.new) > 3){
        X.coef.plot.y = l.y$coef[y.predictors.new[!y.predictors.new %in% c(t, m)]]/sd.y * apply(model.matrix(l.y)[, y.predictors.new[!y.predictors.new %in% c(t, m)]], 2, sd) # X.coef.plot are standardized coefficients
      } else {
        X.coef.plot.y = l.y$coef[y.predictors.new[!y.predictors.new %in% c(t, m)]]/sd.y * sd(model.matrix(l.y)[, y.predictors.new[!y.predictors.new %in% c(t, m)]]) # X.coef.plot are standardized coefficients
      }
    } else {
      if(length(y.predictors.new) > 4){
        X.coef.plot.y = l.y$coef[y.predictors.new[!y.predictors.new %in% c(t, m, paste0(treatment, ":", mediator))]]/sd.y * apply(model.matrix(l.y)[, y.predictors.new[!y.predictors.new %in% c(t, m, paste0(treatment, ":", mediator))]], 2, sd) # X.coef.plot are standardized coefficients
      } else {
        X.coef.plot.y = l.y$coef[y.predictors.new[!y.predictors.new %in% c(t, m, paste0(treatment, ":", mediator))]]/sd.y * sd(model.matrix(l.y)[, y.predictors.new[!y.predictors.new %in% c(t, m, paste0(treatment, ":", mediator))]]) # X.coef.plot are standardized coefficients
      }
    }
  } else {
    X.coef.plot.y = NULL
    if(!paste0(treatment, ":", mediator) %in% y.predictors){
      for(k in which(!y.predictors %in% c(t, m))){
        y.moderators.new[[k + 1]] = 1
        formula.y = paste(c(y.predictors, y.moderators.new[[1]]), collapse = " + ")
        for(i in 2:length(y.moderators.new))
          formula.y = paste(formula.y, "+", paste(y.moderators.new[[i]], ":", y.predictors[i - 1], collapse = " + "))
        if(y.scale == "continuous"){
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data)
          sd.y = sd(data[, outcome])
        }
        if(y.scale == "binary"){
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), family = binomial(link = "probit"), data = data)
          sd.y = sqrt(var(predict(l.y, type = "link")) + 1) # The variance of the latent index is the sum of the predictor variance and the assumed error variance (McKelvey and Zavoina, 1975). 
        }
        if(length(colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]) > 1){
          X.coef.plot.y = c(X.coef.plot.y, l.y$coef[colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]]/sd.y * apply(model.matrix(l.y)[, colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]], 2, sd)) # X.coef.plot are standardized coefficients
        } else {
          X.coef.plot.y = c(X.coef.plot.y, l.y$coef[colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]]/sd.y * sd(model.matrix(l.y)[, colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]])) # X.coef.plot are standardized coefficients
        }
      }
    } else {
      for(k in which(!y.predictors %in% c(t, m, paste0(treatment, ":", mediator)))){
        y.moderators.new[[k + 1]] = 1
        formula.y = paste(c(y.predictors, y.moderators.new[[1]]), collapse = " + ")
        for(i in 2:length(y.moderators.new))
          formula.y = paste(formula.y, "+", paste(y.moderators.new[[i]], ":", y.predictors[i - 1], collapse = " + "))
        if(y.scale == "continuous"){
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data)
          sd.y = sd(data[, outcome])
        }
        if(y.scale == "binary"){
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), family = binomial(link = "probit"), data = data)
          sd.y = sqrt(var(predict(l.y, type = "link")) + 1) # The variance of the latent index is the sum of the predictor variance and the assumed error variance (McKelvey and Zavoina, 1975). 
        }
        if(length(colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]) > 1){
          X.coef.plot.y = c(X.coef.plot.y, l.y$coef[colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]]/sd.y * apply(model.matrix(l.y)[, colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]], 2, sd)) # X.coef.plot are standardized coefficients
        } else {
          X.coef.plot.y = c(X.coef.plot.y, l.y$coef[colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]]/sd.y * sd(model.matrix(l.y)[, colnames(model.matrix(as.formula(paste("~", y.predictors[k])), data = data))[-1]])) # X.coef.plot are standardized coefficients
        }
      }
    }
  }
  X.coef.plot = cbind(X.coef.plot.m = as.numeric(X.coef.plot.m), X.coef.plot.y = as.numeric(X.coef.plot.y[names(X.coef.plot.m)]))
  rownames(X.coef.plot) = names(X.coef.plot.m)
  
  if(is.null(range.b.m))
    range.b.m = c(-max(abs(X.coef.plot[, 1])) * 2, max(abs(X.coef.plot[, 1])) * 2)
  if(is.null(range.b.y))
    range.b.y = c(-max(abs(X.coef.plot[, 2])) * 2, max(abs(X.coef.plot[, 2])) * 2)
  
  if(grid.b.m > 1){
    b.m.all = seq(range.b.m[1], range.b.m[2], length.out = grid.b.m)
  } else {
    b.m.all = range.b.m
  }
  if(grid.b.y > 1){
    b.y.all = seq(range.b.y[1], range.b.y[2], length.out = grid.b.y)
  } else {
    b.y.all = range.b.y
  }
  
  args.full$seed = NULL
  args.full$is.U = TRUE
  args.full.ori = args.full
  
  sd.y.ori = sd.y
  sd.m.ori = sd.m
  grid = NULL
  for(i in 1:grid.b.m){
    for(j in 1:grid.b.y){
      grid = rbind(grid, c(i, j))
    }
  }
  
  cl = makeCluster(ncore)
  on.exit(stopCluster(cl))
  registerDoSNOW(cl)
  
  if(interactive()){
    # Initializes the progress bar
    pb = txtProgressBar(min = 0,      # Minimum value of the progress bar
                        max = nrow(grid) * nsim, # Maximum value of the progress bar
                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                        # width = 50,   # Progress bar width. Defaults to getOption("width")
                        char = "=")   # Character used to create the bar
  }
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  resij = foreach(s.total = 1:(nrow(grid) * nsim), .options.snow = opts) %dopar% {
    modmed = function(
    data,
    treatment,
    mediator,
    outcome,
    covariates.disc = NULL,
    covariates.cont = NULL,
    moderators.disc = NULL,
    moderators.cont = NULL,
    m.model, 
    y.model,
    comp.treatment.value = 1,
    ref.treatment.value = 0,
    comp.mod.disc.values = NULL,
    ref.mod.disc.values = NULL,
    comp.mod.cont.values = NULL,
    ref.mod.cont.values = NULL,
    m.scale = "continuous",
    y.scale = "continuous",
    method = "mc",
    nmc = 1000, 
    nboot = 1000,
    conf.level = 0.95,
    seed = NULL, 
    object = NULL,
    is.U = FALSE,
    sens.effect = NULL,
    b.m = NULL,
    b.y = NULL,
    plot.effect = NULL
    ){
      if(sum(is.na(data)))
        stop("Users need to impute any missing values in the data before running the function.")
      if(is.character(unique(data[, treatment])))
        stop("If a treatment is binary, please code it as 0 and 1.")
      data[, treatment] = as.numeric(as.character(data[, treatment])) # Otherwise, if the treatment is categorical, its name may change after model.matrix(), which will cause problems. Without as.character, the values would change.
      Unmeasure = NULL # Otherwise no visible binding for global variable
      
      confounders = c(covariates.disc, covariates.cont)
      moderators = c(moderators.disc, moderators.cont)
      if(any(moderators %in% confounders))
        stop("Please do not include moderators in covariates.disc or covariates.cont")
      ref.mod.values = c(ref.mod.disc.values, ref.mod.cont.values)
      comp.mod.values = c(comp.mod.disc.values, comp.mod.cont.values)
      
      if(!is.U){
        m.predictors = names(m.model)
        if(!"intercept" %in% m.predictors)
          stop("The element names in m.model must include 'Intercept'.")  
        m.predictors = m.predictors[-which(m.predictors == "intercept")]
        if(treatment %in% m.predictors)
          stop("Please specify the element name that corresponds to the treatment in m.model as 'treatment' rather than its actual variable name.")
        if(!"treatment" %in% m.predictors)
          stop("The element names in m.model must include 'treatment'. In other words, predictors in the main model of the mediator must include the treatment.")  
        m.predictors[which(m.predictors == "treatment")] = treatment
        
        y.predictors = names(y.model)
        if(!"intercept" %in% y.predictors)
          stop("The element names in y.model must include 'Intercept'.") 
        y.predictors = y.predictors[-which(y.predictors == "intercept")]
        if(treatment %in% y.predictors)
          stop("Please specify the element name that corresponds to the treatment in y.model as 'treatment' rather than its actual variable name.")  
        if(!"treatment" %in% y.predictors)
          stop("The element names in y.model must include 'treatment'. In other words, predictors in the main model of the outcome must include the treatment.")
        if(mediator %in% y.predictors)
          stop("Please specify the element name that corresponds to the mediator in y.model as 'mediator' rather than its actual variable name.") 
        if(!"mediator" %in% y.predictors)
          stop("The element names in y.model must include 'mediator'. In other words, predictors in the main model of the outcome must include the mediator.")
        y.predictors[which(y.predictors == "treatment")] = treatment
        y.predictors[which(y.predictors == "mediator")] = mediator
        
        if(any(moderators %in% m.predictors)|any(moderators %in% y.predictors)){
          stop("The element names in m.model or y.model cannot include moderators. In other words, please do not include moderators in the main model.")
        } else {
          if(any(!m.predictors %in% c(treatment, confounders)))
            stop("The pretreatment covariates that are not moderators in the mediator model (i.e., the element names that are not intercept or treatment in m.model) should be contained in the union of covariates.disc and covariates.cont.")
          if(any(!y.predictors %in% c(treatment, mediator, confounders))){
            if(any(!y.predictors %in% c(treatment, mediator, "tm", confounders)))
              stop("The pretreatment covariates that are not moderators in the outcome model (i.e., the element names that are not intercept, treatment, or mediator in y.model) should be contained in the union of covariates.disc and covariates.cont. If there is treatment-by-mediator interaction, it must be named as 'tm'.")
          }
        }
        
        if("tm" %in% y.predictors){
          y.predictors[which(y.predictors == "tm")] = paste0(treatment, ":", mediator)
        }   
        
        m.model = append(m.model[which(names(m.model) == "intercept")], m.model)
        m.model = m.model[-which(names(m.model) == "intercept")[2]]
        m.moderators = m.model
        y.model = append(y.model[which(names(y.model) == "intercept")], y.model)
        y.model = y.model[-which(names(y.model) == "intercept")[2]]
        y.moderators = y.model
        m.moderators.ori = unique(unlist(m.moderators))
        y.moderators.ori = unique(unlist(y.moderators))
        
        if(any(!unique(c(m.moderators.ori, y.moderators.ori)) %in% moderators)|any(!moderators %in% unique(c(m.moderators.ori, y.moderators.ori))))
          stop('The union of the two sets of moderators in m.model and y.model should be the same as the union of moderators.disc and moderators.cont.')
        
        if(!is.null(c(moderators.disc, moderators.cont))){
          for(i in 1:length(moderators)){
            if(!moderators[i] %in% unique(unlist(m.moderators)) & !moderators[i] %in% unique(unlist(y.moderators)))
              stop("A moderator specified in moderators.disc and moderators.cont must moderate at least one slope in either the mediator model or the outcome model.")
          }
        } else {
          if(!is.null(c(unique(unlist(m.moderators)), unique(unlist(y.moderators)))))
            stop("A moderator in the model needs to be specified in moderators.disc/moderators.cont.")
        }
        
        
        if(!all(is.null(unlist(m.moderators)))){
          if(!all(unique(unlist(m.moderators[-which(names(m.moderators) == "intercept")])) %in% m.moderators[[which(names(m.moderators) == "intercept")]]))
            stop("The moderators of the intercept of m.model must include all the moderators in the mediator model.")
        }
        if(!all(is.null(unlist(y.moderators)))){
          if(!all(unique(unlist(y.moderators[-which(names(y.moderators) == "intercept")])) %in% y.moderators[[which(names(y.moderators) == "intercept")]]))
            stop("The moderators of the intercept of y.model must include all the moderators in the outcome model.")
        }
        
        # if(is.null(ref.mod.disc.values) & !is.null(moderators.disc)){
        #   stop("If moderators.disc is not NULL, ref.mod.disc.values should not be NULL. The length and order of ref.mod.disc.values should be the same as moderators.disc. If one does not want to condition some moderators on specific values, one may specify their values to be NA.")
        # }
        if(!is.null(ref.mod.disc.values) & is.null(moderators.disc)){
          stop("If moderators.disc is NULL, ref.mod.disc.values should also be NULL. The length and order of ref.mod.disc.values should be the same as moderators.disc.")
        }
        # if(is.null(ref.mod.cont.values) & !is.null(moderators.cont)){
        #   stop("If moderators.cont is not NULL, ref.mod.cont.values should not be NULL. The length and order of ref.mod.cont.values should be the same as moderators.cont. If one does not want to condition some moderators on specific values, one may specify their values to be NA.")
        # }
        if(!is.null(ref.mod.cont.values) & is.null(moderators.cont)){
          stop("If moderators.cont is NULL, ref.mod.cont.values should also be NULL. The length and order of ref.mod.cont.values should be the same as moderators.cont.")
        }
        
        if(!is.null(ref.mod.disc.values)){
          if(length(ref.mod.disc.values) != length(moderators.disc))
            stop("Please make sure that the length and order of ref.mod.disc.values should be the same as moderators.disc") 
        }
        if(!is.null(ref.mod.cont.values)){
          if(length(ref.mod.cont.values) != length(moderators.cont))
            stop("Please make sure that the length and order of ref.mod.cont.values should be the same as moderators.cont")
        }
        
        if(is.null(ref.mod.disc.values) & !is.null(comp.mod.disc.values))
          stop("Please make sure to specify ref.mod.disc.values before specifying comp.mod.disc.values")
        if(is.null(ref.mod.cont.values) & !is.null(comp.mod.cont.values))
          stop("Please make sure to specify ref.mod.cont.values before specifying comp.mod.cont.values")
        
        if(!is.null(comp.mod.disc.values) & is.null(moderators.disc)){
          stop("If moderators.disc is NULL, comp.mod.disc.values should also be NULL. The length and order of comp.mod.disc.values should be the same as moderators.disc.")
        }
        if(!is.null(comp.mod.cont.values) & is.null(moderators.cont)){
          stop("If moderators.cont is NULL, comp.mod.cont.values should also be NULL. The length and order of comp.mod.cont.values should be the same as moderators.cont.")
        }
        
        if(!is.null(comp.mod.disc.values)){
          if(length(comp.mod.disc.values) != length(moderators.disc))
            stop("Please make sure that the length and order of comp.mod.disc.values should be the same as moderators.disc")
        }
        if(!is.null(comp.mod.cont.values)){
          if(length(comp.mod.cont.values) != length(moderators.cont))
            stop("Please make sure that the length and order of comp.mod.cont.values should be the same as moderators.cont")
        }
        
        if(!is.null(ref.mod.disc.values) & !any(is.na(ref.mod.disc.values))){
          for(i in 1:length(ref.mod.disc.values)){
            if(!ref.mod.disc.values[i] %in% unique(data[, moderators.disc[i]]))
              stop("Please specify each of ref.mod.disc.values to be one of the values of the corresponding moderator in the data")
            if(!is.null(comp.mod.disc.values)){
              if(!comp.mod.disc.values[i] %in% unique(data[, moderators.disc[i]]))
                stop("Please specify each of comp.mod.disc.values to be one of the values of the corresponding moderator in the data")
            } 
          }
        }
        
        if(!is.null(covariates.disc)){
          for(i in 1:length(covariates.disc)){
            if(!is.factor(data[, covariates.disc[i]]))
              data[, covariates.disc[i]] = as.factor(data[, covariates.disc[i]])
          }
        }  
        if(!is.null(moderators.disc)){
          for(i in 1:length(moderators.disc)){
            if(!is.factor(data[, moderators.disc[i]]))
              data[, moderators.disc[i]] = as.factor(data[, moderators.disc[i]])
          }
        }
        
        m.moderators.new = m.moderators
        y.moderators.new = y.moderators
        if(!is.null(m.moderators.new)){
          for(i in 1:length(m.moderators.new)){
            if(is.null(m.moderators.new[[i]]))
              m.moderators.new[[i]] = 1
          }
        }
        if(!is.null(y.moderators.new)){
          for(i in 1:length(y.moderators.new)){
            if(is.null(y.moderators.new[[i]]))
              y.moderators.new[[i]] = 1
          }
        }
        
        if(!is.null(m.moderators.new)){
          formula.m = paste(c(m.predictors, m.moderators.new[[1]]), collapse = " + ")
          for(i in 2:length(m.moderators.new))
            formula.m = paste(formula.m, "+", paste(m.moderators.new[[i]], ":", m.predictors[i - 1], collapse = " + "))
        } else {
          formula.m = paste(m.predictors, collapse = " + ")
        }
        
        if(!is.null(y.moderators.new)){
          formula.y = paste(c(y.predictors, y.moderators.new[[1]]), collapse = " + ")
          for(i in 2:length(y.moderators.new))
            formula.y = paste(formula.y, "+", paste(y.moderators.new[[i]], ":", y.predictors[i - 1], collapse = " + "))
        } else {
          formula.y = paste(y.predictors, collapse = " + ")
        }
        
        new.moderators.predictors = function(predictors, moderators.ori.list){
          if(!is.null(unique(unlist(moderators.ori.list)))){
            v1 = colnames(model.matrix(as.formula(paste("~", paste(moderators.ori.list[[1]], collapse = " + "))), data = data))[-1]
            k = 1
            for(i in 1:length(moderators.ori.list)){
              if(!is.null(moderators.ori.list[[i]]))
                moderators.ori.list[[i]] = colnames(model.matrix(as.formula(paste("~", paste(moderators.ori.list[[i]], collapse = " + "))), data = data))[-1]
            }
          }
          predictors.new = NULL
          for(i in 1:length(predictors)){
            predictors.new = c(predictors.new, colnames(model.matrix(as.formula(paste("~", paste(predictors[i], collapse = " + "))), data = data))[-1])
            if(!is.null(unique(unlist(moderators.ori.list)))){
              if(!predictors[i] %in% covariates.disc){
                assign(paste0("v", k + 1), moderators.ori.list[[i + 1]])
                k = k + 1
              }
              if(predictors[i] %in% covariates.disc){
                predictors.cat = colnames(model.matrix(as.formula(paste("~", paste(predictors[[i]], collapse = " + "))), data = data))[-1]
                for(j in 1:length(predictors.cat)){
                  assign(paste0("v", k + j), moderators.ori.list[[i + 1]])
                }
                k = k + length(predictors.cat)
              }
            }
          }
          if(!is.null(unique(unlist(moderators.ori.list)))){
            moderators.new = rep(list(NULL), k)
            for(i in 1:k){
              if(!is.null(get(paste0("v", i))))
                moderators.new[[i]] = get(paste0("v", i))
            }
          }
          if(!is.null(unique(unlist(moderators.ori.list)))){
            return(list(moderators.new = moderators.new, predictors.new = predictors.new))
          } else {
            return(list(predictors.new = predictors.new))
          }
        }
        m.new.moderators.predictors = new.moderators.predictors(m.predictors, m.moderators)
        m.moderators.new = m.new.moderators.predictors$moderators.new
        m.predictors.new = m.new.moderators.predictors$predictors.new
        y.new.moderators.predictors = new.moderators.predictors(y.predictors, y.moderators)
        y.moderators.new = y.new.moderators.predictors$moderators.new
        y.predictors.new = y.new.moderators.predictors$predictors.new
        
        ref.m.mod.values = NULL
        comp.m.mod.values = NULL
        if(!is.null(m.moderators.ori)){
          for(i in 1:length(m.moderators.ori)){
            if(!is.null(ref.mod.values))
              ref.m.mod.values = c(ref.m.mod.values, ref.mod.values[which(moderators %in% m.moderators.ori[i])])
            if(!is.null(comp.mod.values))
              comp.m.mod.values = c(comp.m.mod.values, comp.mod.values[which(moderators %in% m.moderators.ori[i])])
          }
        }
        ref.y.mod.values = NULL
        comp.y.mod.values = NULL
        if(!is.null(y.moderators.ori)){
          for(i in 1:length(y.moderators.ori)){
            if(!is.null(ref.mod.values))
              ref.y.mod.values = c(ref.y.mod.values, ref.mod.values[which(moderators %in% y.moderators.ori[i])])
            if(!is.null(comp.mod.values))
              comp.y.mod.values = c(comp.y.mod.values, comp.mod.values[which(moderators %in% y.moderators.ori[i])])
          }
        }
        
        # Function for importing moderators.values to the data for prediction
        predict.data.mod = function(predict.data, predict.moderators, predict.moderators.values){
          if(!is.null(predict.moderators.values) & sum(is.na(predict.moderators.values)) > 0){
            predict.moderators = predict.moderators[-which(is.na(predict.moderators.values))]
            predict.moderators.values = predict.moderators.values[-which(is.na(predict.moderators.values))]
          }
          if(!is.null(predict.moderators.values)){
            for(i in 1:length(predict.moderators)){
              moderator = predict.moderators[i]
              if(moderator %in% colnames(predict.data)){
                predict.data[, moderator] = as.numeric(predict.moderators.values[i])
              } else {
                moderator.value = paste0(predict.moderators[i], predict.moderators.values[i])
                name.moderator.value = colnames(predict.data)[grep(moderator, colnames(predict.data))]
                if(moderator.value %in% colnames(predict.data)){
                  predict.data[, moderator.value] = 1
                  if(length(name.moderator.value) > 1){
                    name.moderator.value.0 = name.moderator.value[which(!name.moderator.value %in% moderator.value)]
                    predict.data[, name.moderator.value.0] = 0
                  }
                } else {
                  name.moderator.value.0 = name.moderator.value[which(!name.moderator.value %in% moderator.value)]
                  predict.data[, name.moderator.value.0] = 0
                }
              }
            }
          }
          
          colnames(predict.data) = gsub("[^[:alnum:]\\:\\s]", "", colnames(predict.data))
          predict.data = as.data.frame(predict.data)
          
          return(predict.data)
        }
        
        predict.m.data = model.matrix(as.formula(paste("~", formula.m)), data = data)
        predict.y.data = model.matrix(as.formula(paste("~", formula.y)), data = data)
        
        predict.m.data.ref = predict.data.mod(predict.m.data, m.moderators.ori, ref.m.mod.values)
        predict.m.data.comp = predict.data.mod(predict.m.data, m.moderators.ori, comp.m.mod.values)
        predict.m.data = predict.data.mod(predict.m.data, m.moderators.ori, NULL)
        predict.y.data.ref = predict.data.mod(predict.y.data, y.moderators.ori, ref.y.mod.values)
        predict.y.data.comp = predict.data.mod(predict.y.data, y.moderators.ori, comp.y.mod.values)
        predict.y.data = predict.data.mod(predict.y.data, y.moderators.ori, NULL)
      } else {
        m.moderators.new = object$m.moderators.new
        y.moderators.new = object$y.moderators.new
        m.predictors.new = object$m.predictors.new
        y.predictors.new = object$y.predictors.new
        formula.m = object$formula.m
        formula.y = object$formula.y
        predict.m.data.ref = object$predict.m.data.ref
        predict.m.data.comp = object$predict.m.data.comp
        predict.m.data = object$predict.m.data
        predict.y.data.ref = object$predict.y.data.ref
        predict.y.data.comp = object$predict.y.data.comp
        predict.y.data = object$predict.y.data
      }
      
      # Fit mediator and outcome model based on the original data set first
      if(m.scale == "binary" & y.scale == "continuous"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data, family = binomial(link = "probit"))
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data)
        }
        if(is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data, family = binomial(link = "probit"))
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data)
        }
        
        if(any(is.na(coef(l.m)))){
          print(coef(l.m))
          stop("NA in coefficients of the mediator model. Please double check.")
        }
        if(any(is.na(coef(l.y)))){
          print(coef(l.y))
          stop("NA in coefficients of the outcome model. Please double check.")
        }
      }
      
      if(m.scale == "continuous" & y.scale == "continuous"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data)
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data)
        }
        if(is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data)
          l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data)
        }
        
        if(any(is.na(coef(l.m)))){
          print(coef(l.m))
          stop("NA in coefficients of the mediator model. Please double check.")
        }
        if(any(is.na(coef(l.y)))){
          print(coef(l.y))
          stop("NA in coefficients of the outcome model. Please double check.")
        }
      }
      
      if(m.scale == "binary" & y.scale == "binary"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data, family = binomial(link = "probit"))
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data, family = binomial(link = "probit"))
        }
        if(is.U){
          l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data, family = binomial(link = "probit"))
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data, family = binomial(link = "probit"))
        }
        
        if(any(is.na(coef(l.m)))){
          print(coef(l.m))
          stop("NA in coefficients of the mediator model. Please double check.")
        }
        if(any(is.na(coef(l.y)))){
          print(coef(l.y))
          stop("NA in coefficients of the outcome model. Please double check.")
        }
      }
      
      if(m.scale == "continuous" & y.scale == "binary"){
        # Step 1. Fit mediator and outcome models
        if(!is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data)
          sd.m = sigma(l.m) 
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data, family = binomial(link = "probit"))
        }
        if(is.U){
          l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data)
          sd.m = sigma(l.m) 
          l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data, family = binomial(link = "probit"))
        }
        
        if(any(is.na(coef(l.m)))){
          print(coef(l.m))
          stop("NA in coefficients of the mediator model. Please double check.")
        }
        if(any(is.na(coef(l.y)))){
          print(coef(l.y))
          stop("NA in coefficients of the outcome model. Please double check.")
        }
      }
      
      set.seed(seed)
      
      # Monte Carlo method
      if(method == 'mc'){
        if(m.scale == "binary" & y.scale == "continuous"){
          # Step 2. Simulate model parameters from their sampling distribution.
          coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
          coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
          
          # Step 3. Predict potential outcomes
          predict.m = function(predict.data, t){
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
            if(!is.U)
              return(pnorm(tcrossprod(predict.data, coef.m.sim)))
            if(is.U)
              return(pnorm((tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)))
          }
          
          predict.y = function(predict.data, t, m) {
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
            predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
            if(!is.U)
              return(tcrossprod(predict.data, coef.y.sim))
            if(is.U)
              return(tcrossprod(predict.data, coef.y.sim) + b.y * data$Unmeasure)
          }
          
          y1m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y1m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
          }
          y0m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y0m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
          }
        }
        
        if(m.scale == "continuous" & y.scale == "continuous"){
          # Step 2. Simulate model parameters from their sampling distribution.
          coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
          coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
          
          # Step 3. Predict potential outcomes
          predict.m = function(predict.data, t){
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
            if(!is.U)
              return(tcrossprod(predict.data, coef.m.sim))
            if(is.U)
              return(tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)
          }
          
          predict.y = function(predict.data, t, m) {
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            j.mediator = which(colnames(predict.data) == gsub("[^[:alnum:]\\:\\s]", "", mediator))
            predict.results = matrix(NA, nrow(predict.data), nmc)
            for(i in 1:nmc){
              predict.data[, j.mediator] = m[, i]
              if(!is.U)
                predict.results[, i] = tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ]))
              if(is.U)
                predict.results[, i] = tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ])) + b.y * data$Unmeasure
            }
            return(predict.results)
          }
          
          y1m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y1m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
          }
          y0m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y0m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
          }
        }
        
        if(m.scale == "binary" & y.scale == "binary"){
          # Step 2. Simulate model parameters from their sampling distribution.
          coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
          coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
          
          # Step 3. Predict potential outcomes
          predict.m = function(predict.data, t){
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
            if(!is.U)
              return(pnorm(tcrossprod(predict.data, coef.m.sim)))
            if(is.U)
              return(pnorm((tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)))
          }
          
          predict.y = function(predict.data, t, m) {
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
            predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
            if(!is.U){
              return(pnorm(tcrossprod(predict.data, coef.y.sim)))
            }
            if(is.U){
              return(pnorm((tcrossprod(predict.data, coef.y.sim) + b.y * data$Unmeasure)))
            }
          }
          
          y1m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y1m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
          }
          y0m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y0m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
          }
        }
        
        if(m.scale == "continuous" & y.scale == "binary"){
          # Step 2. Simulate model parameters from their sampling distribution.
          coef.m.sim = rmvnorm(nmc, mean = coef(l.m), sigma = vcov(l.m))
          coef.y.sim = rmvnorm(nmc, mean = coef(l.y), sigma = vcov(l.y))
          
          # Step 3. Predict potential outcomes
          predict.m = function(predict.data, t){
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
            if(!is.U)
              return(tcrossprod(predict.data, coef.m.sim))
            if(is.U)
              return(tcrossprod(predict.data, coef.m.sim) + b.m * data$Unmeasure)
          }
          
          predict.y = function(predict.data, t, m) {
            predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
            j.mediator = which(colnames(predict.data) == gsub("[^[:alnum:]\\:\\s]", "", mediator))
            predict.results = matrix(NA, nrow(predict.data), nmc)
            for(i in 1:nmc){
              denominator = sqrt(sd.m^2 * (sum(coef.y.sim[i, colnames(coef.y.sim)[grepl(mediator, colnames(coef.y.sim))]]))^2 + 1)
              predict.data[, j.mediator] = m[, i]
              if(!is.U)
                predict.results[, i] = pnorm(tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ]))/denominator)
              if(is.U)
                predict.results[, i] = pnorm((tcrossprod(model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data), t(coef.y.sim[i, ])) + b.y * data$Unmeasure)/denominator)
            }
            return(predict.results)
          }
          
          y1m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y1m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
          }
          y0m1 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
          }
          y0m0 = function(predict.m.data, predict.y.data){
            return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
          }
        }
        
        # Step 4. Calculate the final effects for each sample
        est.y1m1 = y1m1(predict.m.data, predict.y.data)
        est.y1m0 = y1m0(predict.m.data, predict.y.data)
        est.y0m0 = y0m0(predict.m.data, predict.y.data)
        est.y0m1 = y0m1(predict.m.data, predict.y.data)
        
        if(is.null(plot.effect)){
          if(is.U){
            sens.results = NULL
            if("TIE" %in% sens.effect){
              TIE = apply(est.y1m1 - est.y1m0, 2, mean)
              sens.results = cbind(sens.results, TIE)
            }
            if("PDE" %in% sens.effect){
              PDE = apply(est.y1m0 - est.y0m0, 2, mean)
              sens.results = cbind(sens.results, PDE)
            }
            if("PIE" %in% sens.effect){
              PIE = apply(est.y0m1 - est.y0m0, 2, mean)
              sens.results = cbind(sens.results, PIE)
            }
            if("TDE" %in% sens.effect){
              TDE = apply(est.y1m1 - est.y0m1, 2, mean)
              sens.results = cbind(sens.results, TDE)
            }
            if("INT" %in% sens.effect){
              INT = apply(est.y1m1 - est.y1m0, 2, mean) - apply(est.y0m1 - est.y0m0, 2, mean)
              sens.results = cbind(sens.results, INT)
            }
          } else {
            TIE = apply(est.y1m1 - est.y1m0, 2, mean)
            PDE = apply(est.y1m0 - est.y0m0, 2, mean)
            PIE = apply(est.y0m1 - est.y0m0, 2, mean)
            INT = TIE - PIE
            TDE = PDE + INT
            results = cbind(TIE, PIE, PDE, TDE, INT)
          }
        }
        
        if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
          if(is.null(confounders)|
             m.scale == "continuous" & y.scale == "continuous"|
             m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
            if(!is.null(ref.mod.values)){
              if(is.null(plot.effect)){
                if(is.U){
                  if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  }
                  if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  }
                  if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  }
                  if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  }
                } else {
                  y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  TIE.ref = y1m1.ref - y1m0.ref
                  PDE.ref = y1m0.ref - y0m0.ref
                  PIE.ref = y0m1.ref - y0m0.ref
                  INT.ref = TIE.ref - PIE.ref
                  TDE.ref = PDE.ref + INT.ref
                }
              } else {
                if(plot.effect == "TE"){
                  y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  results = y1m1.ref - y0m0.ref
                }
                if(plot.effect == "TIE"){
                  y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  results = y1m1.ref - y1m0.ref
                }
                if(plot.effect == "PIE"){
                  y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  results = y0m1.ref - y0m0.ref
                }
                if(plot.effect == "TDE"){
                  y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  results = y1m1.ref - y0m1.ref
                }
                if(plot.effect == "PDE"){
                  y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  results = y1m0.ref - y0m0.ref
                }
                if(plot.effect == "INT"){
                  y1m1.ref = apply(y1m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y1m0.ref = apply(y1m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m0.ref = apply(y0m0(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  y0m1.ref = apply(y0m1(predict.m.data.ref, predict.y.data.ref), 2, mean)
                  results = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
                }
              }
            }  
            if(!is.null(comp.mod.values)){
              if(is.null(plot.effect)){
                if(is.U){
                  if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m1.dif = apply(y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref, 2, mean)
                  }
                  if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m0.dif = apply(y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref, 2, mean)
                  }
                  if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m0.dif = apply(y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref, 2, mean)
                  }
                  if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m1.dif = apply(y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref, 2, mean)
                  }
                } else {
                  y1m1.dif = apply(y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref, 2, mean)
                  y1m0.dif = apply(y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref, 2, mean)
                  y0m0.dif = apply(y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref, 2, mean)
                  y0m1.dif = apply(y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref, 2, mean)
                  TIE.dif = y1m1.dif - y1m0.dif
                  PDE.dif = y1m0.dif - y0m0.dif
                  PIE.dif = y0m1.dif - y0m0.dif
                  INT.dif = TIE.dif - PIE.dif
                  TDE.dif = PDE.dif + INT.dif
                }
              }
            }
          }
        }
      }
      
      # Bootstrap method
      if(method == 'boot'){
        est.y1m1 = est.y1m0 = est.y0m0 = est.y0m1 = y1m1.ref = y1m0.ref = y0m0.ref = y0m1.ref = y1m1.dif = y1m0.dif = y0m0.dif = y0m1.dif = NULL
        for(b in 1:nboot){
          data.boot = data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
          if(m.scale == "binary" & y.scale == "continuous"){
            # Step 1. Fit mediator and outcome models
            if(!is.U){
              l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data.boot, family = binomial(link = "probit"))
              l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data.boot)
            }
            if(is.U){
              l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot, family = binomial(link = "probit"))
              l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot)
            }
            
            # Step 2. model coefficient estimates
            coef.m = coef(l.m)
            coef.y = coef(l.y)
            
            # Step 3. Predict potential outcomes
            predict.m = function(predict.data, t){
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U)
                return(pnorm(tcrossprod(predict.data, t(coef.m))))
              if(is.U)
                return(pnorm((tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)))
            }
            
            predict.y = function(predict.data, t, m) {
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U)
                return(tcrossprod(predict.data, t(coef.y)))
              if(is.U)
                return(tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)
            }
            
            y1m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y1m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
            }
            y0m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y0m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
            }
          }
          
          if(m.scale == "continuous" & y.scale == "continuous"){
            # Step 1. Fit mediator and outcome models
            if(!is.U){
              l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data.boot)
              l.y = lm(as.formula(paste(outcome, "~", formula.y)), data = data.boot)
            }
            if(is.U){
              l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot)
              l.y = lm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot)
            }
            
            # Step 2. Simulate model parameters from their sampling distribution.
            coef.m = coef(l.m)
            coef.y = coef(l.y)
            
            # Step 3. Predict potential outcomes
            predict.m = function(predict.data, t){
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U)
                return(tcrossprod(predict.data, t(coef.m)))
              if(is.U)
                return(tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)
            }
            
            predict.y = function(predict.data, t, m) {
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U)
                return(tcrossprod(predict.data, t(coef.y)))
              if(is.U)
                return(tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)
            }
            
            y1m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y1m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
            }
            y0m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y0m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
            }
          }
          
          if(m.scale == "binary" & y.scale == "binary"){
            # Step 1. Fit mediator and outcome models
            if(!is.U){
              l.m = glm(as.formula(paste(mediator, "~", formula.m)), data = data.boot, family = binomial(link = "probit"))
              l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data.boot, family = binomial(link = "probit"))
            }
            if(is.U){
              l.m = glm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot, family = binomial(link = "probit"))
              l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot, family = binomial(link = "probit"))
            }
            
            # Step 2. Simulate model parameters from their sampling distribution.
            coef.m = coef(l.m)
            coef.y = coef(l.y)
            
            # Step 3. Predict potential outcomes
            predict.m = function(predict.data, t){
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U)
                return(pnorm(tcrossprod(predict.data, t(coef.m))))
              if(is.U)
                return(pnorm((tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)))
            }
            
            predict.y = function(predict.data, t, m) {
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U){
                return(pnorm(tcrossprod(predict.data, t(coef.y))))
              }
              if(is.U){
                return(pnorm((tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)))
              }
            }
            
            y1m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y1m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = comp.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
            }
            y0m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = comp.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y0m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = 1) * predict.m(predict.m.data, t = ref.treatment.value) + predict.y(predict.y.data, t = ref.treatment.value, m = 0) * (1 - predict.m(predict.m.data, t = ref.treatment.value)))
            }
          }
          
          if(m.scale == "continuous" & y.scale == "binary"){
            # Step 1. Fit mediator and outcome models
            if(!is.U){
              l.m = lm(as.formula(paste(mediator, "~", formula.m)), data = data.boot)
              sd.m = sigma(l.m) 
              l.y = glm(as.formula(paste(outcome, "~", formula.y)), data = data.boot, family = binomial(link = "probit"))
            }
            if(is.U){
              l.m = lm(as.formula(paste(mediator, "~", formula.m)), offset = b.m * Unmeasure, data = data.boot)
              sd.m = sigma(l.m) 
              l.y = glm(as.formula(paste(outcome, "~", formula.y)), offset = b.y * Unmeasure, data = data.boot, family = binomial(link = "probit"))
            }
            
            # Step 2. Simulate model parameters from their sampling distribution.
            coef.m = coef(l.m)
            coef.y = coef(l.y)
            
            # Step 3. Predict potential outcomes
            predict.m = function(predict.data, t){
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U)
                return(tcrossprod(predict.data, t(coef.m)))
              if(is.U)
                return(tcrossprod(predict.data, t(coef.m)) + b.m * data$Unmeasure)
            }
            
            predict.y = function(predict.data, t, m) {
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", treatment)] = t
              predict.data[, gsub("[^[:alnum:]\\:\\s]", "", mediator)] = m
              predict.data = model.matrix(as.formula(paste("~", paste(colnames(predict.data)[-1], collapse = " + "))), data = predict.data)
              if(!is.U){
                return(pnorm(tcrossprod(predict.data, t(coef.y))))
              }
              if(is.U){
                return(pnorm((tcrossprod(predict.data, t(coef.y)) + b.y * data$Unmeasure)))
              }
            }
            
            y1m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y1m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = comp.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
            }
            y0m1 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = comp.treatment.value)))
            }
            y0m0 = function(predict.m.data, predict.y.data){
              return(predict.y(predict.y.data, t = ref.treatment.value, m = predict.m(predict.m.data, t = ref.treatment.value)))
            }
          }
          
          # Step 4. Calculate the final effects for each sample
          est.y1m1 = cbind(est.y1m1, y1m1(predict.m.data, predict.y.data))
          est.y1m0 = cbind(est.y1m0, y1m0(predict.m.data, predict.y.data))
          est.y0m0 = cbind(est.y0m0, y0m0(predict.m.data, predict.y.data))
          est.y0m1 = cbind(est.y0m1, y0m1(predict.m.data, predict.y.data))
          
          if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
            if(is.null(confounders)|
               m.scale == "continuous" & y.scale == "continuous"|
               m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
              if(!is.null(ref.mod.values)){
                if(is.null(plot.effect)){
                  if(is.U){
                    if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                    }
                    if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                    }
                    if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                    }
                    if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
                    }
                  } else {
                    y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                    y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                    y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                    y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
                  }
                } else {
                  if(plot.effect == "TE"){
                    y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                    y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                  }
                  if(plot.effect == "TIE"){
                    y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                    y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                  }
                  if(plot.effect == "PIE"){
                    y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
                    y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                  }
                  if(plot.effect == "TDE"){
                    y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                    y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
                  }
                  if(plot.effect == "PDE"){
                    y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                    y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                  }
                  if(plot.effect == "INT"){
                    y1m1.ref = cbind(y1m1.ref, y1m1(predict.m.data.ref, predict.y.data.ref))
                    y1m0.ref = cbind(y1m0.ref, y1m0(predict.m.data.ref, predict.y.data.ref))
                    y0m0.ref = cbind(y0m0.ref, y0m0(predict.m.data.ref, predict.y.data.ref))
                    y0m1.ref = cbind(y0m1.ref, y0m1(predict.m.data.ref, predict.y.data.ref))
                  }
                }
                
              }  
              if(!is.null(comp.mod.values)){
                if(is.null(plot.effect)){
                  if(is.U){
                    if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y1m1.dif = cbind(y1m1.dif, y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref[, b])
                    }
                    if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y1m0.dif = cbind(y1m0.dif, y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref[, b])
                    }
                    if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y0m0.dif = cbind(y0m0.dif, y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref[, b])
                    }
                    if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                      y0m1.dif = cbind(y0m1.dif, y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref[, b])
                    }
                  } else {
                    y1m1.dif = cbind(y1m1.dif, y1m1(predict.m.data.comp, predict.y.data.comp) - y1m1.ref[, b])
                    y1m0.dif = cbind(y1m0.dif, y1m0(predict.m.data.comp, predict.y.data.comp) - y1m0.ref[, b])
                    y0m0.dif = cbind(y0m0.dif, y0m0(predict.m.data.comp, predict.y.data.comp) - y0m0.ref[, b])
                    y0m1.dif = cbind(y0m1.dif, y0m1(predict.m.data.comp, predict.y.data.comp) - y0m1.ref[, b])
                  }
                }
              }
            }
          }
        }
        
        if(is.null(plot.effect)){
          if(is.U){
            sens.results = NULL
            if("TIE" %in% sens.effect){
              TIE = apply(est.y1m1 - est.y1m0, 2, mean)
              sens.results = cbind(sens.results, TIE)
            }
            if("PDE" %in% sens.effect){
              PDE = apply(est.y1m0 - est.y0m0, 2, mean)
              sens.results = cbind(sens.results, PDE)
            }
            if("PIE" %in% sens.effect){
              PIE = apply(est.y0m1 - est.y0m0, 2, mean)
              sens.results = cbind(sens.results, PIE)
            }
            if("TDE" %in% sens.effect){
              TDE = apply(est.y1m1 - est.y0m1, 2, mean)
              sens.results = cbind(sens.results, TDE)
            }
            if("INT" %in% sens.effect){
              INT = apply(est.y1m1 - est.y1m0, 2, mean) - apply(est.y0m1 - est.y0m0, 2, mean)
              sens.results = cbind(sens.results, INT)
            }
          } else {
            TIE = apply(est.y1m1 - est.y1m0, 2, mean)
            PDE = apply(est.y1m0 - est.y0m0, 2, mean)
            PIE = apply(est.y0m1 - est.y0m0, 2, mean)
            INT = TIE - PIE
            TDE = PDE + INT
            results = cbind(TIE, PIE, PDE, TDE, INT)
          }
        }
        
        if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
          if(is.null(confounders)|
             m.scale == "continuous" & y.scale == "continuous"|
             m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
            if(!is.null(ref.mod.values)){
              if(is.null(plot.effect)){
                if(is.U){
                  if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m1.ref = apply(y1m1.ref, 2, mean)
                  }
                  if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m0.ref = apply(y1m0.ref, 2, mean)
                  }
                  if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m0.ref = apply(y0m0.ref, 2, mean)
                  }
                  if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m1.ref = apply(y0m1.ref, 2, mean)
                  }
                } else {
                  TIE.ref = apply(y1m1.ref - y1m0.ref, 2, mean)
                  PDE.ref = apply(y1m0.ref - y0m0.ref, 2, mean)
                  PIE.ref = apply(y0m1.ref - y0m0.ref, 2, mean)
                  INT.ref = TIE.ref - PIE.ref
                  TDE.ref = PDE.ref + INT.ref
                }
              } else {
                if(plot.effect == "TE"){
                  results = apply(y1m1.ref - y0m0.ref, 2, mean)
                }
                if(plot.effect == "TIE"){
                  results = apply(y1m1.ref - y1m0.ref, 2, mean)
                }
                if(plot.effect == "PIE"){
                  results = apply(y0m1.ref - y0m0.ref, 2, mean)
                }
                if(plot.effect == "TDE"){
                  results = apply(y1m1.ref - y0m1.ref, 2, mean)
                }
                if(plot.effect == "PDE"){
                  results = apply(y1m0.ref - y0m0.ref, 2, mean)
                }
                if(plot.effect == "INT"){
                  results = apply((y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref), 2, mean)
                }
              }
            }  
            if(!is.null(comp.mod.values)){
              if(is.null(plot.effect)){
                if(is.U){
                  if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m1.dif = apply(y1m1.dif, 2, mean)
                  }
                  if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y1m0.dif = apply(y1m0.dif, 2, mean)
                  }
                  if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m0.dif = apply(y0m0.dif, 2, mean)
                  }
                  if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                    y0m1.dif = apply(y0m1.dif, 2, mean)
                  }
                } else {
                  TIE.dif = apply(y1m1.dif - y1m0.dif, 2, mean)
                  PDE.dif = apply(y1m0.dif - y0m0.dif, 2, mean)
                  PIE.dif = apply(y0m1.dif - y0m0.dif, 2, mean)
                  INT.dif = TIE.dif - PIE.dif
                  TDE.dif = PDE.dif + INT.dif
                }
              }
            }
          }
        }
      }
      
      eff.earth.fun = function(effect = NULL, plot.effect = NULL, sens.effect = NULL){
        if(!is.null(ref.mod.values)){
          if(any(is.na(ref.mod.values))){
            earth.moderators = colnames(model.matrix(as.formula(paste0("~", paste0(moderators[-which(is.na(ref.mod.values))], collapse = "+"))), data = data))[-1]
          } else {
            earth.moderators = unique(c(m.moderators.new[[1]], y.moderators.new[[1]]))
          }
          earth.moderators = gsub("[^[:alnum:]\\:\\s]", "", earth.moderators)
          earth.data = cbind(predict.m.data, predict.y.data)[, earth.moderators]
          earth.ref = cbind(predict.m.data.ref, predict.y.data.ref)[1, earth.moderators]
          ref.fun = function(x){
            smooth = earth(earth.data, x)
            ref = as.numeric(predict(smooth, earth.ref))
            return(ref)
          }
          if(is.null(plot.effect)){
            if(is.null(effect)){
              if("TIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m1.ref = apply(est.y1m1, 2, ref.fun)
              }
              if("TIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m0.ref = apply(est.y1m0, 2, ref.fun)
              }
              if("PIE.ref" %in% sens.effect|"PDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m0.ref = apply(est.y0m0, 2, ref.fun)
              }
              if("PIE.ref" %in% sens.effect|"TDE.ref" %in% sens.effect|"INT.ref" %in% sens.effect|"PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m1.ref = apply(est.y0m1, 2, ref.fun)
              }
              results = NULL
              if("TIE.ref" %in% sens.effect){
                TIE.ref = y1m1.ref - y1m0.ref
                results = cbind(results, TIE.ref)
              }
              if("PDE.ref" %in% sens.effect){
                PDE.ref = y1m0.ref - y0m0.ref
                results = cbind(results, PDE.ref)
              }
              if("PIE.ref" %in% sens.effect){
                PIE.ref = y0m1.ref - y0m0.ref
                results = cbind(results, PIE.ref)
              }
              if("TDE.ref" %in% sens.effect){
                TDE.ref = y1m1.ref - y0m1.ref
                results = cbind(results, TDE.ref)
              }
              if("INT.ref" %in% sens.effect){
                INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
                results = cbind(results, INT.ref)
              }
            } else {
              y1m1.ref = apply(est.y1m1, 2, ref.fun)
              y1m0.ref = apply(est.y1m0, 2, ref.fun)
              y0m0.ref = apply(est.y0m0, 2, ref.fun)
              y0m1.ref = apply(est.y0m1, 2, ref.fun)
              TIE.ref = y1m1.ref - y1m0.ref
              PDE.ref = y1m0.ref - y0m0.ref
              PIE.ref = y0m1.ref - y0m0.ref
              INT.ref = TIE.ref - PIE.ref
              TDE.ref = PDE.ref + INT.ref
              results = cbind(TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref)
            }
          } else {
            if(plot.effect == "TE"){
              y1m1.ref = apply(est.y1m1, 2, ref.fun)
              y0m0.ref = apply(est.y0m0, 2, ref.fun)
              results = y1m1.ref - y0m0.ref
            }
            if(plot.effect == "TIE"){
              y1m1.ref = apply(est.y1m1, 2, ref.fun)
              y1m0.ref = apply(est.y1m0, 2, ref.fun)
              results = y1m1.ref - y1m0.ref
            }
            if(plot.effect == "PIE"){
              y0m0.ref = apply(est.y0m0, 2, ref.fun)
              y0m1.ref = apply(est.y0m1, 2, ref.fun)
              results = y0m1.ref - y0m0.ref
            }
            if(plot.effect == "TDE"){
              y1m1.ref = apply(est.y1m1, 2, ref.fun)
              y0m1.ref = apply(est.y0m1, 2, ref.fun)
              results = y1m1.ref - y0m1.ref
            }
            if(plot.effect == "PDE"){
              y1m0.ref = apply(est.y1m0, 2, ref.fun)
              y0m0.ref = apply(est.y0m0, 2, ref.fun)
              results = y1m0.ref - y0m0.ref
            }
            if(plot.effect == "INT"){
              y1m1.ref = apply(est.y1m1, 2, ref.fun)
              y1m0.ref = apply(est.y1m0, 2, ref.fun)
              y0m0.ref = apply(est.y0m0, 2, ref.fun)
              y0m1.ref = apply(est.y0m1, 2, ref.fun)
              results = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
            }
            return(results)
          }
        }
        
        if(is.null(plot.effect)){
          if(!is.null(comp.mod.values)){
            if(any(is.na(comp.mod.values))){
              earth.moderators = colnames(model.matrix(as.formula(paste0("~", paste0(moderators[-which(is.na(comp.mod.values))], collapse = "+"))), data = data))[-1]
            } else {
              earth.moderators = unique(c(m.moderators.new[[1]], y.moderators.new[[1]]))
            }
            earth.moderators = gsub("[^[:alnum:]\\:\\s]", "", earth.moderators)
            earth.data = cbind(predict.m.data, predict.y.data)[, earth.moderators]
            earth.comp = cbind(predict.m.data.comp, predict.y.data.comp)[1, earth.moderators]
            comp.fun = function(x){
              smooth = earth(earth.data, x)
              comp = as.numeric(predict(smooth, earth.comp))
              return(comp)
            }
            if(is.null(effect)){
              if("TIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m1.dif = apply(est.y1m1, 2, comp.fun) - y1m1.ref
              }
              if("TIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y1m0.dif = apply(est.y1m0, 2, comp.fun) - y1m0.ref
              }
              if("PIE.dif" %in% sens.effect|"PDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m0.dif = apply(est.y0m0, 2, comp.fun) - y0m0.ref
              }
              if("PIE.dif" %in% sens.effect|"TDE.dif" %in% sens.effect|"INT.dif" %in% sens.effect){
                y0m1.dif = apply(est.y0m1, 2, comp.fun) - y0m1.ref
              }
              if("TIE.dif" %in% sens.effect){
                TIE.dif = y1m1.dif - y1m0.dif
                results = cbind(results, TIE.dif)
              }
              if("PDE.dif" %in% sens.effect){
                PDE.dif = y1m0.dif - y0m0.dif
                results = cbind(results, PDE.dif)
              }
              if("PIE.dif" %in% sens.effect){
                PIE.dif = y0m1.dif - y0m0.dif
                results = cbind(results, PIE.dif)
              }
              if("TDE.dif" %in% sens.effect){
                TDE.dif = y1m1.dif - y0m1.dif
                results = cbind(results, TDE.dif)
              }
              if("INT.dif" %in% sens.effect){
                INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
                results = cbind(results, INT.dif)
              }
            } else {
              y1m1.dif = apply(est.y1m1, 2, comp.fun) - y1m1.ref
              y1m0.dif = apply(est.y1m0, 2, comp.fun) - y1m0.ref
              y0m0.dif = apply(est.y0m0, 2, comp.fun) - y0m0.ref
              y0m1.dif = apply(est.y0m1, 2, comp.fun) - y0m1.ref
              TIE.dif = y1m1.dif - y1m0.dif
              PDE.dif = y1m0.dif - y0m0.dif
              PIE.dif = y0m1.dif - y0m0.dif
              INT.dif = TIE.dif - PIE.dif
              TDE.dif = PDE.dif + INT.dif
              results = cbind(results, TIE.dif, PIE.dif, PDE.dif, TDE.dif, INT.dif)
            }
          }
        }
        if(!is.null(effect))
          return(results[, effect])
        if(!is.null(sens.effect))
          return(results)
      }
      
      if(!any(is.na(c(ref.mod.disc.values, ref.mod.cont.values)))){
        if(is.null(confounders)){
          if(is.U){
            if("TIE.ref" %in% sens.effect){
              TIE.ref = y1m1.ref - y1m0.ref
              sens.results = cbind(sens.results, TIE.ref)
            }
            if("PDE.ref" %in% sens.effect){
              PDE.ref = y1m0.ref - y0m0.ref
              sens.results = cbind(sens.results, PDE.ref)
            }
            if("PIE.ref" %in% sens.effect){
              PIE.ref = y0m1.ref - y0m0.ref
              sens.results = cbind(sens.results, PIE.ref)
            }
            if("TDE.ref" %in% sens.effect){
              TDE.ref = y1m1.ref - y0m1.ref
              sens.results = cbind(sens.results, TDE.ref)
            }
            if("INT.ref" %in% sens.effect){
              INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
              sens.results = cbind(sens.results, INT.ref)
            }
            if("TIE.dif" %in% sens.effect){
              TIE.dif = y1m1.dif - y1m0.dif
              sens.results = cbind(sens.results, TIE.dif)
            }
            if("PDE.dif" %in% sens.effect){
              PDE.dif = y1m0.dif - y0m0.dif
              sens.results = cbind(sens.results, PDE.dif)
            }
            if("PIE.dif" %in% sens.effect){
              PIE.dif = y0m1.dif - y0m0.dif
              sens.results = cbind(sens.results, PIE.dif)
            }
            if("TDE.dif" %in% sens.effect){
              TDE.dif = y1m1.dif - y0m1.dif
              sens.results = cbind(sens.results, TDE.dif)
            }
            if("INT.dif" %in% sens.effect){
              INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
              sens.results = cbind(sens.results, INT.dif)
            }
          } else {
            if(!is.null(ref.mod.values)){
              if(!is.null(comp.mod.values)){
                results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref, TIE.dif, PIE.dif, PDE.dif, TDE.dif, INT.dif)
              } else {
                if(is.null(plot.effect)){
                  results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref)
                } else {
                  return(results)
                }
              }
            }
          }
        } else if(m.scale == "continuous" & y.scale == "continuous"){
          if(!paste0(treatment, ":", mediator) %in% y.predictors){
            if(is.U){
              if("TIE.ref" %in% sens.effect){
                TIE.ref = y1m1.ref - y1m0.ref
                sens.results = cbind(sens.results, TIE.ref)
              }
              if("PDE.ref" %in% sens.effect){
                PDE.ref = y1m0.ref - y0m0.ref
                sens.results = cbind(sens.results, PDE.ref)
              }
              if("PIE.ref" %in% sens.effect){
                PIE.ref = y0m1.ref - y0m0.ref
                sens.results = cbind(sens.results, PIE.ref)
              }
              if("TDE.ref" %in% sens.effect){
                TDE.ref = y1m1.ref - y0m1.ref
                sens.results = cbind(sens.results, TDE.ref)
              }
              if("INT.ref" %in% sens.effect){
                INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
                sens.results = cbind(sens.results, INT.ref)
              }
              if("TIE.dif" %in% sens.effect){
                TIE.dif = y1m1.dif - y1m0.dif
                sens.results = cbind(sens.results, TIE.dif)
              }
              if("PDE.dif" %in% sens.effect){
                PDE.dif = y1m0.dif - y0m0.dif
                sens.results = cbind(sens.results, PDE.dif)
              }
              if("PIE.dif" %in% sens.effect){
                PIE.dif = y0m1.dif - y0m0.dif
                sens.results = cbind(sens.results, PIE.dif)
              }
              if("TDE.dif" %in% sens.effect){
                TDE.dif = y1m1.dif - y0m1.dif
                sens.results = cbind(sens.results, TDE.dif)
              }
              if("INT.dif" %in% sens.effect){
                INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
                sens.results = cbind(sens.results, INT.dif)
              }
            } else {
              if(!is.null(ref.mod.values)){
                if(!is.null(comp.mod.values)){
                  results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref, TIE.dif, PIE.dif, PDE.dif, TDE.dif, INT.dif)
                } else {
                  if(is.null(plot.effect)){
                    results = cbind(results, TIE.ref, PIE.ref, PDE.ref, TDE.ref, INT.ref)
                  } else {
                    return(results)
                  }
                }
              }
            }
          } else {
            if(is.U){
              if("TIE.ref" %in% sens.effect){
                TIE.ref = y1m1.ref - y1m0.ref
                sens.results = cbind(sens.results, TIE.ref)
              }
              if("PIE.ref" %in% sens.effect){
                PIE.ref = y0m1.ref - y0m0.ref
                sens.results = cbind(sens.results, PIE.ref)
              }
              if("INT.ref" %in% sens.effect){
                INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
                sens.results = cbind(sens.results, INT.ref)
              }
              if("TIE.dif" %in% sens.effect){
                TIE.dif = y1m1.dif - y1m0.dif
                sens.results = cbind(sens.results, TIE.dif)
              }
              if("PIE.dif" %in% sens.effect){
                PIE.dif = y0m1.dif - y0m0.dif
                sens.results = cbind(sens.results, PIE.dif)
              }
              if("INT.dif" %in% sens.effect){
                INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
                sens.results = cbind(sens.results, INT.dif)
              }
              if(any(sens.effect %in% c("TDE.ref", "PDE.ref", "TDE.dif", "PDE.dif"))){
                if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TIE.ref", "PIE.ref", "INT.ref", "TIE.dif", "PIE.dif", "INT.dif"))){
                  sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TIE.ref", "PIE.ref", "INT.ref", "TIE.dif", "PIE.dif", "INT.dif"))]))
                } else {
                  sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
                }
              }
            } else {
              if(!is.null(ref.mod.values)){
                if(!is.null(comp.mod.values)){
                  results = cbind(results, TIE.ref, PIE.ref, INT.ref, TIE.dif, PIE.dif, INT.dif, eff.earth.fun(effect = c("PDE.ref", "TDE.ref", "PDE.dif", "TDE.dif")))
                } else {
                  if(is.null(plot.effect)){
                    results = cbind(results, TIE.ref, PIE.ref, INT.ref, eff.earth.fun(effect = c("PDE.ref", "TDE.ref")))
                  } else {
                    if(plot.effect == "TIE"|plot.effect == "PIE"|plot.effect == "INT"){
                      return(results)
                    } else {
                      return(eff.earth.fun(plot.effect = plot.effect))
                    }
                  }
                }
              }
            }
          }
        } else if(m.scale == "binary" & y.scale == "continuous" & !paste0(treatment, ":", mediator) %in% y.predictors){
          if(is.U){
            if("PDE.ref" %in% sens.effect){
              PDE.ref = y1m0.ref - y0m0.ref
              sens.results = cbind(sens.results, PDE.ref)
            }
            if("TDE.ref" %in% sens.effect){
              TDE.ref = y1m1.ref - y0m1.ref
              sens.results = cbind(sens.results, TDE.ref)
            }
            if("INT.ref" %in% sens.effect){
              INT.ref = (y1m1.ref - y1m0.ref) - (y0m1.ref - y0m0.ref)
              sens.results = cbind(sens.results, INT.ref)
            }
            if("PDE.dif" %in% sens.effect){
              PDE.dif = y1m0.dif - y0m0.dif
              sens.results = cbind(sens.results, PDE.dif)
            }
            if("TDE.dif" %in% sens.effect){
              TDE.dif = y1m1.dif - y0m1.dif
              sens.results = cbind(sens.results, TDE.dif)
            }
            if("INT.dif" %in% sens.effect){
              INT.dif = (y1m1.dif - y1m0.dif) - (y0m1.dif - y0m0.dif)
              sens.results = cbind(sens.results, INT.dif)
            }
            if(any(sens.effect %in% c("TIE.ref", "PIE.ref", "TIE.dif", "PIE.dif"))){
              if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TDE.ref", "PDE.ref", "INT.ref", "TDE.dif", "PDE.dif", "INT.dif"))){
                sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT", "TDE.ref", "PDE.ref", "INT.ref", "TDE.dif", "PDE.dif", "INT.dif"))]))
              } else {
                sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
              }
            }
          } else {
            if(!is.null(ref.mod.values)){
              if(!is.null(comp.mod.values)){
                results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "TIE.dif", "PIE.dif")), PDE.ref, TDE.ref, INT.ref, PDE.dif, TDE.dif, INT.dif)
              } else {
                if(is.null(plot.effect)){
                  results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref")), PDE.ref, TDE.ref, INT.ref)
                } else {
                  if(plot.effect == "TDE"|plot.effect == "PDE"|plot.effect == "INT"){
                    return(results)
                  } else {
                    return(eff.earth.fun(plot.effect = plot.effect))
                  }
                }
              }
            }
          }
        } else {
          if(is.U){
            if(any(sens.effect %in% c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif"))){
              if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))){
                sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))]))
              } else {
                sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
              }
            }
          } else {
            if(!is.null(ref.mod.values)){
              if(!is.null(comp.mod.values)){
                results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif")))
              } else {
                if(is.null(plot.effect)){
                  results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref")))
                } else {
                  return(eff.earth.fun(plot.effect = plot.effect))
                }
              }
            }
          }
        }
      } else {
        if(is.U){
          if(any(sens.effect %in% c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif"))){
            if(any(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))){
              sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect[-which(sens.effect %in% c("TIE", "PIE", "TDE", "PDE", "INT"))]))
            } else {
              sens.results = cbind(sens.results, eff.earth.fun(sens.effect = sens.effect))
            }
          }
        } else {
          if(!is.null(ref.mod.values)){
            if(!is.null(comp.mod.values)){
              results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", "INT.dif")))
            } else {
              if(is.null(plot.effect)){
                results = cbind(results, eff.earth.fun(effect = c("TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref")))
              } else {
                return(eff.earth.fun(plot.effect = plot.effect))
              }
            }
          }
        }
      }
      
      if(is.null(plot.effect)){
        if(is.U){
          return(sens.results[, sens.effect])
        } else {
          # Step 5. Compute summary statistics such as point estimates and confidence intervals.
          results.ori = results
          TE = results.ori[, "TIE"] + results.ori[, "PDE"]
          results = cbind(TE, results.ori[, 1:5])
          if(!is.null(ref.mod.values)){
            TE.ref = results.ori[, "TIE.ref"] + results.ori[, "PDE.ref"]
            results = cbind(results, TE.ref, results.ori[, 6:10])
            if(!is.null(comp.mod.values)){
              TE.dif = results.ori[, "TIE.dif"] + results.ori[, "PDE.dif"]
              results = cbind(results, TE.dif, results.ori[, 11:15])
            }
          }
          
          est.results = apply(results, 2, mean)
          se.results = apply(results, 2, sd)
          ci.results = apply(results, 2, quantile, probs = c((1 - conf.level)/2, (1 + conf.level)/2))
          
          
          # Output
          # Mediator and outcome models
          summary.table = function(moderators, predictors, l){
            if(is.null(moderators)){
              return(summary(l))
            }
            
            if(!is.null(moderators)){
              model.coef = names(coef(l))
              
              if(all(predictors %in% y.predictors.new) & all(y.predictors.new %in% predictors)){
                if(paste0(treatment, ":", mediator) %in% model.coef)
                  predictors[which(predictors == paste0(treatment, ":", mediator))] = paste0(treatment, ":", mediator)
                if(paste0(mediator, ":", treatment) %in% model.coef)
                  predictors[which(predictors == paste0(treatment, ":", mediator))] = paste0(mediator, ":", treatment)
              }
              
              main = rep("Intercept", 1 + length(moderators[[1]]))
              moderation = c("Intercept", moderators[[1]])
              for(i in 1:length(predictors)){
                moderation = c(moderation, "Intercept", moderators[[i + 1]])
                main = c(main, rep(predictors[i], 1 + length(moderators[[i + 1]])))
              }
              
              coef = NULL
              for(i in 1:length(main)){
                if(main[i] == "Intercept"){
                  coef = c(coef, paste0(moderation[i]))
                } else {
                  if(moderation[i] == "Intercept"){
                    coef = c(coef, paste0(main[i]))
                  } else {
                    main.i = unlist(strsplit(main[i], ":"))
                    if(length(main.i) == 1){
                      if(paste0(main[i], ":", moderation[i]) %in% model.coef)
                        coef = c(coef, paste0(main[i], ":", moderation[i]))
                      if(paste0(moderation[i], ":", main[i]) %in% model.coef)
                        coef = c(coef, paste0(moderation[i], ":", main[i]))
                    }
                    if(length(main.i) == 2){
                      if(paste0(main.i[1], ":", main.i[2], ":", moderation[i]) %in% model.coef)
                        coef = c(coef, paste0(main.i[1], ":", main.i[2], ":", moderation[i]))
                      if(paste0(main.i[1], ":", moderation[i], ":", main.i[2]) %in% model.coef)
                        coef = c(coef, paste0(main.i[1], ":", moderation[i], ":", main.i[2]))
                      if(paste0(main.i[2], ":", main.i[1], ":", moderation[i]) %in% model.coef)
                        coef = c(coef, paste0(main.i[2], ":", main.i[1], ":", moderation[i]))
                      if(paste0(main.i[2], ":", moderation[i], ":", main.i[1]) %in% model.coef)
                        coef = c(coef, paste0(main.i[2], ":", moderation[i], ":", main.i[1]))
                      if(paste0(moderation[i], ":", main.i[1], ":", main.i[2]) %in% model.coef)
                        coef = c(coef, paste0(moderation[i], ":", main.i[1], ":", main.i[2]))
                      if(paste0(moderation[i], ":", main.i[2], ":", main.i[1]) %in% model.coef)
                        coef = c(coef, paste0(moderation[i], ":", main.i[2], ":", main.i[1]))
                    }
                  }
                }
              }
              
              main = c("Intercept", rep("", length(moderators[[1]])))
              for(i in 1:length(predictors)){
                main = c(main, predictors[i], rep("", length(moderators[[i + 1]])))
              }
              
              coef[which(coef == "Intercept")] = "(Intercept)"
              summary.table = as.data.frame(matrix(NA, length(main), 6))
              summary.table[, 1] = main
              summary.table[, 2] = moderation
              summary.table[, 3:6] = as.data.frame(summary(l)$coefficients[coef, ])
              mod_summary_sign = summary.table[, 6]  
              mod_summary_stars = NA                             # Named vector with significance stars
              mod_summary_stars[mod_summary_sign < 0.1] = "."
              mod_summary_stars[mod_summary_sign < 0.05] = "*"
              mod_summary_stars[mod_summary_sign < 0.01] = "**"
              mod_summary_stars[mod_summary_sign < 0.001] = "***"
              mod_summary_stars[is.na(mod_summary_stars)] = ""
              summary.table = cbind(summary.table, mod_summary_stars)
              colnames(summary.table) = c("main", "moderation", colnames(summary(l)$coefficients), "")
              rownames(summary.table) = coef
              return(summary.table)
            }
          }
          
          summary.m = summary.table(m.moderators.new, m.predictors.new, l.m)
          summary.y = summary.table(y.moderators.new, y.predictors.new, l.y)
          
          # Causal effects
          summary.effects = cbind(est.results, se.results, ci.results[1, ], ci.results[2, ])
          
          colnames(summary.effects) = c("Estimate", "Std. Error", paste(conf.level * 100, "% CI Lower ", (1-conf.level)/2 *100, "%", sep=""),
                                        paste(conf.level * 100, "% CI Upper ", (1-conf.level)/2 *100, "%", sep=""))
          
          args = as.list(match.call())[2:length(as.list(match.call()))]
          args.full = mget(names(formals()), sys.frame(sys.nframe()))
          
          if(y.scale == "continuous")
            return(list(effects = summary.effects, m.model = summary.m, y.model = summary.y, results = results, args = args, args.full = args.full, l.m = l.m, l.y = l.y, formula.m = formula.m, formula.y = formula.y, m.predictors = m.predictors, y.predictors = y.predictors, m.moderators = m.moderators, y.moderators = y.moderators, m.predictors.new = m.predictors.new, y.predictors.new = y.predictors.new, m.moderators.new = m.moderators.new, y.moderators.new = y.moderators.new, predict.m.data.ref = predict.m.data.ref, predict.m.data.comp = predict.m.data.comp, predict.m.data = predict.m.data, predict.y.data.ref = predict.y.data.ref, predict.y.data.comp = predict.y.data.comp, predict.y.data = predict.y.data, data = data))
          if(y.scale == "binary")
            return(list(effects = summary.effects, m.model = summary.m, y.model = summary.y, results = results, args = args, args.full = args.full, l.m = l.m, l.y = l.y, formula.m = formula.m, formula.y = formula.y, m.predictors = m.predictors, y.predictors = y.predictors, m.moderators = m.moderators, y.moderators = y.moderators, m.predictors.new = m.predictors.new, y.predictors.new = y.predictors.new, m.moderators.new = m.moderators.new, y.moderators.new = y.moderators.new, predict.m.data.ref = predict.m.data.ref, predict.m.data.comp = predict.m.data.comp, predict.m.data = predict.m.data, predict.y.data.ref = predict.y.data.ref, predict.y.data.comp = predict.y.data.comp, predict.y.data = predict.y.data, data = data))
        }
      }
    }
    
    s.total = s.total
    i = grid[ceiling(s.total/nsim), 1]
    j = grid[ceiling(s.total/nsim), 2]
    args.full = args.full.ori
    args.full$b.m = b.m.all[i]
    args.full$b.y = b.y.all[j]
    if(!is.null(b.t))
      b.t = b.t.ori 
    if(y.scale == "binary")
      sd.y = sd.y.ori/sqrt(1 - args.full$b.y^2) # When U is adjusted for, the latent index on the original scale changes by args.full$b.y * sqrt(var.y/var.u) * U, where args.full$b.y is standardized and transformed back to the original coefficient after multiplying with sqrt(var.y/var.u). Hence, var.y = var.y.ori + args.full$b.y^2 * var.y. Hence, var.y = var.y.ori/(1- args.full$b.y^2).
    if(m.scale == "binary")
      sd.m = sd.m.ori/sqrt(1 - args.full$b.m^2) # When U is adjusted for, the latent index on the original scale changes by args.full$b.m * sqrt(var.m/var.u) * U, where args.full$b.m is standardized and transformed back to the original coefficient after multiplying with sqrt(var.m/var.u). Hence, var.m = var.m.ori + args.full$b.m^2 * var.m. Hence, var.m = var.m.ori/(1- args.full$b.m^2).
    if(t.rand == TRUE){
      if(U.scale == "binary"){
        args.full$b.y = args.full$b.y/sqrt(p.u * (1 - p.u)) * sd.y
        args.full$b.m = args.full$b.m/sqrt(p.u * (1 - p.u)) * sd.m
        args.full$data$Unmeasure = genU(b.y = args.full$b.y, b.m = args.full$b.m)$Unmeasure # In genU, everything is on the original scale, and thus the sensitivity parameter is transformed back to the original scale.
      }
      if(U.scale == "continuous"){
        args.full$b.y = args.full$b.y/sigma.u * sd.y
        args.full$b.m = args.full$b.m/sigma.u * sd.m
        args.full$data$Unmeasure = genU(b.y = args.full$b.y, b.m = args.full$b.m)$Unmeasure # In genU, everything is on the original scale, and thus the sensitivity parameter is transformed back to the original scale.
      }
    }
    if(t.rand == FALSE){
      if(t.scale == "continuous")
        sd.t = sd(data[, t])
      if(t.scale == "binary"){
        l.t = glm(as.formula(paste(t, "~", formula.t)), family = binomial(link = "probit"), data = data)
        sd.t = sqrt(var(predict(l.t, type = "link")) + 1)/sqrt(1 - b.t^2)
      }
      if(U.scale == "binary"){
        args.full$b.y = args.full$b.y/sqrt(p.u * (1 - p.u)) * sd.y
        args.full$b.m = args.full$b.m/sqrt(p.u * (1 - p.u)) * sd.m
        b.t = b.t/sqrt(p.u * (1 - p.u)) * sd.t
        args.full$data$Unmeasure = genU(b.y = args.full$b.y, b.m = args.full$b.m, b.t = b.t)$Unmeasure# In genU, everything is on the original scale, and thus the sensitivity parameter is transformed back to the original scale.
      }
      if(U.scale == "continuous"){
        args.full$b.y = args.full$b.y/sigma.u * sd.y
        args.full$b.m = args.full$b.m/sigma.u * sd.m
        b.t = b.t/sigma.u * sd.t
        args.full$data$Unmeasure = genU(b.y = args.full$b.y, b.m = args.full$b.m, b.t = b.t)$Unmeasure# In genU, everything is on the original scale, and thus the sensitivity parameter is transformed back to the original scale.
      }
    }
    return(do.call(modmed, args.full))
  }
  
  if(interactive()){
    close(pb)
  }
  
  vals = matrix(NA, grid.b.m, grid.b.y)
  rownames(vals) = round(b.m.all, 2)
  colnames(vals) = round(b.y.all, 2)
  if("TIE" %in% sens.effect){
    TIE.all = CIL_TIE.all = CIU_TIE.all = vals
  }
  if("PDE" %in% sens.effect){
    PDE.all = CIL_PDE.all = CIU_PDE.all = vals
  }
  if("PIE" %in% sens.effect){
    PIE.all = CIL_PIE.all = CIU_PIE.all = vals
  }
  if("TDE" %in% sens.effect){
    TDE.all = CIL_TDE.all = CIU_TDE.all = vals
  }
  if("INT" %in% sens.effect){
    INT.all = CIL_INT.all = CIU_INT.all = vals
  }
  if("TIE.ref" %in% sens.effect){
    TIE.ref.all = CIL_TIE.ref.all = CIU_TIE.ref.all = vals
  }
  if("PDE.ref" %in% sens.effect){
    PDE.ref.all = CIL_PDE.ref.all = CIU_PDE.ref.all = vals
  }
  if("PIE.ref" %in% sens.effect){
    PIE.ref.all = CIL_PIE.ref.all = CIU_PIE.ref.all = vals
  }
  if("TDE.ref" %in% sens.effect){
    TDE.ref.all = CIL_TDE.ref.all = CIU_TDE.ref.all = vals
  }
  if("INT.ref" %in% sens.effect){
    INT.ref.all = CIL_INT.ref.all = CIU_INT.ref.all = vals
  }
  if("TIE.dif" %in% sens.effect){
    TIE.dif.all = CIL_TIE.dif.all = CIU_TIE.dif.all = vals
  }
  if("PDE.dif" %in% sens.effect){
    PDE.dif.all = CIL_PDE.dif.all = CIU_PDE.dif.all = vals
  }
  if("PIE.dif" %in% sens.effect){
    PIE.dif.all = CIL_PIE.dif.all = CIU_PIE.dif.all = vals
  }
  if("TDE.dif" %in% sens.effect){
    TDE.dif.all = CIL_TDE.dif.all = CIU_TDE.dif.all = vals
  }
  if("INT.dif" %in% sens.effect){
    INT.dif.all = CIL_INT.dif.all = CIU_INT.dif.all = vals
  }
  
  for(s.grid in 1:nrow(grid)){
    results = matrix(NA, nsim * args.full$nmc, length(sens.effect))
    for(s.nsim in 1:nsim){
      results[(args.full$nmc * (s.nsim - 1) + 1):(args.full$nmc * s.nsim), ] = resij[[(s.grid - 1) * nsim + s.nsim]]
    }
    colnames(results) = sens.effect
    results.s = c(apply(results, 2, mean), CIL = apply(results, 2, quantile, probs = c((1 - conf.level)/2)), CIU = apply(results, 2, quantile, probs = c((1 + conf.level)/2)))
    
    i = grid[s.grid, 1]
    j = grid[s.grid, 2]
    
    if("TIE" %in% sens.effect){
      TIE.all[i, j] = results.s["TIE"]
      CIL_TIE.all[i, j] = results.s["CIL.TIE"]
      CIU_TIE.all[i, j] = results.s["CIU.TIE"]
    }
    if("PIE" %in% sens.effect){
      PIE.all[i, j] = results.s["PIE"]
      CIL_PIE.all[i, j] = results.s["CIL.PIE"]
      CIU_PIE.all[i, j] = results.s["CIU.PIE"]
    }
    if("PDE" %in% sens.effect){
      PDE.all[i, j] = results.s["PDE"]
      CIL_PDE.all[i, j] = results.s["CIL.PDE"]
      CIU_PDE.all[i, j] = results.s["CIU.PDE"]
    }
    if("TDE" %in% sens.effect){
      TDE.all[i, j] = results.s["TDE"]
      CIL_TDE.all[i, j] = results.s["CIL.TDE"]
      CIU_TDE.all[i, j] = results.s["CIU.TDE"]
    }
    if("INT" %in% sens.effect){
      INT.all[i, j] = results.s["INT"]
      CIL_INT.all[i, j] = results.s["CIL.INT"]
      CIU_INT.all[i, j] = results.s["CIU.INT"]
    }
    if("TIE.ref" %in% sens.effect){
      TIE.ref.all[i, j] = results.s["TIE.ref"]
      CIL_TIE.ref.all[i, j] = results.s["CIL.TIE.ref"]
      CIU_TIE.ref.all[i, j] = results.s["CIU.TIE.ref"]
    }
    if("PIE.ref" %in% sens.effect){
      PIE.ref.all[i, j] = results.s["PIE.ref"]
      CIL_PIE.ref.all[i, j] = results.s["CIL.PIE.ref"]
      CIU_PIE.ref.all[i, j] = results.s["CIU.PIE.ref"]
    }
    if("PDE.ref" %in% sens.effect){
      PDE.ref.all[i, j] = results.s["PDE.ref"]
      CIL_PDE.ref.all[i, j] = results.s["CIL.PDE.ref"]
      CIU_PDE.ref.all[i, j] = results.s["CIU.PDE.ref"]
    }
    if("TDE.ref" %in% sens.effect){
      TDE.ref.all[i, j] = results.s["TDE.ref"]
      CIL_TDE.ref.all[i, j] = results.s["CIL.TDE.ref"]
      CIU_TDE.ref.all[i, j] = results.s["CIU.TDE.ref"]
    }
    if("INT.ref" %in% sens.effect){
      INT.ref.all[i, j] = results.s["INT.ref"]
      CIL_INT.ref.all[i, j] = results.s["CIL.INT.ref"]
      CIU_INT.ref.all[i, j] = results.s["CIU.INT.ref"]
    }
    if("TIE.dif" %in% sens.effect){
      TIE.dif.all[i, j] = results.s["TIE.dif"]
      CIL_TIE.dif.all[i, j] = results.s["CIL.TIE.dif"]
      CIU_TIE.dif.all[i, j] = results.s["CIU.TIE.dif"]
    }
    if("PIE.dif" %in% sens.effect){
      PIE.dif.all[i, j] = results.s["PIE.dif"]
      CIL_PIE.dif.all[i, j] = results.s["CIL.PIE.dif"]
      CIU_PIE.dif.all[i, j] = results.s["CIU.PIE.dif"]
    }
    if("PDE.dif" %in% sens.effect){
      PDE.dif.all[i, j] = results.s["PDE.dif"]
      CIL_PDE.dif.all[i, j] = results.s["CIL.PDE.dif"]
      CIU_PDE.dif.all[i, j] = results.s["CIU.PDE.dif"]
    }
    if("TDE.dif" %in% sens.effect){
      TDE.dif.all[i, j] = results.s["TDE.dif"]
      CIL_TDE.dif.all[i, j] = results.s["CIL.TDE.dif"]
      CIU_TDE.dif.all[i, j] = results.s["CIU.TDE.dif"]
    }
    if("INT.dif" %in% sens.effect){
      INT.dif.all[i, j] = results.s["INT.dif"]
      CIL_INT.dif.all[i, j] = results.s["CIL.INT.dif"]
      CIU_INT.dif.all[i, j] = results.s["CIU.INT.dif"]
    }
  }
  
  results.new = NULL
  if("TIE" %in% sens.effect){
    results.new[["TIE"]] = TIE.all
    results.new[["CIL.TIE"]] = CIL_TIE.all
    results.new[["CIU.TIE"]] = CIU_TIE.all
  }
  if("PIE" %in% sens.effect){
    results.new[["PIE"]] = PIE.all
    results.new[["CIL.PIE"]] = CIL_PIE.all
    results.new[["CIU.PIE"]] = CIU_PIE.all
  }
  if("PDE" %in% sens.effect){
    results.new[["PDE"]] = PDE.all
    results.new[["CIL.PDE"]] = CIL_PDE.all
    results.new[["CIU.PDE"]] = CIU_PDE.all
  }
  if("TDE" %in% sens.effect){
    results.new[["TDE"]] = TDE.all
    results.new[["CIL.TDE"]] = CIL_TDE.all
    results.new[["CIU.TDE"]] = CIU_TDE.all
  }
  if("INT" %in% sens.effect){
    results.new[["INT"]] = INT.all
    results.new[["CIL.INT"]] = CIL_INT.all
    results.new[["CIU.INT"]] = CIU_INT.all
  }
  if("TIE.ref" %in% sens.effect){
    results.new[["TIE.ref"]] = TIE.ref.all
    results.new[["CIL.TIE.ref"]] = CIL_TIE.ref.all
    results.new[["CIU.TIE.ref"]] = CIU_TIE.ref.all
  }
  if("PIE.ref" %in% sens.effect){
    results.new[["PIE.ref"]] = PIE.ref.all
    results.new[["CIL.PIE.ref"]] = CIL_PIE.ref.all
    results.new[["CIU.PIE.ref"]] = CIU_PIE.ref.all
  }
  if("PDE.ref" %in% sens.effect){
    results.new[["PDE.ref"]] = PDE.ref.all
    results.new[["CIL.PDE.ref"]] = CIL_PDE.ref.all
    results.new[["CIU.PDE.ref"]] = CIU_PDE.ref.all
  }
  if("TDE.ref" %in% sens.effect){
    results.new[["TDE.ref"]] = TDE.ref.all
    results.new[["CIL.TDE.ref"]] = CIL_TDE.ref.all
    results.new[["CIU.TDE.ref"]] = CIU_TDE.ref.all
  }
  if("INT.ref" %in% sens.effect){
    results.new[["INT.ref"]] = INT.ref.all
    results.new[["CIL.INT.ref"]] = CIL_INT.ref.all
    results.new[["CIU.INT.ref"]] = CIU_INT.ref.all
  }
  if("TIE.dif" %in% sens.effect){
    results.new[["TIE.dif"]] = TIE.dif.all
    results.new[["CIL.TIE.dif"]] = CIL_TIE.dif.all
    results.new[["CIU.TIE.dif"]] = CIU_TIE.dif.all
  }
  if("PIE.dif" %in% sens.effect){
    results.new[["PIE.dif"]] = PIE.dif.all
    results.new[["CIL.PIE.dif"]] = CIL_PIE.dif.all
    results.new[["CIU.PIE.dif"]] = CIU_PIE.dif.all
  }
  if("PDE.dif" %in% sens.effect){
    results.new[["PDE.dif"]] = PDE.dif.all
    results.new[["CIL.PDE.dif"]] = CIL_PDE.dif.all
    results.new[["CIU.PDE.dif"]] = CIU_PDE.dif.all
  }
  if("TDE.dif" %in% sens.effect){
    results.new[["TDE.dif"]] = TDE.dif.all
    results.new[["CIL.TDE.dif"]] = CIL_TDE.dif.all
    results.new[["CIU.TDE.dif"]] = CIU_TDE.dif.all
  }
  if("INT.dif" %in% sens.effect){
    results.new[["INT.dif"]] = INT.dif.all
    results.new[["CIL.INT.dif"]] = CIL_INT.dif.all
    results.new[["CIU.INT.dif"]] = CIU_INT.dif.all
  }
  
  if(!is.null(b.t))
    b.t = b.t.ori 
  
  return(list(b.t = b.t, X.coef.plot = X.coef.plot, range.b.m = range.b.m, range.b.y = range.b.y, b.y.all = b.y.all, b.m.all = b.m.all, results.new = results.new))
}

#' Simulation-Based Sensitivity Analysis Plot for Causal Moderated Mediation Analysis
#' 
#' @param object Output from the \code{modmed} function.
#' @param sens.results An output from the \code{modmed.sens} function.
#' @param effect The name of the effect whose sensitivity analysis results are to be plotted (string). Only one effect is plotted at a time. It can be specified as "TIE", "PIE", "PDE", "TDE", "INT", "TIE.ref", "PIE.ref", "PDE.ref", "TDE.ref", "INT.ref", "TIE.dif", "PIE.dif", "PDE.dif", "TDE.dif", or "INT.dif". It must be included in sens.effect when running the \code{modmed.sens} function.
#' @return Sensitivity analysis plots for the causal effects in the causal moderated mediation analysis.
#' @author Xu Qin and Fan Yang
#' @references Qin, X., & Yang, F. (2020). Simulation-Based Sensitivity Analysis for Causal Mediation Studies.
#' @export
#' @importFrom graphics plot contour abline text 
#' @importFrom stats loess
#' @importFrom reshape2 melt
#' @examples
#' \donttest{
#' data(newws)
#' modmed.results = modmed(data = newws, treatment = "treat", mediator = "emp", outcome = "depression", covariates.disc = c("emp_prior", "nevmar", "hispanic", "nohsdip"), covariates.cont = c("workpref", "attitude", "depress_prior"), moderators.disc = "CHCNT", moderators.cont = "ADCPC", m.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), y.model = list(intercept = c("ADCPC", "CHCNT"), treatment = c("ADCPC", "CHCNT"), mediator = c("ADCPC", "CHCNT"), tm = c("ADCPC", "CHCNT"), emp_prior = NULL, nevmar = NULL, hispanic = NULL, nohsdip = NULL, workpref = NULL, attitude = NULL, depress_prior = NULL), comp.mod.disc.values = 3, ref.mod.disc.values = 2, comp.mod.cont.values = 5050, ref.mod.cont.values = 5050, m.scale = "binary", y.scale = "continuous", seed = 1) 
#' sens.results = modmed.sens(modmed.results, sens.effect = "TIE.ref", U.scale = "binary", grid.b.m = 2, grid.b.y = 2, iter = 2, nsim = 2, ncore = 1)
#' sens.plot(modmed.results, sens.results, "TIE.ref")
#' }
sens.plot = function(object, sens.results, effect){
  args.full = object$args.full
  if(!effect %in% names(sens.results$results.new))
    stop("effect must be included in sens.effect of the modmed.sens function.")
  m.scale = args.full$m.scale
  y.scale = args.full$y.scale
  
  if(m.scale == "binary")
    xlab = "Standardized partial effect of U on the latent index of M"
  if(m.scale =="continuous")
    xlab = "Standardized partial effect of U on M"
  if(y.scale == "binary")
    ylab = "Standardized partial effect of U on the latent index of Y"
  if(y.scale =="continuous")
    ylab = "Standardized partial effect of U on Y"
  
  main = paste("Sensitivity Analysis for", effect)
  
  plot(sens.results$X.coef.plot[, 1], sens.results$X.coef.plot[, 2], pch = 19, xlim = sens.results$range.b.m, ylim = sens.results$range.b.y, xlab = xlab, ylab = ylab, main = main)
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)
  text(sens.results$X.coef.plot[, 1], sens.results$X.coef.plot[, 2], 1:nrow(sens.results$X.coef.plot), pos = 4, cex = 0.6)
  text(0, 0, paste0(round(object$effects[effect, 1], 3), " \n(unadjusted)"), cex = 0.8, col ="blue")
  
  # Smooth contours using loess
  a = melt(sens.results$results.new[[effect]])
  b = loess(value ~ Var1 + Var2, data = a, span = 0.8) # higher values of span will deliver more smoothing
  smooth.effect = predict(b, newdata = data.frame(Var1 = a$Var1, Var2 = a$Var2)) # estimate the smoothed values on our grid
  smooth.effect = array(smooth.effect, dim = dim(sens.results$results.new[[effect]]))
  colnames(smooth.effect) = colnames(sens.results$results.new[[effect]])
  rownames(smooth.effect) = rownames(sens.results$results.new[[effect]])
  
  a = melt(sens.results$results.new[[paste0("CIL.", effect)]])
  b = loess(value ~ Var1 + Var2, data = a, span = 0.8) # higher values of span will deliver more smoothing
  smooth.effect.CIL = predict(b, newdata = data.frame(Var1 = a$Var1, Var2 = a$Var2)) # estimate the smoothed values on our grid
  smooth.effect.CIL = array(smooth.effect.CIL, dim = dim(sens.results$results.new[[paste0("CIL.", effect)]]))
  colnames(smooth.effect.CIL) = colnames(sens.results$results.new[[paste0("CIL.", effect)]])
  rownames(smooth.effect.CIL) = rownames(sens.results$results.new[[paste0("CIL.", effect)]])
  
  a = melt(sens.results$results.new[[paste0("CIU.", effect)]])
  b = loess(value ~ Var1 + Var2, data = a, span = 0.8) # higher values of span will deliver more smoothing
  smooth.effect.CIU = predict(b, newdata = data.frame(Var1 = a$Var1, Var2 = a$Var2)) # estimate the smoothed values on our grid
  smooth.effect.CIU = array(smooth.effect.CIU, dim = dim(sens.results$results.new[[paste0("CIU.", effect)]]))
  colnames(smooth.effect.CIU) = colnames(sens.results$results.new[[paste0("CIU.", effect)]])
  rownames(smooth.effect.CIU) = rownames(sens.results$results.new[[paste0("CIU.", effect)]])
  
  contour(sens.results$b.m.all, sens.results$b.y.all, smooth.effect, add = T, labcex = 1)
  # Add curves that represent the contour along which the treatment effect estimate is reduced to zero
  contour(sens.results$b.m.all, sens.results$b.y.all, smooth.effect, levels = 0, add = T, col = "red", lwd = 2, lty = 2, labcex = 1)
  # Add curves for significance change
  contour(sens.results$b.m.all, sens.results$b.y.all, round(smooth.effect.CIL, 3), levels = 0, labels = "Sig.Change", add = T, col = "blue", lwd = 2, lty = 3, labcex = 1)
  contour(sens.results$b.m.all, sens.results$b.y.all, round(smooth.effect.CIU, 3), levels = 0, labels = "Sig.Change", add = T, col = "blue", lwd = 2, lty = 3, labcex = 1)
  
  message(paste(paste0(1:nrow(sens.results$X.coef.plot), " is comparable to ", rownames(sens.results$X.coef.plot), collapse = ", ")))
  if(!is.null(sens.results$b.t))
    message(paste0("This is the sensitivity plot when b.t = ", sens.results$b.t, "."))
}  
