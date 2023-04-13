phyCor=function(featureMat, trait, A, method='spearman'){
	trait = trait[rownames(featureMat)]
	x = A %*% trait
	y = A %*% featureMat
	coeff_out = cor(y, x, method=method)
	coeff_out
}



getCorrelationCoefficient=function(featureMat, trait, method='kendall'){
	trait = trait[rownames(featureMat)]
	coeff_out = cor(featureMat, trait, method=method)
	coeff_out
}



getCorrelationCoefficient_oldversion=function(featureMat, trait, method='kendall'){
	trait = trait[rownames(featureMat)]
	corrcoeff = rep(0, ncol(featureMat))
	for (i in 1:ncol(featureMat)){
		corrout = cor.test(featureMat[,i], trait, method=method)
		corrcoeff[i] = corrout$estimate
	}
	names(corrcoeff) = colnames(featureMat)
	corrcoeff
}



getCoeffsPGLS=function(featureMat, trait, A){
	# reorder trait names to match featureMat rownames
	trait = trait[rownames(featureMat)]
	mod=model.matrix(~1+trait)
	corrout = lm.glsMat.A(t(featureMat), mod, A=A, do.pval = F)
	coeff = corrout$coeff[2,]
	coeff
}

getPermulatedPhenotypesContinuous=function(trait, tree, numperms){
	nulltraits = NULL
	
	trait_reps = lapply(1:numperms, return_object, x=trait)
	nulltraits = lapply(trait_reps, simpermvec, tree=tree)
		
	names(nulltraits) = paste0("null", seq(1,numperms,1))
	nulltraits
}

return_object=function(x_idx, x){
	  return(x)
}


computePermulatedScoresApply=function(motif_idx, coeff_obs_all, coeff_null_all){
	nullscores = unlist(coeff_null_all[motif_idx,])
	nullscores = nullscores[!is.na(nullscores)]
	obs_score = coeff_obs_all[motif_idx]
	out = computePermulatedScores(obs_score, nullscores)
	out
}


computePermulatedScores=function(obs_score, nullscores){
	if (is.na(obs_score)){
		out = data.frame("permPval"=NA, "score"=NA)
	} else {
		nullscores = nullscores[!is.na(nullscores)]
		midpoint = median(nullscores)
		if (obs_score <= midpoint){
			one_sided_null_scores = nullscores[which(nullscores <= midpoint)]
			ind_extreme = which(one_sided_null_scores <= obs_score)
		} else if (obs_score > midpoint){
			one_sided_null_scores = nullscores[which(nullscores >= midpoint)]
			ind_extreme = which(one_sided_null_scores >= obs_score)
		}
		permPval = (length(ind_extreme)+1)/(length(one_sided_null_scores)+1)
		score = -log10(permPval)*sign(obs_score - midpoint)

		out = data.frame("permPval"=permPval, "score"=score)	
	}
	out
}


computePermulatedScoresV2=function(obs_score, nullscores){
	if (is.na(obs_score)){
		out = data.frame("permPval"=NA, "score"=NA, "obs_score"=NA, "relative_obs_score"=NA)
	} else {
		nullscores = nullscores[!is.na(nullscores)]
		midpoint = median(nullscores)
		relative_obs_score = obs_score-midpoint
		
		if (sign(relative_obs_score)==sign(obs_score)){
			if (obs_score <= midpoint){
				one_sided_null_scores = nullscores[which(nullscores <= midpoint)]
				ind_extreme = which(one_sided_null_scores <= obs_score)
			} else if (obs_score > midpoint){
				one_sided_null_scores = nullscores[which(nullscores >= midpoint)]
				ind_extreme = which(one_sided_null_scores >= obs_score)
			}
			permPval = (length(ind_extreme)+1)/(length(one_sided_null_scores)+1)
			score = -log10(permPval)*sign(obs_score - midpoint)
		} else {
			permPval = NA
			score = NA
		}
		out = data.frame("permPval"=permPval, "score"=score, "obs_score"=obs_score, "relative_obs_score"=relative_obs_score)
	}
	out
}

computePermulatedScoresApplyV2=function(motif_idx, coeff_obs_all, coeff_null_all){
	nullscores = unlist(coeff_null_all[motif_idx,])
	nullscores = nullscores[!is.na(nullscores)]
	obs_score = coeff_obs_all[motif_idx]
	out = computePermulatedScoresV2(obs_score, nullscores)
	out
}



